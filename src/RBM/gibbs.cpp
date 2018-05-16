#include "solver.h"
#include "gibbs.h"

using namespace arma;

Gibbs::Gibbs(          double s_omega,
                       double s_rho,
                       int s_mc,
                       int s_N,
                       int s_dim,
                       double s_dt,
                       double sig,
                       int s_H,
                       double s_M)
:
    Solver(s_omega, s_rho, s_mc, s_N, s_dim, s_dt, sig, s_H, s_M)
{}

double Gibbs::sample_gibbs(const vec &a, const vec &b, const mat &w,std::ofstream &myfile, ofstream &myfile2){
    double energy = energy_analytic();
    myfile << "# dim = " << dim << ", N = " << N << ", dt = " << dt << " and mc = " << mc << endl << endl;
    myfile << scientific << "# Theoretical Energy = " << energy << endl << endl;
    double sumE = 0;
    start=clock();
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(0.,0.5);
    static uniform_real_distribution<double> doubleRNG(0,1);

    int j; double Pj; double E_LGibbs;
    //mat w = init_w();
    vec hid = init_h_bool();
    int i;
    double accept = 0;
    double newE;

    //a = init_a();
    //vec b = init_b()*0.001;
    vec X = init_X();

    vec sum_d_wf = zeros(M); vec sum_d_wf_E = zeros(M);
    vec sum_d_wf_b = zeros(H); vec sum_d_wf_E_b = zeros(H);
    vec sum_d_wf_a = zeros(M); vec sum_d_wf_E_a = zeros(M);
    mat sum_d_wf_w = zeros(M,H); mat sum_d_wf_E_w = zeros(M,H);
    vec dwfa; vec dwfb; mat dwfw;
    for(int k = 0; k < mc; k++){

        //X = init_X_gaus();
        newE = 0;
        //a = init_a();

        for(i = 0; i < M; i++){ // finding the ideal positions X(i)
            double mu = get_mu(i,hid,w);
            //cout << "lol" << endl;
            X(i) = random_mu_std(mu);

        }

        for(j = 0; j < H; j++){
            Pj = prob(X,b(j),w.col(j));
            //cout << "lol" << endl;
            if(doubleRNG(genMT64)<=Pj){
                hid(j) = 1;
                accept += 1;
            } else {
                hid(j) = 0;
            }
        }
        E_LGibbs = ELGibbs(a,b,w,X);
        dwfa = grad_ai(X,a);
        sum_d_wf_a += dwfa;
        sum_d_wf_E_a += dwfa*E_LGibbs;

        dwfb = grad_bj(b,X,w);
        sum_d_wf_b += dwfb;
        sum_d_wf_E_b += dwfb*E_LGibbs;

        dwfw = grad_wij(b,X,w);
        sum_d_wf_w += dwfw;
        sum_d_wf_E_w += dwfw*E_LGibbs;
        //cout << E_LGibbs << endl;
        sumE += E_LGibbs;
        newE += E_LGibbs;
        myfile2 << scientific << E_LGibbs << endl;
    }
    vec mean_d_wf_a = sum_d_wf_a/(M*mc);
    vec mean_d_wf_E_a = sum_d_wf_E_a/(M*mc);
    vec mean_d_wf_b = sum_d_wf_b/(M*mc);
    vec mean_d_wf_E_b = sum_d_wf_E_b/(M*mc);
    mat mean_d_wf_w = sum_d_wf_w/(M*mc);
    mat mean_d_wf_E_w = sum_d_wf_E_w/(M*mc);
    calcg1(mean_d_wf_a,mean_d_wf_b,mean_d_wf_w);
    calcg2(mean_d_wf_E_a,mean_d_wf_E_b,mean_d_wf_E_w);
    end=clock();
    myfile << "# Energy" << "     " << "Acceptance" << "   " << "CPU time" << "        " << "Solver" << endl;
    myfile << scientific << sumE/(mc) << " " << scientific << accept/(mc*M) << " " << scientific << ((double)end-(double)start)/CLOCKS_PER_SEC << "    " << 2 << "  # gibbs" << endl;
    return sumE/mc;
}

vec Gibbs::init_h_bool(){
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static uniform_int_distribution<long> boolRNG(0,1);
    vec hid = zeros(H);
    for(int j = 0; j < H; j++){
        hid(j) = boolRNG(genMT64);
    }
    return hid;
}

double Gibbs::prob(const vec& X, double bj, const mat &wj){
    return 1/(1+exp(u(bj,X,wj)));
}


double Gibbs::get_mu(int i,vec const &hid, const mat &w){
    double mu = 0;
    for(int j = 0; j < H; j++){
        mu += hid(j)*w(i,j);
    }
    return mu;
}

double Gibbs::random_mu_std(double mu){
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(mu,sigma);
    return gaussianRNG(genMT64);
}

rowvec Gibbs::best_params(std::ofstream &myfile, ofstream &myfile2, double gamma, int lol){// gamma is learning rate <3
    ofstream afile; ofstream afile2;

    vec b = init_b();
    mat w = init_w();
    vec a = init_a();

    int MHMH = M+H+M*H;
    afile.open("gibbs_best_N2_D2.dat");
    afile2.open("gibbs_converge_N2_D2.dat");
    mat alphamat = zeros(lol,MHMH);
    mat startalpha = mat(init_alpha(a,b,w));
    startalpha.print();
    alphamat.row(0) = startalpha;

    double mean_EL;
    rowvec g1; rowvec g2; rowvec alphanow;
    for(int r=0;r<lol-1;r++){
        mean_EL = sample_gibbs(a, b, w, myfile, myfile2);

        g1 = getG1();
        g2 = getG2();
        alphamat.row(r+1) = alphamat.row(r) - gamma*2*(g2 - mean_EL*g1);
        alphanow = alphamat.row(r+1);
        afile2 << setprecision(12) << mean_EL << endl;

        //need to reconstruct
        int k = 0;
        for(int i=0;i<M;i++){
            for(int j=0;j<H;j++){
                b(j) = alphanow(j+M);
                w(i,j) = alphanow(M+H+k);
                k++;
            }
        }
        cout << "_" << endl;
    }
    afile << a << endl;
    afile << b << endl;
    afile << w << endl;
    afile2.close();
    afile.close();

    return alphanow;
}

