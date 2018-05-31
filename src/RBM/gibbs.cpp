#include "solver.h"
#include "gibbs.h"

using namespace arma;

Gibbs::Gibbs(double s_omega,
             double s_rho,
             int s_mc,
             int s_N,
             int s_dim,
             double s_dt,
             double sig,
             int s_H,
             bool s_interact,
             double s_spread)
:
    Solver(s_omega, s_rho, s_mc, s_N, s_dim, s_dt, sig, s_H, s_interact, s_spread)
{M = s_N*s_dim;}

double Gibbs::sample_gibbs(const vec &a, const vec &b, const mat &w,std::ofstream &myfile, ofstream &myfile2){
    double sumE = 0;
    double sumEsq = 0;
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
/*
    vec sum_d_wf(M); vec sum_d_wf_E(M);
    vec sum_d_wf_b(H); vec sum_d_wf_E_b(H);
    vec sum_d_wf_a(M); vec sum_d_wf_E_a(M);
    mat sum_d_wf_w(M,H); mat sum_d_wf_E_w(M,H);
*/

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
            X(i) = random_mu_std(mu + a(i));

        }

        for(j = 0; j < H; j++){
            Pj = prob(X,b(j),w.col(j));
            //cout << "lol" << endl;
            if(doubleRNG(genMT64)<=Pj){
                hid(j) = 0;
//                accept += 1;
            } else {
                hid(j) = 1;
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
            sumEsq += E_LGibbs*E_LGibbs;
            newE += E_LGibbs;
        }



//        myfile2 << scientific << E_LGibbs << endl;
    }



    vec mean_d_wf_a = sum_d_wf_a/(mc*H);
    vec mean_d_wf_E_a = sum_d_wf_E_a/(mc*H);
    vec mean_d_wf_b = sum_d_wf_b/(mc*H);
    vec mean_d_wf_E_b = sum_d_wf_E_b/(mc*H);
    mat mean_d_wf_w = sum_d_wf_w/(mc*H);
    mat mean_d_wf_E_w = sum_d_wf_E_w/(mc*H);
    calcg1(mean_d_wf_a,mean_d_wf_b,mean_d_wf_w);
    calcg2(mean_d_wf_E_a,mean_d_wf_E_b,mean_d_wf_E_w);
    end=clock();

    double mean_E = sumE/(mc*H);
    double mean_E_sq = sumEsq/(mc*H);
    double var = (mean_E_sq - mean_E*mean_E)/(mc*H);


    //myfile << "# Energy" <<"        " << "Variance" << "   " << "CPU time" << endl; //sanity
    myfile << scientific << sumE/(mc*H) << " " << scientific << var << " " << scientific << ((double)end-(double)start)/CLOCKS_PER_SEC << " " << endl;
    return sumE/(mc*H);
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
    static normal_distribution<double> gaussianRNG(mu,sigma*sigma);
    return gaussianRNG(genMT64);
}

rowvec Gibbs::best_params(std::ofstream &myfile, ofstream &myfile2, double gamma, vec a, vec b, mat w, int gdc){// gamma is learning rate <3
    cout << "All about that Gibbs..." << endl;
    ofstream afile; ofstream afile2;
    double energy = energy_analytic();

//    vec b = init_b();
//    mat w = init_w();
//    vec a = init_a();

    int MHMH = M+H+M*H;

    string filename ="N" + std::to_string(N)+ "_d" + std::to_string(dim)+ "_gam" + std::to_string(gamma) + "_H" + std::to_string(H) + "_sig"+std::to_string(sigma);
    afile.open("gibbs_params_" + filename + ".dat");
    //afile2.open("gibbs_energy_" + filename + ".dat");
    myfile << "# dim" << "  N " << "  mc  " << " sigma "<< " Gibbs " << endl;
    myfile << "  " << dim << "    " << N << " " << mc << " " << sigma  << endl;
    myfile << "#" << endl;
    myfile << scientific << "# Theoretical Energy = " << energy << endl;
    myfile << "# Energy" <<"     " << "Variance"<<"     " << "CPU time" << endl;

    mat alphamat = zeros(gdc,MHMH);
    mat startalpha = mat(init_alpha(a,b,w));
    //startalpha.print();
    alphamat.row(0) = startalpha;

    double mean_EL;
    rowvec g1; rowvec g2; rowvec alphanow;
    for(int r=0;r<gdc-1;r++){
        mean_EL = sample_gibbs(a, b, w, myfile, myfile2);

        g1 = 0.5*getG1();
        g2 = 0.5*getG2();
        alphamat.row(r+1) = alphamat.row(r) - gamma*2*(g2 - mean_EL*g1);
        alphanow = alphamat.row(r+1);
        //afile2 << setprecision(12) << mean_EL << endl;

        //need to reconstruct
        int k = 0;
        for(int i=0;i<M;i++){
            for(int j=0;j<H;j++){
                b(j) = alphanow(j+M);
                w(i,j) = alphanow(M+H+k);
                k++;
            }
        }
        cout << r << endl;
    }
    afile << a << endl;
    afile << b << endl;
    afile << w << endl;
    //afile2.close();
    afile.close();

    return alphanow;
}

