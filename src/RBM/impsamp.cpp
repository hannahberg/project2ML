#include "solver.h"
#include "impsamp.h"

using namespace arma;

Impsamp::Impsamp(double s_omega,
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


void Impsamp::go_imp(std::ofstream &myfile, ofstream &myfile2){
    vec a = init_a(spread);
    vec b = init_b(spread);
    mat w = init_w(spread);
    vec X = init_X();

    langevin(a,b,w,X,myfile, myfile2);
}

double Impsamp::langevin(const vec &a, const vec &b, const mat &w,const vec &Xin,std::ofstream &myfile, ofstream &myfile2){

    start=clock();
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(0.,0.5);
    static uniform_real_distribution<double> doubleRNG(0,1);
    double sdt = sqrt(dt);
    double newE = 0; double sumE;
    double totsumE = 0;
    double totsumEsq = 0;
    vec X = Xin;
    vec Xnew = X;
    int i; int j;
    double accept = 0;
    double greens;
    vec F = drift(b,X,w,a);
    vec Fnew = F;
    double D = 0.5; // diffusion coefficient
    double Ddt = D*dt;
    double Ddt05 = 0.5*D*dt;
    double sigma_2 = 2/(sigma*sigma);
    double bajs; double A;
    vec sum_d_wf = zeros(M); vec sum_d_wf_E = zeros(M);
    vec sum_d_wf_b = zeros(H); vec sum_d_wf_E_b = zeros(H);
    vec sum_d_wf_a = zeros(M); vec sum_d_wf_E_a = zeros(M);
    mat sum_d_wf_w = zeros(M,H); mat sum_d_wf_E_w = zeros(M,H);
    vec dwfa; vec dwfb; mat dwfw;
    for(i = 0; i < mc; i++){
        sumE = 0;

        for(j = 0; j < M; j++){
            Xnew(j) = X(j) + Ddt*F(j) + gaussianRNG(genMT64)*sdt;
            Fnew(j) = sigma_2*(-X(j)+a(j) + drifti(b,X,w,j));

            A = wavefunc(a,b,w,Xnew)/(wavefunc(a,b,w,X));
            greens = 0.5*(F(j)+Fnew(j))*(Ddt05*(F(j)-Fnew(j))-Xnew(j)+X(j));
            greens = exp(greens);
            A *= A;
            A = A*greens;
            // test if new position is more probable than random number between 0 and 1.
            if((A > 1) || (A > doubleRNG(genMT64))){
                X(j) = Xnew(j); //accept new position
                F(j) = Fnew(j);
                accept += 1;
            } else {
                Xnew(j) = X(j);
                Fnew(j) = F(j);
            }


            bajs = E_L(a,b,w,X); // local energy
            newE += bajs; // calculate change in energy
            sumE += bajs;

            dwfa = grad_ai(X,a);
            sum_d_wf_a += dwfa;
            sum_d_wf_E_a += dwfa*bajs;

            dwfb = grad_bj(b,X,w);
            sum_d_wf_b += dwfb;
            sum_d_wf_E_b += dwfb*bajs;

            dwfw = grad_wij(b,X,w);
            sum_d_wf_w += dwfw;
            sum_d_wf_E_w += dwfw*bajs;
       }

    myfile2 << scientific << sumE/M << endl;
    totsumE += sumE;
    totsumEsq += sumE*sumE;
    }
    double mean_E = totsumE/(M*mc);
    double mean_E_sq = totsumEsq/(M*mc);
    double var = mean_E_sq - mean_E*mean_E;

    //cout << "Impsamp finished! The end is near <3" << endl;

    /*
    double energy = energySum/(mc * N);
    double totalenergy = energySum/mc;
    double energySquared = energySquaredSum/(mc * N);
    */
    double E_ = newE/(mc*M);
    vec mean_d_wf_a = sum_d_wf_a/(M*mc);
    vec mean_d_wf_E_a = sum_d_wf_E_a/(M*mc);
    vec mean_d_wf_b = sum_d_wf_b/(M*mc);
    vec mean_d_wf_E_b = sum_d_wf_E_b/(M*mc);
    mat mean_d_wf_w = sum_d_wf_w/(M*mc);
    mat mean_d_wf_E_w = sum_d_wf_E_w/(M*mc);
    calcg1(mean_d_wf_a,mean_d_wf_b,mean_d_wf_w);
    calcg2(mean_d_wf_E_a,mean_d_wf_E_b,mean_d_wf_E_w);
    end=clock();
    cout << E_ << endl;
    //myfile << "# Energy" <<"        " << "Variance" << "   " << "CPU time" << "Acceptance" << endl; //sanity
    myfile << scientific << E_ << " " << scientific << var << " " << scientific << ((double)end-(double)start)/CLOCKS_PER_SEC << " " << scientific << accept/(mc*M) << endl;
    return E_;
}

rowvec Impsamp::best_params(std::ofstream &myfile, ofstream &myfile2, double gamma, vec a, vec b, mat w, vec X, int gdc){
    // gamma is learning rate <3
    ofstream afile; ofstream afile2;
    double energy = energy_analytic();

//    vec b = init_b();
//    mat w = init_w();
//    vec X = init_X();
//    vec a = init_a();
    int MHMH = M+H+M*H;

    string filename ="_N" + std::to_string(N)+ "_d" + std::to_string(dim)+ "_gam" + std::to_string(gamma) + "_H" + std::to_string(H)+"_dt"+std::to_string(dt);
    afile.open("imp_params" + filename + ".dat");
    //afile2.open("imp_energy" + filename + ".dat");

    myfile << "# dim" << "  N " << "  mc  " << " dt "<< " Impsamp " << endl;
    myfile << "  " << dim << "    " << N << " " << mc << " " << dt  << endl;
    myfile << "#" << endl;
    myfile << scientific << "# Theoretical Energy = " << energy << endl;
    myfile << "# Energy" <<"     " << "Variance"<<"     " << "CPU time" << "     Acceptance" << endl;


    mat alphamat = zeros(gdc,MHMH);
    mat startalpha = mat(init_alpha(a,b,w));
//    startalpha.print();
    alphamat.row(0) = startalpha;

    double mean_EL;
    rowvec g1; rowvec g2; rowvec alphanow;
    for(int r=0;r<gdc-1;r++){
        mean_EL = langevin(a, b, w, X, myfile, myfile2);

        g1 = getG1();
        g2 = getG2();
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
    afile.close();
    //afile2.close();
    return alphanow;
}

