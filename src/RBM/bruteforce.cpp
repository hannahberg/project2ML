#include "solver.h"
#include "bruteforce.h"

using namespace arma;


Bruteforce::Bruteforce(double s_omega,
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
    Solver(s_omega, s_rho, s_mc, s_N, s_dim, s_dt, sig, s_H, s_interact, s_spread){
    M = s_N*s_dim;
    //Mmc = 1.0/(M*mc);
}

void Bruteforce::go_brute(std::ofstream &myfile, ofstream &myfile2){
    vec a = init_a(spread);
    vec b = init_b(spread);
    mat w = init_w(spread);
    vec X = init_X();

    solve(a,b,w,X,myfile, myfile2);
}

double Bruteforce::solve(const vec &a, const vec &b, const mat &w,const vec &X, std::ofstream &myfile, std::ofstream &myfile2){

    start=clock();
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static uniform_real_distribution<double> doubleRNG(0,1);
    double sumE = 0;

    //double newE = 0;
    vec Xnew = X;

    vec Xflex = X;

    int i; int j;
    double accept = 0;
    double localenergy;
    double A;
    vec sum_d_wf = zeros(M); vec sum_d_wf_E = zeros(M);
    vec sum_d_wf_b = zeros(H); vec sum_d_wf_E_b = zeros(H);
    vec sum_d_wf_a = zeros(M); vec sum_d_wf_E_a = zeros(M);
    mat sum_d_wf_w = zeros(M,H); mat sum_d_wf_E_w = zeros(M,H);
    vec dwfa; vec dwfb; mat dwfw;
    double totsumE = 0;
    double totsumEsq = 0;
    for(i=0;i<mc;i++){
        sumE = 0;
        for(j=0;j<M;j++){

            Xnew(j) = Xflex(j) + (doubleRNG(genMT64) - 0.5)*rho;
            A = wavefunc(a,b,w,Xnew)/(wavefunc(a,b,w,Xflex));
            A *= A;

            if((A > 1) || (A > doubleRNG(genMT64))){// test if new position is more probable than random number between 0 and 1.
                Xflex(j) = Xnew(j); //accept new position
                accept += 1;
            } else { // reset
                Xnew(j) = Xflex(j);
            }

            localenergy = E_L(a,b,w,Xflex);
            sumE += localenergy; // calculate change in energy

            dwfa = grad_ai(Xflex,a);
            sum_d_wf_a += dwfa;
            sum_d_wf_E_a += dwfa*localenergy;

            dwfb = grad_bj(b,Xflex,w);
            sum_d_wf_b += dwfb;
            sum_d_wf_E_b += dwfb*localenergy;

            dwfw = grad_wij(b,Xflex,w);
            sum_d_wf_w += dwfw;
            sum_d_wf_E_w += dwfw*localenergy;

        }
        totsumE += sumE;
        totsumEsq += sumE*sumE;
        //myfile4 << scientific << bommelom/N << endl;
        //myfile2 << scientific << sumE/M << endl;
    }

    double mean_E = totsumE/(M*mc);
    double mean_E_sq = totsumEsq/(M*mc);
    double var = mean_E_sq - mean_E*mean_E;

    vec mean_d_wf_a = sum_d_wf_a*Mmc;
    vec mean_d_wf_E_a = sum_d_wf_E_a*Mmc;
    vec mean_d_wf_b = sum_d_wf_b*Mmc;
    vec mean_d_wf_E_b = sum_d_wf_E_b*Mmc;
    mat mean_d_wf_w = sum_d_wf_w*Mmc;
    mat mean_d_wf_E_w = sum_d_wf_E_w*Mmc;

    calcg1(mean_d_wf_a,mean_d_wf_b,mean_d_wf_w);
    calcg2(mean_d_wf_E_a,mean_d_wf_E_b,mean_d_wf_E_w);
    myfile2 << scientific << mean_E << endl;



    end=clock();
    //myfile << "# Energy" <<"        " << "Variance" << "   " << "CPU time" << "Acceptance" << endl; //sanity
    myfile << scientific << mean_E << " " << scientific << var << " " << scientific << ((double)end-(double)start)/CLOCKS_PER_SEC << " " << scientific << accept/(mc*M) << endl;
    return mean_E;
}

rowvec Bruteforce::best_params(std::ofstream &myfile, ofstream &myfile2, double gamma, vec a, vec b, mat w, vec X, int gdc){
    ofstream afile; ofstream afile2;
    double energy = energy_analytic();
  
    string filename ="N" + std::to_string(N)+ "_d" + std::to_string(dim)+ "_gam" + std::to_string(gamma) + "_H" + std::to_string(H)+"_rho"+std::to_string(rho);
    afile.open("brute_params_" + filename + ".dat");
    afile2.open("brute_energy_" + filename + ".dat");
    myfile << "# dim" << "  N " << "  mc  " << " rho "<< " Brute " << endl;
    myfile << "  " << dim << "    " << N << " " << mc << " " << rho  << endl;
    myfile << "#" << endl;
    myfile << scientific << "# Theoretical Energy = " << energy << endl;
    myfile << "# Energy" <<"     " << "Variance"<<"     " << "CPU time" << "     Acceptance" << endl;
    rowvec alpha_best;
    int MHMH = M+H +M*H;

    mat alphamat = zeros(gdc,MHMH);

    mat startalpha = mat(init_alpha(a,b,w));
    //startalpha.print();
    alphamat.row(0) = startalpha;

    double mean_EL;
    rowvec g1; rowvec g2; rowvec alphanow;
    for(int r=0;r<gdc-1;r++){
        mean_EL = solve(a, b, w, X, myfile, myfile2);

        g1 = getG1();
        g2 = getG2();
        alphamat.row(r+1) = alphamat.row(r)- gamma*2*(g2 - mean_EL*g1);
        alphanow = alphamat.row(r+1);
        //afile2 << setprecision(12) << mean_EL << endl;
        cout << mean_EL << endl;

        //need to reconstruct
        //alphanow = alphamat.row(r+1);
        int k = 0;
        for(int i=0;i<M;i++){
            for(int j=0;j<H;j++){
                b(j) = alphanow(j+M);
                w(i,j) = alphanow(M+H+k);
                k++;
            }
        }

            //afile << i << endl; //checking if it oscillates in the bottom
        }
    afile2.close();
    alpha_best = alphanow;
    afile <<  setprecision(12)  << a << endl;
    afile <<  setprecision(12) << b << endl;
    afile <<  setprecision(12) << w << endl;
    afile.close();
    return alpha_best;
}
