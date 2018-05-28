#include "solver.h"
#include "bruteforce.h"

using namespace arma;

Bruteforce::Bruteforce(double s_hbar,
                       double mass,
                       double s_omega,
                       double s_rho,
                       int s_mc,
                       int s_N,
                       int s_dim,
                       double s_dt,
                       double sig,
                       int s_H,
                       double s_M,
                       bool s_interact)
:
    Solver(s_hbar, mass,s_omega, s_rho, s_mc, s_N, s_dim, s_dt, sig, s_H, s_M,s_interact)
{}

double Bruteforce::energy_analytic(){
    return 0.5 * N * dim;
}

void Bruteforce::go_brute(std::ofstream &myfile, ofstream &myfile2){
    vec a = init_a();
    vec b = init_b();
    mat w = init_w();
    vec X = init_X();

    solve(a,b,w,X,myfile, myfile2);
}

double Bruteforce::solve(const vec &a, const vec &b, const mat &w,const vec &X, std::ofstream &myfile, std::ofstream &myfile2){
    double energy = energy_analytic();
    myfile << "# dim = " << dim << ", N = " << N << ", dt = " << dt << " and mc = " << mc << endl << endl;
    myfile << scientific << "# Theoretical Energy = " << energy << endl << endl;

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
    for(i=0;i<mc;i++){
        //double bajsen = 0;
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
        //myfile4 << scientific << bommelom/N << endl;
        //myfile2 << scientific << sumE/M << endl;
    }

    double mean_E = totsumE/(M*mc);
    vec mean_d_wf_a = sum_d_wf_a/(M*mc);
    vec mean_d_wf_E_a = sum_d_wf_E_a/(M*mc);
    vec mean_d_wf_b = sum_d_wf_b/(M*mc);
    vec mean_d_wf_E_b = sum_d_wf_E_b/(M*mc);
    mat mean_d_wf_w = sum_d_wf_w/(M*mc);
    mat mean_d_wf_E_w = sum_d_wf_E_w/(M*mc);

    calcg1(mean_d_wf_a,mean_d_wf_b,mean_d_wf_w);
    calcg2(mean_d_wf_E_a,mean_d_wf_E_b,mean_d_wf_E_w);
    myfile2 << scientific << mean_E << endl;

    //cout << "Brute force finished! Hang in there <3" << endl;

    /*
    double energy = energySum/(mc * N);
    double totalenergy = energySum/mc;
    double energySquared = energySquaredSum/(mc * N);
    */

    end=clock();
    myfile << "# Energy" << "     " << "Acceptance" << "   " << "CPU time" << "        " << "Solver" << endl;
    myfile << scientific << mean_E << " " << scientific << accept/(mc*M) << " " << scientific << ((double)end-(double)start)/CLOCKS_PER_SEC << "    " << 0 << "  # Analytic" << endl;
    return mean_E;
}

rowvec Bruteforce::best_params(std::ofstream &myfile, ofstream &myfile2, double gamma, vec a, vec b, mat w, vec X){
    ofstream afile; ofstream afile2;

    string filename ="gam" + std::to_string(gamma) + "_N" + std::to_string(N)+ "_d" + std::to_string(dim)+ "_H" + std::to_string(H)+"_rho"+std::to_string(rho);
    afile.open("params_" + filename + ".dat");
    afile2.open("energy_" + filename + ".dat");
    rowvec alpha_best;
    int lol = 100;
    int MHMH = M+H +M*H;

    mat alphamat = zeros(lol,MHMH);

    mat startalpha = mat(init_alpha(a,b,w));
    //startalpha.print();
    alphamat.row(0) = startalpha;

    double mean_EL;
    rowvec g1; rowvec g2; rowvec alphanow;
    for(int r=0;r<lol-1;r++){
        mean_EL = solve(a, b, w, X, myfile, myfile2);

        g1 = getG1();
        g2 = getG2();
        alphamat.row(r+1) = alphamat.row(r)- gamma*2*(g2 - mean_EL*g1);
        alphanow = alphamat.row(r+1);
        afile2 << setprecision(12) << mean_EL << endl;
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
