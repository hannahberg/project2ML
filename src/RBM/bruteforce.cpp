#include "solver.h"
#include "bruteforce.h"

using namespace arma;

Bruteforce::Bruteforce(
                       double s_hbar,
                       double mass,
                       double s_omega,
                       double s_rho,
                       int s_mc,
                       int s_N,
                       int s_dim,
                       double s_h,
                       double s_dt,
                       double sig,
                       double s_H,
                       double s_M)
:
    Solver(s_hbar, mass,s_omega, s_rho, s_mc, s_N, s_dim, s_h, s_dt, sig, s_H, s_M)
{}

double Bruteforce::energy_analytic(){
    return 0.5 * N * dim;
}

void Bruteforce::solve(std::ofstream &myfile, ofstream &myfile2){
    double energy = energy_analytic();
    myfile << "# dim = " << dim << ", N = " << N << ", dt = " << dt << " and mc = " << mc << endl << endl;
    myfile << scientific << "# Theoretical Energy = " << energy << endl << endl;

    start=clock();
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static uniform_real_distribution<double> doubleRNG(0,1);

    double newE = 0;
    vec a = init_a();
    vec b = init_b();
    mat w = init_w();
    vec X = init_X();
    vec Xnew = X;
    int i; int j; int q;
    double accept = 0;

    for(i=0;i<mc;i++){
        double bajsen = 0;
        for(j=0;j<N;j++){
            //for(q=0;q<dim;q++){
                Xnew(j) = X(j) + (doubleRNG(genMT64) - 0.5)*rho;
           // }

            double A = (wavefunc(a,b,w,X))/wavefunc(a,b,w,Xnew);
            A *= A;

            // test if new position is more probable than random number between 0 and 1.
            if((A > 1) || (A > doubleRNG(genMT64))){
                X(j) = Xnew(j); //accept new position
                accept += 1;
            } else {
                Xnew(j) = X(j);
            }
            double bajs = E_L(a,b,w,X);
            newE += bajs; // calculate change in energy
            bajsen += bajs;
       }

    myfile2 << scientific << bajsen/N << endl;
    }

    cout << "Brute force finished! Hang in there <3" << endl;

    /*
    double energy = energySum/(mc * N);
    double totalenergy = energySum/mc;
    double energySquared = energySquaredSum/(mc * N);
    */

    end=clock();
    myfile << "# Energy" << "     " << "Acceptance" << "   " << "CPU time" << "        " << "Solver" << endl;
    myfile << scientific << newE/(mc*N) << " " << scientific << accept/(mc*N) << " " << scientific << ((double)end-(double)start)/CLOCKS_PER_SEC << "    " << 0 << "  # Analytic" << endl;
}
