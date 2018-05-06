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

void Bruteforce::go_brute(std::ofstream &myfile, ofstream &myfile2){
    vec a = init_a();
    vec b = init_b();
    mat w = init_w();
    vec X = init_X();

    solve(a,b,w,X,myfile, myfile2);
}

vec Bruteforce::solve(const vec &a, const vec &b, const mat &w,const vec &X, std::ofstream &myfile, std::ofstream &myfile2){
    double energy = energy_analytic();
    myfile << "# dim = " << dim << ", N = " << N << ", dt = " << dt << " and mc = " << mc << endl << endl;
    myfile << scientific << "# Theoretical Energy = " << energy << endl << endl;

    start=clock();
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static uniform_real_distribution<double> doubleRNG(0,1);
    double sumE = 0;
    double sum_d_wf = 0;
    double sum_d_wf_E = 0;
    double newE = 0;
    vec Xnew = X;
    vec Xflex = X;
    int i; int j; int q;
    double accept = 0;

    for(i=0;i<mc;i++){
        double bajsen = 0;
        for(j=0;j<M;j++){
            //for(q=0;q<dim;q++){
                Xnew(j) = Xflex(j) + (doubleRNG(genMT64) - 0.5)*rho;
           // }

            double A = (wavefunc(a,b,w,Xflex))/wavefunc(a,b,w,Xnew);

            A *= A;

            // test if new position is more probable than random number between 0 and 1.
            if((A > 1) || (A > doubleRNG(genMT64))){
                Xflex(j) = Xnew(j); //accept new position
                accept += 1;
            } else {
                Xnew(j) = Xflex(j);
            }


            double localenergy = E_L(a,b,w,Xflex);
            sumE += localenergy; // calculate change in energy
            double dwf = grad_ai(Xflex,a);
            sum_d_wf += dwf;
            sum_d_wf_E += dwf*localenergy;
            }
    //
        //myfile4 << scientific << bommelom/N << endl;
    }
    double mean_E = sumE/(M*mc);
    double mean_d_wf = sum_d_wf/(M*mc);
    double mean_d_wf_E = sum_d_wf_E/(M*mc);

    vec mean_values = zeros(3);
    mean_values(0) = mean_E;
    mean_values(1) = mean_d_wf;
    mean_values(2) = mean_d_wf_E;

    myfile2 << scientific << sumE/M << endl;

    cout << "Brute force finished! Hang in there <3" << endl;

    /*
    double energy = energySum/(mc * N);
    double totalenergy = energySum/mc;
    double energySquared = energySquaredSum/(mc * N);
    */

    end=clock();
    myfile << "# Energy" << "     " << "Acceptance" << "   " << "CPU time" << "        " << "Solver" << endl;
    myfile << scientific << sumE/(mc*M) << " " << scientific << accept/(mc*M) << " " << scientific << ((double)end-(double)start)/CLOCKS_PER_SEC << "    " << 0 << "  # Analytic" << endl;
    return mean_values;
}

vec Bruteforce::best_a(std::ofstream &myfile, ofstream &myfile2){
    ofstream afile;ofstream afile2;

    vec b = init_b();
    mat w = init_w();
    vec X = init_X();

    afile.open("a_start.dat");
    afile2.open("a_converge.dat");
    vec a_best;
    vec mean_values;
    int lol = 1000;
    mat a_ = zeros(lol,M);
    a_.row(0) = init_a();
    double gamma = 0.2; //hilsen alocias 0.02
    double tol = 0.01;
    for(int i=0;i<lol-1;i++){
        mean_values = solve(a_.row(i), b, w, X, myfile, myfile2);
        double mean_EK = mean_values(0);
        double mean_d_wf = mean_values(1);
        double mean_d_wf_E = mean_values(2);
        a_.row(i+1) = a_.row(i)- gamma*2*(mean_d_wf_E - mean_EK*mean_d_wf);
        afile2 << setprecision(12) << a_(i) << endl;

            //afile << i << endl; //checking if it oscillates in the bottom
        }
    a_best = a_(lol);
    afile.close();
    afile2.close();
    return a_best;
}
