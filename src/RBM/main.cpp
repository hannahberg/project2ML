#include <iostream>
#include "solver.h"
#include "bruteforce.h"
#include "impsamp.h"
#include "gibbs.h"
#include <limits>
using namespace std;
using namespace arma;


int main(){
    double rho = 1.0;
    int numpart = 2;
    int howmanyDs = 2;
    int mc = (1048576 + 1000) / (numpart*howmanyDs); // Monte Carlo Cycles
    double omega = 1;
    double sig = 1;
    int hidden = 2;
    double dt = 0.1;
    double gamma = 0.1;
    int gdc = 5000; // Gradient Decent Cycles
    bool interactionswitch = true;
    double spread = 0.2;

    for(int i=0;i<1;i++){

        Solver S(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread); // initialize Solver class
        string filename ="N" + std::to_string(numpart)+ "_d" + std::to_string(howmanyDs)+"gam" + std::to_string(gamma) + "_H" + std::to_string(hidden);
        vec X = S.init_X(); 
        vec a = S.init_a(spread);
        vec b = S.init_b(spread);
        mat w = S.init_w(spread);

        filename ="N" + std::to_string(numpart)+ "_d" + std::to_string(howmanyDs)+"gam" + std::to_string(gamma) + "_H" + std::to_string(hidden);
        ofstream brutefile;
        ofstream brutefile2;
        brutefile.open("brute_" + filename + "_rho"+std::to_string(rho) +  ".dat");
        Bruteforce* B = new Bruteforce(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
        B->best_params(brutefile,brutefile2,gamma,a,b,w,X,gdc);
        brutefile.close();

        filename ="N" + std::to_string(numpart)+ "_d" + std::to_string(howmanyDs)+"gam" + std::to_string(gamma) + "_H" + std::to_string(hidden);
        ofstream impfile;
        ofstream impfile2;
        impfile.open("imp_" + filename + "_dt"+std::to_string(dt) + ".dat");
        Impsamp* I = new Impsamp(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
        I->best_params(impfile,impfile2,gamma,a,b,w,X,gdc);
        impfile.close();

        sig = 0.7;
        filename ="N" + std::to_string(numpart)+ "_d" + std::to_string(howmanyDs)+"gam" + std::to_string(gamma) + "_H" + std::to_string(hidden);
        Solver S2(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);

        ofstream gibfile;
        ofstream gibfile2;
        gibfile.open("gibbs_" + filename + "_sig"+std::to_string(sig) + ".dat");
        Gibbs* G = new Gibbs(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
        G->best_params(gibfile,gibfile2,gamma,a,b,w,gdc);
        gibfile.close();


        delete G;
        delete I;
        delete B;

}
}
