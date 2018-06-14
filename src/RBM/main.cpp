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

    //vec gammavec = {0.01, 0.05, 0.1, 0.5, 1.};
    //vec rhovec = {0.1,0.3,0.6,0.8,1.0};
    //vec dtvec = {0.001,0.005,0.01,0.1,1};
    //vec hidvec = {0,1,2,3,5,10};
    //vec numvec = {1,2};
    //vec sigmavec = {0.6,0.7,0.8};



    for(int i=0;i<1;i++){
        //rho = rhovec(i);
        //gamma = gammavec(i);
        //hidden = hidvec(i);
        //sig = sigmavec(i);
        //numpart = numvec(i);

        Solver S(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread); // initialize Solver class
        string filename ="N" + std::to_string(numpart)+ "_d" + std::to_string(howmanyDs)+"gam" + std::to_string(gamma) + "_H" + std::to_string(hidden);
        vec X = S.init_X(); //GAUSSIAN IN IMPSAMP?
        vec a = S.init_a(spread);
        vec b = S.init_b(spread);
        mat w = S.init_w(spread);

        gamma = 0.01;

        filename ="N" + std::to_string(numpart)+ "_d" + std::to_string(howmanyDs)+"gam" + std::to_string(gamma) + "_H" + std::to_string(hidden);
        ofstream brutefile;
        ofstream brutefile2;
        brutefile.open("brute_" + filename + "_rho"+std::to_string(rho) +  ".dat");
        Bruteforce* B = new Bruteforce(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
        B->best_params(brutefile,brutefile2,gamma,a,b,w,X,gdc);
        brutefile.close();

        gamma = 0.01;
        filename ="N" + std::to_string(numpart)+ "_d" + std::to_string(howmanyDs)+"gam" + std::to_string(gamma) + "_H" + std::to_string(hidden);
        ofstream impfile;
        ofstream impfile2;
        impfile.open("imp_" + filename + "_dt"+std::to_string(dt) + ".dat");
        Impsamp* I = new Impsamp(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
        I->best_params(impfile,impfile2,gamma,a,b,w,X,gdc);
        impfile.close();

        gamma = 0.01;
        sig = 0.7;
        filename ="N" + std::to_string(numpart)+ "_d" + std::to_string(howmanyDs)+"gam" + std::to_string(gamma) + "_H" + std::to_string(hidden);
        Solver S2(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);

        //b = S2.init_b(spread);
        //w = S2.init_w(spread);
        ofstream gibfile;
        ofstream gibfile2;
        gibfile.open("gibbs_" + filename + "_sig"+std::to_string(sig) + ".dat");
        Gibbs* G = new Gibbs(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
        G->best_params(gibfile,gibfile2,gamma,a,b,w,gdc);
        gibfile.close();


        delete G;
        delete I;
        delete B;
/*

    }
    hidden = 1;
    for(int k=0;k<5;k++){
        gamma = gammavec(k);

        Solver S3(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread); // initialize Solver class
        string filename ="N" + std::to_string(numpart)+ "_d" + std::to_string(howmanyDs)+"gam" + std::to_string(gamma) + "_H" + std::to_string(hidden);
        ofstream impfile3;
        ofstream impfile4;
        impfile3.open("imp_" + filename + "_dt"+std::to_string(dt) + ".dat");
        Impsamp* I = new Impsamp(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
        I->best_params(impfile3,impfile4,gamma,a,b,w,X,gdc);
        impfile3.close();

        delete I;
    }


*/
//    double exact_vint = 1./sqrt(2) + 1./sqrt(18)+1./sqrt(8);
//    vec Xtry = {0, 0, 1, 1, 3, 3, 7,7};
//    cout << scientific << exact_vint << endl;
//    cout << scientific << S.calc_interaction(Xtry) << endl;


    //}
    //}
    // thermalize the grad
    // vary spread
    // vary learning rate
    // vary hidden nodes, in non-interaction case H = 0 is best
    // interaction - 5 nodes?
}
}
