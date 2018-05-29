#include <iostream>
#include "solver.h"
#include "bruteforce.h"
#include "impsamp.h"
#include "gibbs.h"
//#include "interact.h"
#include <limits>
using namespace std;
using namespace arma;

int main(){
    double rho = 0.4;
    int numpart = 2;
    int mc = 100000;//(1048576 + 1000) / numpart; // Monte Carlo Cycles
    int howmanyDs = 2;
    double omega = 1;
    double sig = 1;
    int hidden = 4;
    double dt = 0.01;
    double gamma = 0.1;
    int gdc = 100; // Gradient Decent Cycles
    bool interactionswitch = true;

    //vec gammavec = {0.5, 1, 2, 3, 4};
    //vec dtvec = logspace<vec>(-4,0,51);


    for(int i=0;i<1;i++){
        Solver S(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch); // initialize Solver class
        string filename ="_N" + std::to_string(N)+ "_d" + std::to_string(dim)+ "gam" + std::to_string(gamma) + "_H" + std::to_string(H);

        vec b = S.init_b();
        mat w = S.init_w();
        vec X = S.init_X(); //GAUSSIAN IN IMPSAMP?
        vec a = S.init_a();

        ofstream brutefile;
        ofstream brutefile2;
        brutefile.open("brute_" + filename + ".dat");
        brutefile2.open("brute_" + filename + "_energy.dat");
        Bruteforce* B = new Bruteforce(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch);
        B->best_params(brutefile,brutefile2,gamma,a,b,w,X,gdc);
        brutefile2.close();
        brutefile.close();

        ofstream impfile;
        ofstream impfile2;
        impfile.open("imp_" + filename + "_dt"+std::to_string(dt) + ".dat");
        impfile2.open("imp_" + filename + "_dt"+std::to_string(dt) + "_energy.dat");
        Impsamp* I = new Impsamp(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch);
        I->best_params(impfile,impfile2,gamma,a,b,w,X,gdc);
        impfile2.close();
        impfile.close();

        ofstream gibfile;
        ofstream gibfile2;
        gibfile.open("gibbs_" + filename + ".dat");
        gibfile2.open("gibbs_" + filename + "_energy.dat");
        Gibbs* G = new Gibbs(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch);
        G->best_params(gibfile,gibfile2,gamma,a,b,w,gdc);
        gibfile2.close();
        gibfile.close();

//       B->solve(a, b, w, X, myfile, myfile2);

        delete G;
        delete I;
        delete B;
        delete S;

        /*
        vec b = {0.000667904,0.00171753};
        mat w = {{0.00117447, 0.000381313},{-0.00166561, -0.000559385},{-0.000244628, -0.00130079}, {0.00103996, -0.00180347}};
        vec X = {1,2,3,4};
        vec X = {-0.436988, -0.243587,0.22773,0.146376};
        vec a ={0.000586069,0.000838821,0.000520396,0.000574091};
        */
    }



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
