#include <iostream>
#include "solver.h"
#include "bruteforce.h"
#include "impsamp.h"
#include "gibbs.h"
//#include "interact.h"
using namespace std;
using namespace arma;

int main(){
    vec dtvec = logspace<vec>(-4,0,51);
    double rho = 0.1;
    int numpart = 1; //CHANGE THE NAME!!!!!!!!!!!!!!!!!!!!!!!!!
    int mc = 100000;//(1048576 + 1000) / numpart; // monte carlo cycles
    int howmanyDs = 1;
    double hbar = 1;
    double mass = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.001;
    vec gammavec = {0.5, 1, 2, 3, 4};
    double gamma;

    ofstream myfile;
    //myfile.open("interaction_N10.dat");
    //for(int elem=0; elem<size(dtvec,0); elem++){
    //    if(alphavec(elem) != 0.5){
    //dtvec(elem);
    Solver S(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM); // initialize Solver class
    Bruteforce* B = new Bruteforce(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM);
    Impsamp* Imp = new Impsamp(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM);
    //Interact* Int = new Interact(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM);
    Gibbs* G = new Gibbs(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM);
    ofstream myfile2;
    ofstream myfile3;
    ofstream myfile4;
    ofstream myfile5;
    string filename ="test_N" + std::to_string(numpart)+ "_d" + std::to_string(howmanyDs);
    myfile.open(filename + ".dat");
    myfile2.open(filename + "_energy.dat");

    //B->solve();
    vec b = S.init_b();
    mat w = S.init_w();
    vec X = S.init_X();
    vec a = S.init_a();
    for(int i=0;i<5;i++){
        gamma = gammavec(i);
        B->best_params(myfile,myfile2,gamma,a,b,w,X);
    }
    //B->go_brute(myfile, myfile2);

    //B->solve(myfile, myfile2);
    //Imp->langevin(myfile,myfile2);

    G->sample_gibbs(myfile,myfile2);


    //delete Int;
    delete G;
    delete Imp;
    delete B;
    //}
    myfile.close();
    myfile2.close();
    //}
    //}
}
