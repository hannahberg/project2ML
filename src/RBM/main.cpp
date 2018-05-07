#include <iostream>
#include "solver.h"
#include "bruteforce.h"
#include "impsamp.h"
//#include "interact.h"
using namespace std;
using namespace arma;

int main(){
    vec dtvec = logspace<vec>(-4,0,51);
    double rho = 0.1;
    //double dt = 0.01; // [0.001, 0.01]
    double h = 0.0001;
    int numpart = 1; //CHANGE THE NAME!!!!!!!!!!!!!!!!!!!!!!!!!
    int mc = (1048576 + 1000) / numpart; // monte carlo cycles
    int howmanyDs = 1;
    double hbar = 1;
    double mass = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    double hidden = 2;
    double dt = 0.001;

    ofstream myfile;
    //myfile.open("interaction_N10.dat");
    //for(int elem=0; elem<size(dtvec,0); elem++){
    //    if(alphavec(elem) != 0.5){
    //dtvec(elem);
    Solver S(hbar, mass, omega, rho, mc, numpart, howmanyDs, h, dt, sig, hidden, numM); // initialize Solver class
    Bruteforce* B = new Bruteforce(hbar, mass, omega, rho, mc, numpart, howmanyDs, h, dt, sig, hidden, numM);
    Impsamp* Imp = new Impsamp(hbar, mass, omega, rho, mc, numpart, howmanyDs, h, dt, sig, hidden, numM);
    //Interact* Int = new Interact(hbar, mass, omega, rho, mc, numpart, howmanyDs, h, dt, sig, hidden, numM);

    ofstream myfile2;
    ofstream myfile3;
    ofstream myfile4;
    ofstream myfile5;
    string filename ="test_N" + std::to_string(numpart)+ "_d" + std::to_string(howmanyDs);
    myfile.open(filename + ".dat");
    myfile2.open(filename + "_energy.dat");
    //B->solve(myfile, myfile2);
    Imp->langevin(myfile,myfile2);
    //delete Int;
    delete Imp;
    delete B;
    //}
    myfile.close();
    myfile2.close();
    //}
    //}
}
