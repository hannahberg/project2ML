#include <iostream>
#include "solver.h"
#include "bruteforce.h"
#include "impsamp.h"
#include "interact.h"
using namespace std;
using namespace arma;

int main(){
    vec alphavec = linspace<vec>(0.3, 0.7, 5);
    vec dtvec = logspace<vec>(-4,0,51);
    double rho = 0.1;
    //double dt = 0.01; // [0.001, 0.01]
    double h = 0.0001;
    int numpart = 100; //CHANGE THE NAME!!!!!!!!!!!!!!!!!!!!!!!!!
    int mc = (1048576 + 1000) / numpart; // monte carlo cycles
    //int mc = 100000 / numpart;
    int howmanyDs = 3;
    double beta = 1;
    double hbar = 1;
    double mass = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    double hidden = 3;

    ofstream myfile;
    //myfile.open("interaction_N10.dat");
    //for(int elem=0; elem<size(dtvec,0); elem++){
    //    if(alphavec(elem) != 0.5){
        double alpha = 0.5;
    double dt = 0.01;//dtvec(elem);
    Solver S(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt, sig, numM, hidden); // initialize Solver class
    Bruteforce* B = new Bruteforce(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
    Impsamp* Imp = new Impsamp(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
    Interact* Int = new Interact(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);

    ofstream myfile2;
    ofstream myfile3;
    ofstream myfile4;
    ofstream myfile5;
    string filename ="e_a0_b1_n" + std::to_string(numpart)+ "_d" + std::to_string(howmanyDs);
    myfile.open(filename + ".dat");     //CHANGE MY NAME!!!!!!!!!!!!!  DONT YOU DARE NOT CHANGE ME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //myfile2.open("Ea_" + filename + ".dat");
    //myfile3.open("En_" + filename + ".dat");
    //myfile4.open("Eimp_" + filename + "dt" + std::to_string(dt) + ".dat");
    //myfile4.open(std::to_string(dt) + ".dat");
    myfile5.open("Eint_" + filename  + "_alpha"  + std::to_string(alpha) + ".dat");
    //B->solve(myfile,myfile2);
    //B->solve_num(myfile,myfile3);
    //Imp->langevin(myfile,myfile4,alpha);
    Int->solve_interact(myfile, myfile5, alpha);
    myfile.close();
    //myfile2.close();
    //myfile3.close();
    //myfile4.close();
    myfile5.close();
    //Int->best_alpha();
    delete Int;
    delete Imp;
    delete B;
    //}
    myfile.close();
    //}
    //}
}
