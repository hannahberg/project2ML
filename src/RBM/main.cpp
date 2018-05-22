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
    //double rho = 0.01;
    int numpart = 2; //CHANGE THE NAME!!!!!!!!!!!!!!!!!!!!!!!!!
    int mc = 100000;//(1048576 + 1000) / numpart; // monte carlo cycles
    int howmanyDs = 2;
    double hbar = 1;
    double mass = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    vec gammavec = {0.5, 1, 2, 3, 4};
    double gamma;
    bool interactionswitch = true;
    ofstream myfile;
    //myfile.open("interaction_N10.dat");
    //for(int elem=0; elem<size(dtvec,0); elem++){
    //    if(alphavec(elem) != 0.5){
    //dtvec(elem);
    //Solver S(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM,interactionswitch); // initialize Solver class
    //Bruteforce* B = new Bruteforce(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM,interactionswitch);
    //Impsamp* Imp = new Impsamp(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM,interactionswitch);
    //Interact* Int = new Interact(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM);
    //Gibbs* G = new Gibbs(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM,interactionswitch);
    ofstream myfile2;
    ofstream myfile3;
    ofstream myfile4;
    ofstream myfile5;
    //string filename ="test_N" + std::to_string(numpart)+ "_d" + std::to_string(howmanyDs);
//    myfile.open(filename + ".dat");
//    myfile2.open(filename + "_energy.dat");

    //B->solve();

    int grad_cycle = 1000;
    double learning = 0.2;

    //Imp->best_params(myfile,myfile2,learning,grad_cycle);
    //G->best_params(myfile,myfile2,learning,grad_cycle);
    /*
    vec b = S.init_b();
    mat w = S.init_w();
    vec X = S.init_X();
    vec a = S.init_a();

    vec b = {0.000667904,0.00171753};
    mat w = {{0.00117447, 0.000381313},{-0.00166561, -0.000559385},{-0.000244628, -0.00130079}, {0.00103996, -0.00180347}};
//    vec X = {1,2,3,4};
    vec X = {-0.436988, -0.243587,0.22773,0.146376};
    vec a ={0.000586069,0.000838821,0.000520396,0.000574091};
    */
    /*cout <<"b" << endl;
    b.raw_print();
    cout <<"w" << endl;
    w.raw_print();
    cout <<"X" << endl;
    X.raw_print();
    cout <<"a" << endl;
    a.raw_print();*/
//    B->solve(a, b, w, X,myfile, myfile2);

    for(int i=0;i<1;i++){
        double rho = 2.5;

        string filename ="newintdy_N" + std::to_string(numpart)+ "_d" + std::to_string(howmanyDs) + "_dt"+std::to_string(dt);
        myfile.open(filename + ".dat");
        myfile2.open(filename + "_energy.dat");
        gamma = 0.2;
        Solver S(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM,interactionswitch); // initialize Solver class
        vec b = S.init_b();
        mat w = S.init_w();
        vec X = S.init_X();
        vec a = S.init_a();
        Bruteforce* B = new Bruteforce(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM,interactionswitch);
        B->best_params(myfile,myfile2,gamma,a,b,w,X);
        Impsamp* Imp = new Impsamp(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM,interactionswitch);
//        B->solve(a, b, w, X,myfile, myfile2);
        Imp->best_params(a, b, w, X,myfile,myfile2, gamma, 100);
        myfile2.close();
        myfile.close();
        delete Imp;
    }

//    double exact_vint = 1./sqrt(2) + 1./sqrt(18)+1./sqrt(8);
//    vec Xtry = {0, 0, 1, 1, 3, 3, 7,7};
//    cout << scientific << exact_vint << endl;
//    cout << scientific << S.calc_interaction(Xtry) << endl;

    gamma = 0.2;
    //B->best_params(myfile,myfile2,gamma,a,b,w,X);

    //B->solve(myfile, myfile2);
    //Imp->langevin(myfile,myfile2);

    //G->sample_gibbs(myfile,myfile2);


    //delete Int;
//    delete G;
//    delete Imp;
//    delete B;
    //}


    //}
    //}
}
