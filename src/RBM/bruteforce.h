#ifndef BRUTEFORCE_H
#define BRUTEFORCE_H

#include "solver.h"

class Bruteforce : public Solver {

public:
    Bruteforce(double s_hbar, double mass, double s_omega, double s_rho, int s_mc, int s_N, int s_dim, double s_h, double s_dt);
    double energy_analytic();
    void solve(std::ofstream &myfile, std::ofstream &myfile2);// int mc=10, int N=1
};


#endif
