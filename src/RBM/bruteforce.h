#ifndef BRUTEFORCE_H
#define BRUTEFORCE_H

#include "solver.h"

class Bruteforce : public Solver {

public:
    Bruteforce(double s_omega, double s_rho, int s_mc, int s_N, int s_dim, double s_dt, double sig, int s_H, bool s_interact);
    double energy_analytic();
    double solve(const vec &a, const vec &b, const mat &w, const vec &X, ofstream &myfile, ofstream &myfile2);// int mc=10, int N=1
    rowvec best_params(std::ofstream &myfile, ofstream &myfile2, double gamma, vec a, vec b, mat w, vec X);
    void go_brute(std::ofstream &myfile, std::ofstream &myfile2);
};


#endif
