#ifndef IMPSAMP_H
#define IMPSAMP_H
#include "solver.h"

class Impsamp: public Solver{
public:
    Impsamp(double s_omega, double s_rho, int s_mc, int s_N, int s_dim, double s_dt, double sig, int s_H, bool s_interact);
    double energy_impsamp(const mat &R, double alpha);
    double best_alpha();
    double langevin(const vec &a, const vec &b, const mat &w, const vec &Xin, std::ofstream &myfile, ofstream &myfile2);
    double energy_analytic();
    rowvec best_params(std::ofstream &myfile, ofstream &myfile2, double gamma, vec a, vec b, mat w, vec X, int gdc);
    void go_imp(std::ofstream &myfile, ofstream &myfile2);

};
#endif
