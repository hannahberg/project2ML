#ifndef IMPSAMP_H
#define IMPSAMP_H
#include "solver.h"

class Impsamp: public Solver{
public:
    Impsamp(double s_beta, double s_hbar, double mass, double s_omega, double s_alpha, double s_rho, int s_mc, int s_N, int s_dim, double s_h, double s_dt);
    vec langevin(std::ofstream &myfile, std::ofstream &myfile4, double alphanow);
    double energy_impsamp(const mat &R, double alpha);
    double best_alpha();

};
#endif
