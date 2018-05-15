#ifndef GIBBS_H
#define GIBBS_H
#include "solver.h"

class Gibbs: public Solver{
public:
    Gibbs(double s_hbar, double mass, double s_omega, double s_rho, int s_mc, int s_N, int s_dim, double s_dt, double sig, int s_H, double s_M);
    //vec langevin(std::ofstream &myfile, std::ofstream &myfile4, double alphanow);
    vec init_h_bool();
    double prob(const vec& X, double bj, const mat &wj);
    void sample_gibbs(std::ofstream &myfile, ofstream &myfile2);
    double get_mu(int i, vec const &hid, const mat &w);
    double random_mu_std(double mu);

};
#endif
