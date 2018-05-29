#ifndef GIBBS_H
#define GIBBS_H
#include "solver.h"

class Gibbs: public Solver{
public:
    Gibbs(double s_omega, double s_rho, int s_mc, int s_N, int s_dim, double s_dt, double sig, int s_H, bool s_interact, double s_spread);
    vec init_h_bool();
    double prob(const vec& X, double bj, const mat &wj);
    double sample_gibbs(const vec &a, const vec &b, const mat &w, std::ofstream &myfile, ofstream &myfile2);
    double get_mu(int i, vec const &hid, const mat &w);
    double random_mu_std(double mu);
    rowvec best_params(std::ofstream &myfile, ofstream &myfile2, double gamma, vec a, vec b, mat w, int gdc);

};
#endif
