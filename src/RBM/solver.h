#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <fstream>
#include <random>
#include <armadillo>
#include <cmath>
#include <iomanip>
#include <math.h>
//#include <QProgressBar>
using namespace std;
using namespace arma;

class Solver{
public:
    clock_t start, end;
    Solver(double s_omega, double s_rho, int s_mc, int s_N, int s_dim, double s_dt, double sig, int s_H, bool s_interact, double s_spread);
    int N; //number of particles
    double B;
    double omega;
    double dt;
    double sigma;
    double M;
    int H;
    int mc; //num MC cycles
    int dim;
    double rho; //position update parameter
    bool interact;
    double spread;

    // functions in class
    double wavefunc(const vec &a, const vec &b, const mat &w, const vec &X);
    double d_wavefunc(const mat &R, double alpha_);
    mat init_pos_gaus();
    mat distance_part(const mat &R);
    double energy_num(const mat &R, double alpha);
    mat F(const mat &R_, double alpha_);
    double energy_real(const mat &R, double alpha); //
    vec init_a(double spread);
    vec init_b(double spread);
    vec init_X();
    vec init_X_gaus();
    mat init_w(double spread);
    double energy_analytic();
    double wavefunc_g(vec a, vec b, mat w, vec X);
    rowvec init_alpha(const vec &a, const vec &b, const mat &w);
    double u(double bj, const vec &X, const mat &wj);
    vec grad_ai(const vec &X,const vec &a);
    vec grad_bj(const vec &b,const vec &X, const mat &w);
    mat grad_wij(const vec &b,const vec &X, const mat &w);
    double E_L(const vec &a, const vec &b, const mat &w, const vec &X);
    double ELGibbs(const vec &a, const vec &b, const mat &w, const vec &X);
    const rowvec& getG1();
    const rowvec& getG2();
    void calcg1(const vec &mean_d_wf_a, const vec &mean_d_wf_b,const mat &mean_d_wf_w);
    void calcg2(const vec &mean_d_wf_E_a, const vec &mean_d_wf_E_b,const mat &mean_d_wf_E_w);
    vec drift(const vec &b, const vec &X, const mat &w, const vec &a);
    double drifti(const vec &b, const vec &X, const mat &w, int k);
    double calc_interaction(const vec &X);
    double I(const vec &a, const vec &b, const mat &w,const vec &X);
    double II(const vec &b, const mat &w, const vec &X);
    double Mmc;

private:
    rowvec G1;
    rowvec G2;
    double halfsigma;
    double sig2;
    double halfomega;
    double Msig;
};
#endif
