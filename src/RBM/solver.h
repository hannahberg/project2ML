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
    Solver(double s_hbar, double mass, double s_omega, double s_rho, int s_mc, int s_N, int s_dim, double s_h, double s_dt, double sig, double s_M, double s_H);
    double hbar;
    int N; //number of particles
    double a_h0;
    double B;
    double omega;
    double m;
    double h;
    double h2;
    double dt;
    double sigma;
    double M;
    double H;
    int mc; //num MC cycles
    int dim;
    double rho; //position update parameter

    // functions in class
    double wavefunc(vec a, vec b, mat w, vec X);
    double d_wavefunc(const mat &R, double alpha_);
    mat init_pos_gaus();
    mat distance_part(const mat &R);
    double energy_num(const mat &R, double alpha);
    mat F(const mat &R_, double alpha_);
    double energy_real(const mat &R, double alpha); //
    vec init_a();
    vec init_b();
    vec init_X();
    mat init_w();
    rowvec init_alpha(const vec &a, const vec &b, const mat &w);

    double u(double bj, const vec &X, const mat &wj);

    vec grad_ai(const vec &X,const vec &a);

    vec grad_bj(const vec &b,const vec &X, const mat &w);

    mat grad_wij(const vec &b,const vec &X, const mat &w);

    double E_L(const vec &a, const vec &b, const mat &w, const vec &X);
    const rowvec& getG1();
    const rowvec& getG2();

    void calcg1(const vec &mean_d_wf_a, const vec &mean_d_wf_b,const mat &mean_d_wf_w);
    void calcg2(const vec &mean_d_wf_E_a, const vec &mean_d_wf_E_b,const mat &mean_d_wf_E_w);

private:

    rowvec G1;
    rowvec G2;
};
#endif
