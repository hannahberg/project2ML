#ifndef INTERACT_H
#define INTERACT_H
#include "solver.h"

class Interact : public Solver{
public:
    double a;
    double binsize;
    Interact(double s_beta, double s_hbar, double mass, double s_omega, double s_alpha, double s_rho, int s_mc, int s_N, int s_dim, double s_h, double s_dt);
    mat init_pos_interact();
    mat too_close(mat Rtull);
    double wavefunc_interact(mat R, double alphanow);
    double energy_interact(const mat &R, double alphanow);
    vec solve_interact(std::ofstream &myfile, ofstream &myfile5, double alphanoe);
    double d_wavefunc_interact(const mat &R, double alphanow);

    mat lapphi(const mat &R, double alpha_);
    mat nablaphi(const mat &R, double alpha_);
    mat nablaphinablaF(mat R, mat distR, double alpha_);
    mat nablaf(mat R, mat distR);
    mat doublesum(mat R, mat distanceR);
    mat suma2(const mat &distanceR);
    mat quantumF(const mat &R, double alpha_, mat rij);
    mat nablafsquared(const mat &init_distance, const mat &R);
    mat newnablaf(const mat &init_distance, const mat &R);
    double best_alpha();
    vec dist_origo(const mat& R);
};
#endif
