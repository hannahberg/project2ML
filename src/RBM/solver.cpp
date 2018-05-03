#include "solver.h"

Solver::Solver(
               double s_hbar,
               double mass,
               double s_omega,
               double s_rho,
               int s_mc,
               int s_N,
               int s_dim,
               double s_h,
               double s_dt,
               double sig,
               double s_M,
               double s_H){
    hbar = s_hbar;
    m = mass;
    omega = s_omega;
    a_h0 = sqrt(hbar/(m*omega));
    rho = s_rho;
    mc = s_mc;
    N = s_N;
    dim = s_dim;
    h = s_h;
    h2 = 1.0/(h*h);
    dt = s_dt;
    sigma = sig;
    M = s_M;
    H = s_H;
}

vec Solver::init_a(){
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(0,0.001);
    vec a = zeros(M);
    for(int i=0; i<M; i++){
        a(i) = gaussianRNG(genMT64);
    }
    return a;
}

vec Solver::init_b(){
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(0,0.001);
    vec b = zeros(H);
    for(int i=0; i<H; i++){
        b(i) = gaussianRNG(genMT64);
    }
    return b;
}

mat Solver::init_w(){
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(0,0.001);
    mat w = zeros(M,H);
    for(int i=0; i<M; i++){
        for(int j=0; j<H; i++){
            w(i,j) = gaussianRNG(genMT64);
    }
    return w;
}
}

vec Solver::init_X(){
    random_device rd;
    mt19937_64 genMT64(rd());
    uniform_real_distribution<double> doubleRNG(-0.5,0.5);
    vec X = zeros(H);
    for(int i=0; i<H; i++){
        X(i) = doubleRNG(genMT64);
    }
    return X;

}

double Solver::wavefunc(vec a, vec b, mat w, vec X){
    //bool interact = y/n ??
    int i; int j;
    double g = 0;
    double f = 1;
        for(i=0;i<N;i++){
            for(int l=i+1;l<N;l++){
                if(i!=l){
                    f*=(1 + exp(b(j) + (X(i)*(w(i,j)/(sigma*sigma)))));
                for(j=0;j<dim;j++){
                        g += pow((X(i) - a(i)),2)/(2*sigma*sigma);
                    }
                }
            }
        }

    double psi = exp(g)*f;
    return psi;
}

/*
mat Solver::init_pos_gaus(){
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(0.,1);
    double sdt = sqrt(dt);
    int k; int l;
    mat position = zeros(N,dim);
    for(k=0;k<N;k++){
        for(l=0;l<dim;l++){
            position(k,l) = gaussianRNG(genMT64)*sdt;
        }
    }
    return position;
}
mat Solver::distance_part(const mat &R){
    mat rij = zeros(N,N);
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){
            rij(i,j) = norm(R.row(i) - R.row(j));
            rij(j,i) = rij(i,j);

            //rij(i,j) = absdistance(R(i),R(j));
        }
    }
    return rij;
}
*/

double Solver::u(double bj, const vec &X, const mat &wj){
    double sum = 0;
    for(int i = 0; i < M; i++){
        sum += X(i)*wj(i);
    }
    return bj + sum;
}

double Solver::grad_ai(double Xi,double ai, double sigma2){
    return (Xi - ai)/sigma2;
}

double Solver::grad_bj(double bj, const vec &X, const mat &wj){
    double uj = u(bj,X,wj);
    return 1/(1+exp(-uj));
}

double Solver::grad_wij(double Xi, double sigma2, double bj, const vec &X, double wj){
   return Xi*sigma2*grad_bj(bj, X, &wj);
}






