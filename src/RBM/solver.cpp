#include "solver.h"

Solver::Solver(double s_hbar,
               double mass,
               double s_omega,
               double s_rho,
               int s_mc,
               int s_N,
               int s_dim,
               double s_dt,
               double sig,
               int s_H,
               int s_M){
    hbar = s_hbar;
    m = mass;
    omega = s_omega;
    a_h0 = sqrt(hbar/(m*omega));
    rho = s_rho;
    mc = s_mc;
    N = s_N;
    dim = s_dim;
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
    for(int j=0; j<H; j++){
        b(j) = gaussianRNG(genMT64);
    }
    return b;
}

mat Solver::init_w(){
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(0,0.001);
    mat w = zeros(M,H);
    for(int i=0; i<M; i++){
        for(int j=0; j<H; j++){
            w(i,j) = gaussianRNG(genMT64);
        }
    }
    return w;
}

vec Solver::init_X(){
    random_device rd;
    mt19937_64 genMT64(rd());
    uniform_real_distribution<double> doubleRNG(-0.5,0.5);
    vec X = zeros(M);
    for(int i=0; i<M; i++){
        X(i) = doubleRNG(genMT64);
    }
    return X;
}

vec Solver::init_X_gaus(){
    random_device rd;
    mt19937_64 genMT64(rd());
    normal_distribution<double> gaussianRNG(0.,0.5);
    vec X = zeros(M);
    for(int i=0; i<M; i++){
        X(i) = gaussianRNG(genMT64);
    }
    return X;
}
rowvec Solver::init_alpha(const vec &a, const vec &b, const mat &w){
    rowvec alpha = zeros<rowvec>(M+H+(M*H));
    int i; int j; int k;
    for(i=0;i<M;i++){
        alpha(i) = a(i);
    }
    for(int j=M;j<(M+H);j++){
        alpha(j) = b(j-M);
    }
    k = 0;
    for(i=0;i<M;i++){
        for(j=0;j<H;j++){
            alpha(k+M+H) = w(i,j);
            k++;
        }
    }
    return alpha;
}


double Solver::wavefunc(vec a, vec b, mat w, vec X){
    double g = 0;
    double f = 1;
    for(int i=0;i<M;i++){
        g += pow((X(i) - a(i)),2)/(2*sigma*sigma);
    }
    for(int j=0;j<H;j++){
        f*=(1 + exp(u(b(j), X, w.col(j))));
    }
    double psi = exp(g)*f;
    return psi;
}

double Solver::wavefunc_g(vec a, vec b, mat w, vec X){
    return sqrt(wavefunc(a,b,w,X));
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

vec Solver::grad_ai(const vec &X,const vec &a){
    vec grada = zeros(M);
    for(int i=0; i < M; i++){
        grada(i) = X(i) - a(i);
    }
    return grada/(sigma*sigma);
}

vec Solver::grad_bj(const vec &b, const vec &X, const mat &w){
    vec gradb = zeros(H);
    double uj;
    for (int j=0; j < H; j++){
        uj = u(b(j),X,w.col(j));
        gradb(j) = 1+exp(-uj);
    }
    return 1/gradb;
}


mat Solver::grad_wij(const vec &b,const vec &X, const mat &w){
    mat gradw = zeros(M,H);
    double uj;
    for(int i=0; i < M; i++){
        for (int j=0; j < H; j++){
            uj = u(b(j),X,w.col(j));
            gradw(i,j) = X(i)/(1+exp(-uj));
        }
    }
    return gradw/(sigma*sigma);
}

const rowvec& Solver::getG1(){
    return G1;
}

const rowvec& Solver::getG2(){
    return G2;
}


double Solver::E_L(const vec &a, const vec &b, const mat &w,const vec &X){
    double energysum = 0;
    double r2 = 0;
    double temp_I = 0;
    double u_ = 0;
    double sigma2_ = 1/(sigma*sigma);
    double I = 0;
    double II = 0;
    double temp_II = 0;
    double eu = 0;
    double wmj = 0;
    double Xm = 0;
    for(int m = 0; m < M; m++){
        Xm = X(m);
        I = -(Xm-a(m));
        r2 += Xm*Xm;
        temp_I = 0;
        temp_II = 0;
        II = 0;
        for(int j = 0; j < H; j++){
            wmj = w(m,j);
            u_ = u(b(j), X, w.col(j));
            eu = exp(u_);
            temp_I += wmj/(1+1/eu);
            temp_II += wmj*wmj*eu/((1+eu)*(1+eu));
        }
        II = sigma2_*sigma2_*temp_II;
        I *= I + sigma2_*temp_I;
        energysum += I+II;
    }
    return -0.5*(energysum - M*sigma2_) + 0.5*r2*omega*omega;

}

double Solver::ELGibbs(const vec &a, const vec &b, const mat &w,const vec &X){
    double energysum = 0;
    double r2 = 0;
    double temp_I = 0;
    double u_ = 0;
    double sigma2_ = 1/(sigma*sigma);
    double I = 0;
    double II = 0;
    double temp_II = 0;
    double eu = 0;
    double wmj = 0;
    double Xm = 0;
    for(int m = 0; m < M; m++){
        Xm = X(m);
        I = -(Xm-a(m));
        r2 += Xm*Xm;
        temp_I = 0;
        temp_II = 0;
        II = 0;
        for(int j = 0; j < H; j++){
            wmj = w(m,j);
            u_ = u(b(j), X, w.col(j));
            eu = exp(u_);
            temp_I += wmj/(1+1/eu);
            temp_II += wmj*wmj*eu/((1+eu)*(1+eu));
        }
        II = sigma2_*sigma2_*temp_II;
        I *= I + sigma2_*temp_I;
        energysum += 0.25*I+0.5*II;
    }
    return -0.5*(energysum - 0.5*M*sigma2_) + 0.5*r2*omega*omega;
}


void Solver::calcg1(const vec &mean_d_wf_a, const vec &mean_d_wf_b,const mat &mean_d_wf_w){
    G1 = init_alpha(mean_d_wf_a,mean_d_wf_b,mean_d_wf_w);
}


void Solver::calcg2(const vec &mean_d_wf_E_a, const vec &mean_d_wf_E_b,const mat &mean_d_wf_E_w){
    G2 = init_alpha(mean_d_wf_E_a,mean_d_wf_E_b,mean_d_wf_E_w);
}


vec Solver::drift(const vec &b, const vec &X, const mat &w, const vec &a){
    vec F = zeros(M);
    for(int i = 0; i < M; i++){
        F(i) += - X(i)+ a(i) + drifti(b,X,w,i);
    }
    F = (2/(sigma*sigma))*F;
    return F;
}


double Solver::drifti(const vec &b, const vec &X, const mat &w, int k){
    double sum = 0;
    for(int j = 0; j < H; j++){
        sum += w(k,j)/(1+exp(-u(b(j),X,w.col(j))));
    }
    return sum;
}

double Solver::energy_analytic(){
    return 0.5 * M;
}
