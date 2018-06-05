#include "solver.h"

Solver::Solver(double s_omega,
               double s_rho,
               int s_mc,
               int s_N,
               int s_dim,
               double s_dt,
               double sig,
               int s_H,
               bool s_interact,
               double s_spread){
    omega = s_omega;
    rho = s_rho;
    mc = s_mc;
    N = s_N;
    dim = s_dim;
    dt = s_dt;
    sigma = sig;
    M = s_N*s_dim;
    H = s_H;
    interact = s_interact;
    halfsigma = 1/(2*sigma*sigma);
    sig2 = 1.0/(sigma*sigma);
    halfomega = 0.5*omega*omega;
    Msig = M*sig2;
    spread = s_spread;
    Mmc = 1.0/(M*mc);
}


vec Solver::init_a(double spread){
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(0,spread);
//    vec a = zeros(M);
    vec a(M);
    for(int i=0; i<M; i++){
        a(i) = gaussianRNG(genMT64);
    }
    return a;
}

vec Solver::init_b(double spread){
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(0,spread);
//    vec b = zeros(H);
    vec b(H);
    for(int j=0; j<H; j++){
        b(j) = gaussianRNG(genMT64);
    }
    return b;
}

mat Solver::init_w(double spread){
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(0,spread);
//    mat w = zeros(M,H);
    mat w(M,H);
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


double Solver::wavefunc(const vec &a,const vec &b,const mat &w, const vec &X){
    double g = 0;
    double f = 1;
    for(int i=0;i<M;i++){
        g += (X(i) - a(i))*(X(i) - a(i));
    }
    for(int j=0;j<H;j++){
        f*=(1 + exp(u(b(j), X, w.col(j))));
    }
    return exp(-g*halfsigma)*f;
}

double Solver::wavefunc_g(vec a, vec b, mat w, vec X){
    return sqrt(wavefunc(a,b,w,X));
}

double Solver::u(double bj, const vec &X, const mat &wj){
    double sum = 0;
    for(int i = 0; i < M; i++){
        sum += X(i)*wj(i);
    }
    return bj + sum;
}

double Solver::calc_interaction(const vec &X){
    double sum = 0;
    double dr; double r2;
    for(int i = 0; i < M-dim; i+=dim){
        for(int j = i+dim; j < M; j+=dim){
            //cout << "in ur mama" << endl;
            r2 = 0;
            for(int d = 0; d < dim; d++){
                //cout << "i " << i + d << " " <<"j " << j+d << endl;
//                cout << X(i+d) << " " <<X(j+d) << endl;
                dr = X(i+d)-X(j+d);
                dr *= dr;
                r2 += dr;
            }
            sum += 1.0/sqrt(r2);
        }
    }
    //sum = 1.0/(sqrt(sum));
//    cout << "1/rij " << sum << endl;
    return sum;

}

vec Solver::grad_ai(const vec &X,const vec &a){
//    vec grada = X - a;
//    for(int i=0; i < M; i++){
//        grada(i) = X(i) - a(i);
//    }
    return (X - a)*sig2;
}

vec Solver::grad_bj(const vec &b, const vec &X, const mat &w){
//    vec gradb = zeros(H);
    vec gradb(H);
    double uj;
    for (int j=0; j < H; j++){
        uj = u(b(j),X,w.col(j));
        gradb(j) = 1/(1+exp(-uj));
    }
    return gradb;
}


mat Solver::grad_wij(const vec &b,const vec &X, const mat &w){
//    mat gradw = zeros(M,H);
    mat gradw(M,H);
    double uj;
    for(int i=0; i < M; i++){
        for (int j=0; j < H; j++){
            uj = u(b(j),X,w.col(j));
            gradw(i,j) = X(i)/(1+exp(-uj));
        }
    }
    return gradw*sig2;
}

const rowvec& Solver::getG1(){
    return G1;
}

const rowvec& Solver::getG2(){
    return G2;
}


double Solver::E_L(const vec &a, const vec &b, const mat &w,const vec &X){
    double energysum = 0;
    double Xm2 = 0;
    double temp_I = 0;
    double u_ = 0;
    double sumI = 0;
    double sumII = 0;
    double temp_II = 0;
    double eu = 0;
    double wmj = 0;
    double Xm = 0;
    for(int m = 0; m < M; m++){
        Xm = X(m);
        sumI = -(Xm-a(m))*sig2;
        Xm2 += Xm*Xm;
        temp_I = 0;
        temp_II = 0;
        sumII = 0;
        for(int j = 0; j < H; j++){
            wmj = w(m,j);
            u_ = u(b(j), X, w.col(j));
            eu = exp(u_);
            temp_I += wmj/(1+exp(-u_));
            temp_II += wmj*wmj*eu/((1+eu)*(1+eu));
        }
        sumI += sig2*temp_I;
        sumI *= sumI;
        sumII = sig2*sig2*temp_II;
        energysum += sumI+sumII;
    }
    double E1 = 0;
    if(interact){
        E1 = calc_interaction(X);
        //cout << "interacting" << endl;

    }
    return -0.5*(energysum - Msig) + halfomega*Xm2+ E1;

}

double Solver::I(const vec &a, const vec &b, const mat &w,const vec &X){
    double energysum = 0;
    double temp_I = 0;
    double u_ = 0;
    double sumI = 0;
    double wmj = 0;
    double Xm = 0;
    for(int m = 0; m < M; m++){
        Xm = X(m);
        sumI = -(Xm-a(m))*sig2;
        temp_I = 0;
        for(int j = 0; j < H; j++){
            wmj = w(m,j);
            u_ = u(b(j), X, w.col(j));
            temp_I += wmj/(1+exp(-u_));
        }
        sumI += sig2*temp_I;
        sumI *= sumI;

        energysum += sumI;
       //cout << energysum << endl;
    }

    return energysum;
}


double Solver::II(const vec &b, const mat &w,const vec &X){
    double energysum = 0;
    double u_ = 0;
    double sumII = 0;
    double temp_II = 0;
    double eu = 0;
    double wmj = 0;
    for(int m = 0; m < M; m++){
        temp_II = 0;
        for(int j = 0; j < H; j++){
            wmj = w(m,j);
            u_ = u(b(j), X, w.col(j));
            eu = exp(u_);
            temp_II += wmj*wmj*eu/((1+eu)*(1+eu));
        }

        sumII = sig2*sig2*temp_II;
        energysum += sumII;
    }
    return energysum;
}

double Solver::ELGibbs(const vec &a, const vec &b, const mat &w,const vec &X){
    double energysum = 0;
    double Xm2 = 0;
    double temp_I = 0;
    double u_ = 0;
    double sumI = 0;
    double sumII = 0;
    double temp_II = 0;
    double eu = 0;
    double wmj = 0;
    double Xm = 0;
    for(int m = 0; m < M; m++){
        Xm = X(m);
        sumI = -(Xm-a(m))*sig2;
        Xm2 += Xm*Xm;
        temp_I = 0;
        temp_II = 0;
        sumII = 0;
        for(int j = 0; j < H; j++){
            wmj = w(m,j);
            u_ = u(b(j), X, w.col(j));
            eu = exp(u_);
            temp_I += wmj/(1+exp(-u_));
            temp_II += wmj*wmj*eu/((1+eu)*(1+eu));
        }

        sumI += sig2*temp_I;
        sumI *= sumI;
        sumII = sig2*sig2*temp_II;
        energysum += 0.25*sumI + 0.5*sumII;
    }
    double E1 = 0;
    if(interact){
        E1 = calc_interaction(X);
        //cout << "interacting" << endl;

    }
    return -0.5*(energysum - 0.5*Msig) + halfomega*Xm2 + E1;
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
//    F = (2/(sigma*sigma))*F;
    return 2*sig2*F;
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
