#include "solver.h"

Solver::Solver(double s_beta,
               double s_hbar,
               double mass,
               double s_omega,
               double s_alpha,
               double s_rho,
               int s_mc,
               int s_N,
               int s_dim,
               double s_h,
               double s_dt
               double sig
               double s_M
               double s_H){
    beta = s_beta;
    hbar = s_hbar;
    m = mass;
    omega = s_omega;
    a_h0 = sqrt(hbar/(m*omega));
    alpha = s_alpha;
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
    for(int i>0; i=M; i++){
        a(i) = gaussianRNG(genMT64);
    }
    return a;
}

vec Solver::init_b(){
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(0,0.001);
    vec b = zeros(H);
    for(int i>0; i=H; i++){
        b(i) = gaussianRNG(genMT64);
    }
    return b;
}

mat Solver::init_w(){
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(0,0.001);
    mat w = zeros(M,H);
    for(int i>0; i=M; i++){
        for(int j>0; j=H; i++){
            w(i,j) = gaussianRNG(genMT64);
    }
    return w;
}

vec Solver::init_X(){
    random_device rd;
    mt19937_64 genMT64(rd());
    uniform_real_distribution<double> doubleRNG(-0.5,0.5);
    vec X = zeros(H);
    for(int i>0; i=H; i++){
        X(i) = doubleRNG(genMT64);
    }
    return b;

}

double Solver::wavefunc(const mat &R, double alpha_){
    //bool interact = y/n ??
    int i; int j;
    double g = 0;

        for(i=0;i<N;i++){
            for(int l=i+1;l<N;l++){
                if(i!=l){
                    f*=(1 + exp(b(j) + (X(i)*(w(i,j)/(sigma*sigma)))));
                for(j=0;j<dim;j++){
                    //if(j==2){
                    //    g += beta*beta*newR(i,j)*newR(i,j);
                    //} else{
                        g += pow((X(i) - a(i)),2)/(2*sigma*sigma)
                    }
                }
            }
            }
        }


    double psi = exp(-alpha_*g)*f;
    return psi;





}

double Solver::d_wavefunc(const mat &R, double alpha_){
    //bool interact = y/n ??
    int i; int j;
    double g = 0;
    if(dim==1){
        for(i=0;i<N;i++){
            g += R(i)*R(i); // take Product of Pi(g(Ri)
        }
    } else{
        for(i=0;i<N;i++){
            for(j=0;j<dim;j++){
                g += R(i,j)*R(i,j);//g += dot(R.row(i),R.row(i));
            }
        }
    }
    double f = 1; //no interaction here!!
    double psi = exp(-alpha_*g)*f;
    return psi*(-g);
}
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

double Solver::energy_num(const mat &R, double alphanow){

    double wavefuncnow = wavefunc(R, alphanow);
    double Ek = 0;
    double Vext = 0;
    double r2 = 0;
    // Calculate kinetic energy by numerical derivation
    mat Rplus = R + h;
    mat Rminus = R - h;

    double c = 0.5*m*omega*omega;
    for(int j = 0; j < N; j++) {
        r2 = 0;
        for(int q = 0; q < dim; q++) {
            r2 += R(j,q)*R(j,q);
        }
        Vext += c*r2; //calculate potential energy
    }
    double wavefuncplus = wavefunc(Rplus, alphanow);
    double wavefuncminus = wavefunc(Rminus, alphanow);
    Ek -= (wavefuncplus+wavefuncminus - 2*wavefuncnow);
    Ek = 0.5 * Ek * h2 / wavefuncnow;
    return Ek + Vext;
}

mat Solver::F(const mat &R_, double alpha_) {
    return -4*R_*alpha_;
}

double Solver::energy_real(const mat &R, double alpha){ //done optimization
    int i; int j;
    double energy = 0;
    for(i = 0; i < N; i++){
        for(j = 0; j < dim; j++){
            energy += R(i,j)*R(i,j);
        }
    }
    double en = (0.5*omega*omega - 2*alpha*alpha)*energy + alpha*dim*N;
    return en;

}

