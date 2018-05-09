#include "solver.h"
#include "gibbs.h"

using namespace arma;

Gibbs::Gibbs(
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
                       double s_H,
                       double s_M)
:
    Solver(s_hbar, mass,s_omega, s_rho, s_mc, s_N, s_dim, s_h, s_dt, sig, s_H, s_M)
{}

void Gibbs::sample_gibbs(std::ofstream &myfile, ofstream &myfile2){
    double energy = energy_analytic();
    myfile << "# dim = " << dim << ", N = " << N << ", dt = " << dt << " and mc = " << mc << endl << endl;
    myfile << scientific << "# Theoretical Energy = " << energy << endl << endl;
    double sumE = 0;
    start=clock();
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(0.,0.5);
    static uniform_real_distribution<double> doubleRNG(0,1);

    vec a; vec X; int j; double Pj; double E_LGibbs;
    mat w = init_w();
    vec hid = init_h_bool();
    int i;
    double accept = 0;
    double newE;

    a = init_a();
    vec b = init_b()*0.001;
    X = init_X();
    for(int k = 0; k < mc; k++){

        //X = init_X_gaus();
        newE = 0;
        //a = init_a();

        for(i = 0; i < M; i++){ // finding the ideal positions X(i)
            double mu = get_mu(i,hid,w);
            //cout << "lol" << endl;
            X(i) = random_mu_std(mu);


        }

        for(j = 0; j < H; j++){
            Pj = prob(X,b(j),w.col(j));
            //cout << "lol" << endl;
            if(doubleRNG(genMT64)<=Pj){
                hid(j) = 1;
                accept += 1;
            } else {
                hid(j) = 0;
            }
        }
        E_LGibbs = ELGibbs(a,b,w,X);
        //cout << E_LGibbs << endl;
        sumE += E_LGibbs;
        newE += E_LGibbs;
        myfile2 << scientific << newE << endl;
    }
    end=clock();
    myfile << "# Energy" << "     " << "Acceptance" << "   " << "CPU time" << "        " << "Solver" << endl;
    myfile << scientific << sumE/(mc) << " " << scientific << accept/(mc*M) << " " << scientific << ((double)end-(double)start)/CLOCKS_PER_SEC << "    " << 2 << "  # gibbs" << endl;
}

vec Gibbs::init_h_bool(){
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static uniform_int_distribution<long> boolRNG(0,1);
    vec hid = zeros(H);
    for(int j = 0; j < H; j++){
        hid(j) = boolRNG(genMT64);
    }
    return hid;
}

double Gibbs::prob(const vec& X, double bj, const mat &wj){
    return 1/(1+exp(u(bj,X,wj)));
}


double Gibbs::get_mu(int i,vec const &hid, const mat &w){
    double mu = 0;
    for(int j = 0; j < H; j++){
        mu += hid(j)*w(i,j);
    }
    return mu;
}

double Gibbs::random_mu_std(double mu){
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(mu,sigma);
    return gaussianRNG(genMT64);
}
