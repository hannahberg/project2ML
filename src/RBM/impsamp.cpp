/*
#include "impsamp.h"
#include "solver.h"

Impsamp::Impsamp(
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

vec Impsamp::langevin(std::ofstream &myfile, ofstream &myfile4, double alphanow){
    start=clock();
    double D = 0.5; //diffusion coefficient
    int num_alpha = 0;
    vec alpha_vec = ones(1);
    double Ddt = D*dt;
    double Ddt05 = Ddt*0.5;
    //double current_alpha = alphanow;
    double sumKE = 0;
    double sum_d_wf = 0;
    double sum_d_wf_E = 0;
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(0.,0.5);
    static uniform_real_distribution<double> doubleRNG(0,1);

    // loop over alpha when we try out
    double sdt = sqrt(dt);
    double alpha4 = 0.5*(-4);
    double current_alpha = 0.5;
    // initialize random positions
    mat R3 = init_pos_gaus();
    mat R3new = R3;
    int i; int j; int q;
    mat Fq = F(R3,current_alpha);
    //initialize expectation values
    mat R3plus = zeros(N,dim);
    mat R3minus = zeros(N,dim);
    double accept = 0;
    mat Fqnew = Fq;
    double greens;
    // iterate over MC cycles

    for(i=0;i<mc;i++){
        double bommelom = 0;
        //propose a new position Rnew(boson_j) by moving one boson from position R(boson_j) one at the time
        for(j=0;j<N;j++){
            greens = 0;
            for(q=0;q<dim;q++){
                R3new(j,q) = R3(j,q) + Ddt*Fq(j,q) + gaussianRNG(genMT64)*sdt;
                Fqnew(j,q) = alpha4*R3new(j,q);
                greens += 0.5*(Fq(j,q) + Fqnew(j,q))*(Ddt05*(Fq(j,q)-Fqnew(j,q))-R3new(j,q)+R3(j,q));
            }
            greens = exp(greens);
            double A = (wavefunc(R3new,current_alpha))/wavefunc(R3,current_alpha);
            A *= A;
            A = A * greens;
            // test if new position is more probable than random number between 0 and 1.
            if((A > 1) || (A > doubleRNG(genMT64))){
                R3(j) = R3new(j); //accept new position
                Fq(j) = Fqnew(j);
                accept += 1;
            }else {
                R3new(j) = R3(j);
                Fqnew(j) = Fq(j);
            }
            // calculate change in energy
            double deltakinE = energy_real(R3, current_alpha); // energy_real can also be used? that makes no sense though, since its already perfect..
            //myfile4 << scientific << deltakinE << endl;
            double dwf = d_wavefunc(R3new,current_alpha);
            sumKE += deltakinE;
            bommelom += deltakinE;
            sum_d_wf += dwf;
            sum_d_wf_E += dwf*deltakinE;
            }
    //
        //myfile4 << scientific << bommelom/N << endl;
    }
    num_alpha += 1;
    double mean_KE = sumKE/(N*mc);
    double mean_d_wf = sum_d_wf/(N*mc);
    double mean_d_wf_E = sum_d_wf_E/(N*mc);

    end=clock();
    cout << "Impsamp finished! " << endl;
    vec mean_values = zeros(3);
    mean_values(0) = mean_KE;
    mean_values(1) = mean_d_wf;
    mean_values(2) = mean_d_wf_E;

    //myfile << scientific << dt << " "  <<scientific << mean_KE << " " << scientific << accept/(mc*N) << " " << scientific << ((double)end-(double)start)/CLOCKS_PER_SEC << "    " << 2 << "  # Impsamp" << endl;
    myfile << scientific << mean_KE << endl;
    return mean_values;
}

double Impsamp::energy_impsamp(const mat &R, double alpha){
    double Vext = 0;
    double r2 = 0;
    int i; int j;
    double energy = 0;
    for(i = 0; i < N; i++){
        for(j = 0; j < dim; j++){
            energy += R(i,j)*R(i,j);
        }
    }

    double c = 0.5*m*omega*omega;
    double Ek = (c - 2*alpha*alpha)*energy + alpha*dim*N;
    for(int j = 0; j < N; j++) {
        r2 = 0;
        for(int q = 0; q < dim; q++) {
            r2 += R(j,q)*R(j,q);
        }
        Vext += c*r2; //calculate potential energy
    }
    return Ek + Vext;
}

double Impsamp::best_alpha(){
    ofstream alphafile;ofstream alphafile2;
    alphafile.open("alpha_start08.dat");
    alphafile2.open("alphaconverge_impsamp08.dat");
    double alpha_the_best;
    vec mean_values;
    vec alpha_ = zeros(1000);
    alpha_(0) = 0.8;
    double gamma = 0.2; //hilsen alocias 0.02
    double tol = 0.01;
    for(int i=0;i<999;i++){
        mean_values = langevin(alphafile, alphafile2, alpha_(i));
        double mean_EK = mean_values(0);
        double mean_d_wf = mean_values(1);
        double mean_d_wf_E = mean_values(2);
        alpha_(i+1) = alpha_(i)- gamma*2*(mean_d_wf_E - mean_EK*mean_d_wf);
        alphafile2 << setprecision(12) << alpha_(i) << endl;
        if(abs(alpha_(i+1) - alpha_(i)) < tol){
            alpha_the_best = alpha_(i+1);
            //alphafile << i << endl; //checking if it oscillates in the bottom

        }
    }
    //alphafile << alpha_the_best << endl;
    alphafile.close();
    alphafile2.close();
    return alpha_the_best;
}
*/

#include "solver.h"
#include "impsamp.h"

using namespace arma;

Impsamp::Impsamp(
                       double s_hbar,
                       double mass,
                       double s_omega,
                       double s_rho,
                       int s_mc,
                       int s_N,
                       int s_dim,
                       double s_dt,
                       double sig,
                       int s_H,
                       double s_M,
                       bool s_interact)
:
    Solver(s_hbar, mass,s_omega, s_rho, s_mc, s_N, s_dim, s_dt, sig, s_H, s_M, s_interact)
{}

double Impsamp::energy_analytic(){
    return 0.5 * N * dim;
}
void Impsamp::go_imp(std::ofstream &myfile, ofstream &myfile2){
    vec a = init_a();
    vec b = init_b();
    mat w = init_w();
    vec X = init_X();

    langevin(a,b,w,X,myfile, myfile2);
}

double Impsamp::langevin(const vec &a, const vec &b, const mat &w,const vec &Xin,std::ofstream &myfile, ofstream &myfile2){
    double energy = energy_analytic();
    myfile << "# dim = " << dim << ", N = " << N << ", dt = " << dt << " and mc = " << mc << endl << endl;
    myfile << scientific << "# Theoretical Energy = " << energy << endl << endl;

    start=clock();
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(0.,0.5);
    static uniform_real_distribution<double> doubleRNG(0,1);
    double sdt = sqrt(dt);

    double newE = 0;
    vec X = Xin;
    vec Xnew = X;
    int i; int j;
    double accept = 0;
    double greens;
    vec F = drift(b,X,w,a);
    vec Fnew = F;
    double D = 0.5; // diffusion coefficient
    double Ddt = D*dt;
    double Ddt05 = 0.5*D*dt;
    double sigma_2 = 2/(sigma*sigma);
    double bajsen; double bajs; double A;
    vec sum_d_wf = zeros(M); vec sum_d_wf_E = zeros(M);
    vec sum_d_wf_b = zeros(H); vec sum_d_wf_E_b = zeros(H);
    vec sum_d_wf_a = zeros(M); vec sum_d_wf_E_a = zeros(M);
    mat sum_d_wf_w = zeros(M,H); mat sum_d_wf_E_w = zeros(M,H);
    vec dwfa; vec dwfb; mat dwfw;
    for(i = 0; i < mc; i++){
        bajsen = 0;

        for(j = 0; j < M; j++){
            Xnew(j) = X(j) + Ddt*F(j) + gaussianRNG(genMT64)*sdt;
            Fnew(j) = sigma_2*(-X(j)+a(j) + drifti(b,X,w,j));

            A = (wavefunc(a,b,w,X))/wavefunc(a,b,w,Xnew);
            greens = 0.5*(F(j)+Fnew(j))*(Ddt05*(F(j)-Fnew(j))-Xnew(j)+X(j));
            greens = exp(greens);
            A *= A;
            A = A*greens;
            // test if new position is more probable than random number between 0 and 1.
            if((A > 1) || (A > doubleRNG(genMT64))){
                X(j) = Xnew(j); //accept new position
                F(j) = Fnew(j);
                accept += 1;
            } else {
                Xnew(j) = X(j);
                Fnew(j) = F(j);
            }


            bajs = E_L(a,b,w,X); // local energy
            newE += bajs; // calculate change in energy
            bajsen += bajs;

            dwfa = grad_ai(X,a);
            sum_d_wf_a += dwfa;
            sum_d_wf_E_a += dwfa*bajs;

            dwfb = grad_bj(b,X,w);
            sum_d_wf_b += dwfb;
            sum_d_wf_E_b += dwfb*bajs;

            dwfw = grad_wij(b,X,w);
            sum_d_wf_w += dwfw;
            sum_d_wf_E_w += dwfw*bajs;
       }

    myfile2 << scientific << bajsen/M << endl;
    }

    //cout << "Impsampt finished! The end is near <3" << endl;

    /*
    double energy = energySum/(mc * N);
    double totalenergy = energySum/mc;
    double energySquared = energySquaredSum/(mc * N);
    */
    double E_ = newE/(mc*M);
    vec mean_d_wf_a = sum_d_wf_a/(M*mc);
    vec mean_d_wf_E_a = sum_d_wf_E_a/(M*mc);
    vec mean_d_wf_b = sum_d_wf_b/(M*mc);
    vec mean_d_wf_E_b = sum_d_wf_E_b/(M*mc);
    mat mean_d_wf_w = sum_d_wf_w/(M*mc);
    mat mean_d_wf_E_w = sum_d_wf_E_w/(M*mc);
    calcg1(mean_d_wf_a,mean_d_wf_b,mean_d_wf_w);
    calcg2(mean_d_wf_E_a,mean_d_wf_E_b,mean_d_wf_E_w);
    end=clock();
    myfile << "# Energy" << "     " << "Acceptance" << "   " << "CPU time" << "        " << "Solver" << endl;
    myfile << scientific << E_ << " " << scientific << accept/(mc*M) << " " << scientific << ((double)end-(double)start)/CLOCKS_PER_SEC << "    " << 1 << "  # impsamp" << endl;
    return E_;
}

rowvec Impsamp::best_params(vec a, vec b, mat w, vec X,std::ofstream &myfile, ofstream &myfile2, double gamma, int lol){// gamma is learning rate <3
    ofstream afile; ofstream afile2;

//    vec b = init_b();
//    mat w = init_w();
//    vec X = init_X();
//    vec a = init_a();
    int MHMH = M+H+M*H;
    afile.open("imp_best_N2_D2.dat");
    afile2.open("imp_converge_N2_D2.dat");
    mat alphamat = zeros(lol,MHMH);
    mat startalpha = mat(init_alpha(a,b,w));
//    startalpha.print();
    alphamat.row(0) = startalpha;

    double mean_EL;
    rowvec g1; rowvec g2; rowvec alphanow;
    for(int r=0;r<lol-1;r++){
        mean_EL = langevin(a, b, w, X, myfile, myfile2);

        g1 = getG1();
        g2 = getG2();
        alphamat.row(r+1) = alphamat.row(r) - gamma*2*(g2 - mean_EL*g1);
        alphanow = alphamat.row(r+1);
        afile2 << setprecision(12) << mean_EL << endl;

        //need to reconstruct
        int k = 0;
        for(int i=0;i<M;i++){
            for(int j=0;j<H;j++){
                b(j) = alphanow(j+M);
                w(i,j) = alphanow(M+H+k);
                k++;
            }
        }
        cout << r << endl;
    }
    afile << a << endl;
    afile << b << endl;
    afile << w << endl;
    afile.close();
    afile2.close();
    return alphanow;
}

