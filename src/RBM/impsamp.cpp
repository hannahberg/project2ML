#include "impsamp.h"
#include "solver.h"

Impsamp::Impsamp(double s_beta,
                 double s_hbar,
                 double mass,
                 double s_omega,
                 double s_alpha,
                 double s_rho,
                 int s_mc,
                 int s_N,
                 int s_dim,
                 double s_h,
                 double s_dt)
:
    Solver(s_beta, s_hbar, mass,s_omega, s_alpha, s_rho, s_mc, s_N, s_dim, s_h, s_dt)
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
    double alpha4 = alpha*(-4);
    double current_alpha = alphanow;//alpha_vec(num_alpha);
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
