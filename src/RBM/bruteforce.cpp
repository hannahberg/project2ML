#include "solver.h"
#include "bruteforce.h"

using namespace arma;

Bruteforce::Bruteforce(double s_beta,
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

double Bruteforce::energy_local(){
    return 0.5 * hbar * omega * N * dim;
}

void Bruteforce::solve(std::ofstream &myfile, ofstream &myfile2){
    double energy = energy_local();

    myfile << "# dim = " << dim << ", N = " << N << ", dt = " << dt << ", alpha = " << alpha << " and mc = " << mc << endl << endl;
    myfile << scientific << "# Theoretical Energy = " << energy << endl << endl;

    start=clock();
    random_device rd;
    mt19937_64 genMT64(rd());
    uniform_real_distribution<double> doubleRNG(0,1);

    int num_alpha = 0;
    vec alpha_vec = ones(1);
    double current_alpha;
    double newE = 0;
    mat testR = init_pos();
    mat dis = distance_part(testR);
    cout << testR << endl;
    cout << dis << endl;
    current_alpha = alpha;
    mat R = init_pos(); // initialize random positions
    mat Rnew = R;
    int i; int j; int q;
    double accept = 0;

    for(i=0;i<mc;i++){ // iterate over MC cycles
        double bajsen = 0;
        //propose a new position Rnew(boson_j) by moving one boson from position R(boson_j) one at the time
        for(j=0;j<N;j++){
            for(q=0;q<dim;q++){
                Rnew(j,q) = R(j,q) + (doubleRNG(genMT64) - 0.5)*rho;
            }

            double A = (wavefunc(Rnew,current_alpha))/wavefunc(R, current_alpha);
            A *= A;

            // test if new position is more probable than random number between 0 and 1.
            if((A > 1) || (A > doubleRNG(genMT64))){
                R(j) = Rnew(j); //accept new position
                accept += 1;
            } else {
                Rnew(j) = R(j);
            }
            double bajs = energy_real(R, current_alpha);
            newE += bajs; // calculate change in energy
            bajsen += bajs;
       }

    myfile2 << scientific << bajsen/N << endl;
    }

    num_alpha += 1;
    cout << "Analytical finished! Hang in there <3" << endl;

    /*
    double energy = energySum/(mc * N);
    double totalenergy = energySum/mc;
    double energySquared = energySquaredSum/(mc * N);
    */

    end=clock();
    myfile << "# Energy" << "     " << "Acceptance" << "   " << "CPU time" << "        " << "Solver" << endl;
    myfile << scientific << newE/(mc*N) << " " << scientific << accept/(mc*N) << " " << scientific << ((double)end-(double)start)/CLOCKS_PER_SEC << "    " << 0 << "  # Analytic" << endl;
}

mat Bruteforce::init_pos(){
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static uniform_real_distribution<double> doubleRNG(0,1);
    int k; int l;
    mat position = zeros(N,dim);
    for(k=0;k<N;k++){
        for(l=0;l<dim;l++){
            //position(k,l) = (rando() - 0.5)*rho;
            position(k,l) = (doubleRNG(genMT64) - 0.5)*rho;
        }
    }
    return position;
}
void Bruteforce::solve_num( std::ofstream &myfile, std::ofstream &myfile3){
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static uniform_real_distribution<double> doubleRNG(0,1);
    start=clock();

    // loop over alpha when we try out
    int num_alpha = 0;
    vec alpha_vec = ones(1);
    double current_alpha;
    double sumKE = 0;
    double sumKe2 = 0;
    current_alpha = alpha;//alpha_vec(num_alpha);
    // initialize random positions
    mat R2 = init_pos();
    mat R2new = R2;
    int i; int j; int q;

    //initialize expectation values
    double accept = 0;

        // iterate over MC cycles
    for(i=0;i<mc;i++){
        //propose a new position Rnew(boson_j) by moving one boson from position R(boson_j) one at the time
        for(j=0;j<N;j++){
            for(q=0;q<dim;q++){
                //R2new(j,q) = R2(j,q) + (rando() - 0.5)*rho;
                //cout << R2new(j,q) << endl;
                R2new(j,q) = R2(j,q) + (doubleRNG(genMT64) - 0.5)*rho;
            }

            double A = (wavefunc(R2new,current_alpha))/wavefunc(R2,current_alpha);
            A = abs(A);
            A *= A;
            // test if new position is more probable than random number between 0 and 1.
            if((A > doubleRNG(genMT64))) {

                R2(j) = R2new(j); //accept new position
                accept += 1;
            } else {
                R2new(j) = R2(j);
            }
            // calculate change in energy
            double drit = energy_num(R2, current_alpha);


            //sumKE += energy_num(R2, current_alpha);
            sumKE += drit;
            sumKe2 += drit*drit;

            myfile3 << scientific << drit << " " << wavefunc(R2new,current_alpha) << endl;

        }


    }
    num_alpha += 1;
    end=clock();
    double var = sumKe2/(mc*N) - pow(sumKE/(mc*N),2);
    cout << "Numerical energy is finished, yay!" << endl;
    myfile << scientific << sumKE/(mc*N) << " " << scientific << accept/(mc*N) << " var:" << var << " " << scientific << ((double)end-(double)start)/CLOCKS_PER_SEC << "    " << 1 << "  # Numerical" << endl;
}

