#include "solver.h"
#include "bruteforce.h"

using namespace arma;

Bruteforce::Bruteforce(
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

double Bruteforce::energy_analytic(){
    return 0.5 * N * dim;
}

void Bruteforce::go_brute(std::ofstream &myfile, ofstream &myfile2){
    vec a = init_a();
    vec b = init_b();
    mat w = init_w();
    vec X = init_X();

    solve(a,b,w,X,myfile, myfile2);
}

double Bruteforce::solve(const vec &a, const vec &b, const mat &w,const vec &X, std::ofstream &myfile, std::ofstream &myfile2){
    double energy = energy_analytic();
    myfile << "# dim = " << dim << ", N = " << N << ", dt = " << dt << " and mc = " << mc << endl << endl;
    myfile << scientific << "# Theoretical Energy = " << energy << endl << endl;

    start=clock();
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static uniform_real_distribution<double> doubleRNG(0,1);
    double sumE = 0;
    vec sum_d_wf = zeros(M);
    vec sum_d_wf_E = zeros(M);
    vec sum_d_wf_b = zeros(M);
    vec sum_d_wf_E_b = zeros(M);
    double newE = 0;
    vec Xnew = X;
    vec Xflex = X;
    int i; int j; int q;
    double accept = 0;

    for(i=0;i<mc;i++){
        double bajsen = 0;
        for(j=0;j<M;j++){
            //for(q=0;q<dim;q++){
                Xnew(j) = Xflex(j) + (doubleRNG(genMT64) - 0.5)*rho;
           // }

            double A = (wavefunc(a,b,w,Xflex))/wavefunc(a,b,w,Xnew);

            A *= A;

            // test if new position is more probable than random number between 0 and 1.
            if((A > 1) || (A > doubleRNG(genMT64))){
                Xflex(j) = Xnew(j); //accept new position
                accept += 1;
            } else {
                Xnew(j) = Xflex(j);
            }


            double localenergy = E_L(a,b,w,Xflex);
            sumE += localenergy; // calculate change in energy

            vec dwfa = grad_ai(Xflex,a);
            sum_d_wf_a += dwfa;
            sum_d_wf_E_a += dwfa*localenergy;

            vec dwfb = grad_bj(b,Xflex,w);
            sum_d_wf_b += dwfb;
            sum_d_wf_E_b += dwfb*localenergy;

            mat dwfw = grad_wij(b,Xflex,w);
            sum_d_wf_w += dwfw;
            sum_d_wf_E_w += dwfw*localenergy;

            }
    //
        //myfile4 << scientific << bommelom/N << endl;
    }

    double mean_E = sumE/(M*mc);
    vec mean_d_wf_a = sum_d_wf/(M*mc);
    vec mean_d_wf_E_a = sum_d_wf_E/(M*mc);
    vec mean_d_wf_b = sum_d_wf/(M*mc);
    vec mean_d_wf_E_b = sum_d_wf_E/(M*mc);
    mat mean_d_wf_w = sum_d_wf/(M*mc);
    mat mean_d_wf_E_w = sum_d_wf_E/(M*mc);

    vec G1 = getG1(); //definer disse i header!!
    vec G2 = getG2();

    for(i=0;i<M;i++){
        G1(i) = mean_d_wf_a(i);
        G2(i) = mean_d_wf_E_a(i);
    }
    for(j=M;j<(M+H);j++){
        G1(j) = mean_d_wf_b(j);
        G2(j) = mean_d_wf_E_b(j);
    }
    int k = 0;
    for(i=0;i<M;i++){
        for(j=0;j<H;j++){
            G1(k+M+H) = mean_d_wf_w(i,j);
            G2(k+M+H) = mean_d_wf_E_w(i,j);
            k++
        }
    }

    myfile2 << scientific << sumE/M << endl;

    cout << "Brute force finished! Hang in there <3" << endl;

    /*
    double energy = energySum/(mc * N);
    double totalenergy = energySum/mc;
    double energySquared = energySquaredSum/(mc * N);
    */

    end=clock();
    myfile << "# Energy" << "     " << "Acceptance" << "   " << "CPU time" << "        " << "Solver" << endl;
    myfile << scientific << mean_E << " " << scientific << accept/(mc*M) << " " << scientific << ((double)end-(double)start)/CLOCKS_PER_SEC << "    " << 0 << "  # Analytic" << endl;
    return mean_E;
}

vec Bruteforce::best_params(std::ofstream &myfile, ofstream &myfile2){
    ofstream afile;ofstream afile2;

    vec b = init_b();
    mat w = init_w();
    vec X = init_X();
    vec a = init_a();

    afile.open("start.dat");
    afile2.open("converge.dat");
    vec alpha_best;
    int lol = 100;
    mat alphamat = zeros(M+H+(M*H),lol)
    alphamat.row(0) = init_alpha();

    double gamma = 0.2; // LEARNING RATE <3


    for(int r=0;r<lol-1;r++){

        double mean_EL = solve(a, b, w, X, myfile, myfile2);
        vec G1 = getG1();
        vec G2 = getG2();
        alphamat.row(r+1) = alphamat.row(r)- gamma*2*(G2 - mean_EL_a*G1;
        vec alphanow = alphamat.row(r+1));
        afile2 << setprecision(12) << mean_EL << endl;

        //need to reconstruct
        vec alphanow = alphamat.row(r+1))
        int k = 0;
        for(int i=0;i<M;i++){
            for(int j=0;j<H;j++){
                b(j) = alphanow(j+M);
                w(i,j) = alphanow(M+H+k);
                k++;
            }
        }

            //afile << i << endl; //checking if it oscillates in the bottom
        }
    alpha_best = alpha.row(lol);

    afile.close();
    afile2.close();
    return alpha_best;
}
