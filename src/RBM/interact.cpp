#include "interact.h"
#include "solver.h"

Interact::Interact(double s_beta,
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
{a = 0.0;binsize = 0.0032;}

mat Interact::init_pos_interact(){
    mat position = init_pos_gaus();

    mat comfort_zone = too_close(position);
    //cout << "cf" << comfort_zone << endl;
    //cout << comfort_zone << endl;
    //cout << min(min(distance_part(comfort_zone))) << endl;
    //mat zeropos = zeros(N, dim);
    //mat doesthiswork = too_close(zeropos);
    //cout << "check" << doesthiswork << endl;
    //cout << "is all above a " << distance_part(doesthiswork) << endl;
    return comfort_zone;
    //return doesthiswork;
}

mat Interact::too_close(mat Rtull){
    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(0,0.5);

    mat closeness = distance_part(Rtull);

    double sdt = sqrt(dt);
    int counter = 0;
    bool get_away_you_stink = true;
    while(get_away_you_stink == true){
        //cout << "while" << endl;
        counter += 1;
        get_away_you_stink = false;
        int o; int p; int dude;
        for(o=0;o<N;o++){
            //cout << "for o" << endl;

            for(p=o+1;p<N;p++){
                //cout << "for p" << endl;
                //cout << closeness(o,p) << endl;
                if(closeness(o,p) < a){
                    //cout << "ifclose" << endl;
                    get_away_you_stink = true;
                    for(dude=0;dude<dim;dude++){
                        //cout << "for every dim new position is made" << endl;
                        Rtull(o,dude) = gaussianRNG(genMT64)*sdt;
                        closeness(o,p) = norm(Rtull.row(o)-Rtull.row(p));
                        //cout << Rtull(o,dude)<<endl;
                    }
                }
            }
            //closeness = distance_part(Rtull); //test

        }
    }

    return Rtull;
}

double Interact::wavefunc_interact(mat R, double alpha_){
    int i; int j;
    double g = 0;
    double f = 1;
    mat newR = R;
    if(dim==3){
    newR.col(2) = beta*newR.col(2);
    }
    mat newdist = distance_part(newR); //do we need this?

    if(dim==1){
        for(i=0;i<N;i++){
            g += newR(i)*newR(i); // take Product of Pi(g(Ri)            

        }
    } else{
        for(i=0;i<N;i++){
            for(int l=i+1;l<N;l++){
                if(i!=l){
                    f*=(1 - a/newdist(i,l));
                for(j=0;j<dim;j++){
                    //if(j==2){
                    //    g += beta*beta*newR(i,j)*newR(i,j);
                    //} else{
                        g += newR(i,j)*newR(i,j);//g += dot(R.row(i),R.row(i));
                    }
                }
            }
            }
        }


    double psi = exp(-alpha_*g)*f;
    return psi;
}

double Interact::d_wavefunc_interact(const mat &R, double alpha_){
    mat newR = R;
    if(dim==3){
        newR.col(2) = newR.col(2)*beta;
    }

    int i; int j;
    double g = 0;
    if(dim==1){
        for(i=0;i<N;i++){
            g += newR(i)*newR(i); // take Product of Pi(g(Ri)
        }
    } else{
        for(i=0;i<N;i++){
            for(j=0;j<dim;j++){
                g += newR(i,j)*newR(i,j);//g += dot(R.row(i),R.row(i));
            }
        }
    }
    //double f = 1; //no interaction here!!
    //double psi = exp(-alpha_*g)*f;
    return -g;
}


vec Interact::solve_interact( std::ofstream &myfile, std::ofstream &myfile5, double alphanow){ // make him take alpha as an input
    //myfile << endl << "# Calculation with interaction: " << endl;
    //myfile << "# dim = " << dim << ", N = " << N << ", dt = " << dt << ", alpha = " << alpha << " and mc = " << mc << endl << endl;

    start=clock();
    double D = 0.5; //diffusion coefficient
    double Ddt = D*dt;
    //double Ddt05 = Ddt*0.5;

    static random_device rd;
    static mt19937_64 genMT64(rd());
    static normal_distribution<double> gaussianRNG(0.,0.5);
    static uniform_real_distribution<double> doubleRNG(0,1);

    // loop over alpha when we try out
    int num_alpha = 0;
    vec alpha_vec = ones(1);
    double current_alpha = alphanow;
    double sumKE = 0;
    double sum_d_wf = 0;
    double sumKEsq = 0;
    double sum_d_wf_E = 0;
    double sdt = sqrt(dt);
    double maxR = 3.0;
    int n_bins = maxR/binsize;
    vec histbins = linspace(0,maxR,n_bins);
    vec volumes = zeros(n_bins);
    vec hists = zeros(n_bins);

    volumes(0) = (4/3.)*M_PI*pow(binsize,3);
    double sumvol = volumes(0);
    for(int bin=1;bin<n_bins;bin++){
        volumes(bin) = (4/3.)*M_PI*pow((binsize * (bin+1)),3) - sumvol;
        //volumes(bin) = 1;
        sumvol += volumes(bin);
    }

    double VOLUM = (4/3.)*datum::pi*pow(maxR,3);
    cout << sumvol << endl;
    cout << VOLUM << endl;



    vec hist_total = zeros(n_bins);
    //uvec hist;
    vec hist_TOT = zeros(n_bins);

    //current_alpha = alpha;//alpha_vec(num_alpha);
    // initialize random positions
    mat R4 = init_pos_interact();
    mat R4new = R4;
    mat distancematrix = distance_part(R4);
    //cout << distancematrix << endl;
    int i; int j; int q;
    mat Fq = quantumF(R4, current_alpha,distancematrix);
    //initialize expectation values
    double accept = 0;
    mat Fqnew = Fq;
    double greensnew;
    double greensold = 0;
    mat distR4new = distancematrix;
    for(i=0;i<mc;i++){
        //cout << i << endl;

        //propose a new position Rnew(boson_j) by moving one boson from position R(boson_j) one at the time
        for(j=0;j<N;j++){

            greensnew = 0;
            for(q=0;q<dim;q++){

                R4new(j,q) = R4(j,q) + Ddt*Fq(j,q) + gaussianRNG(genMT64)*sdt;
                //cout << R4(j,q) << endl;
                //cout << R4new(j,q) << endl;
                //Fqnew(j,q) = quantumF(R4, current_alpha,distR4new);//replace
                //cout << j << q << endl;
            }


            distR4new = distance_part(R4new);
            Fqnew = quantumF(R4new, current_alpha,distR4new);// REPLACE THIS!!!!!

            greensnew = norm(R4.row(j) - R4new.row(j) - Ddt*Fqnew.row(j));
            greensnew *= greensnew;
            greensold =norm(R4new.row(j) - R4.row(j) - Ddt*Fq.row(j));
            greensold *= greensold;
            //greens = dot(0.5*(Fq.row(j) + Fqnew.row(j)),(Ddt05*(Fq.row(j)-Fqnew.row(j))-R4new.row(j)+R4.row(j)));

            /*for(p = 0; p < N; p++){

                if(p!=j){
                    distR4new(j,p) = norm(R4new.row(j)- R4new.row(p));
                    distR4new(p,j) = distR4new(j,p);

                }
            }*/

            double greens = exp((greensold-greensnew)/(4*Ddt));

            double A = (wavefunc_interact(R4new,current_alpha))/wavefunc_interact(R4,current_alpha);
            A = abs(A); //not strictly necessary
            A *= A;
            A = greens*A;
            // test if new position is more probable than random number between 0 and 1.
            if((A > 1) || (A > doubleRNG(genMT64))){
                R4(j) = R4new(j); //accept new position
                Fq(j) = Fqnew(j);
                accept += 1;
                distancematrix = distR4new;
            }else {
                R4new(j) = R4(j);
                Fqnew(j) = Fq(j);
                distR4new = distancematrix;
            }
            // calculate change in energy
            double deltakinE = energy_interact(R4, current_alpha);

            sumKE += deltakinE;
            sumKEsq += deltakinE*deltakinE;
            myfile5 << scientific << deltakinE << endl;
            double dwf = d_wavefunc_interact(R4new,current_alpha);
            sum_d_wf += dwf;
            sum_d_wf_E += dwf*deltakinE;

            //Calculating the One-Body Density Matrix
  /*
            vec origodist = dist_origo(R4);

            for(int h=1;h<n_bins;h++){
                for(int l=0; l<N; l++) {
                    if (histbins(h-1) <= origodist(l) < histbins(h)) {
                        hist_TOT(h-1) += 1/volumes(h-1);
                    }

                }
            }
            */
            }

    }

    vec meanhist = hist_TOT/(mc*N);
    //myfile << meanhist << endl;
    num_alpha += 1;
    double mean_KE = sumKE/(mc*N);
    double mean_d_wf = sum_d_wf/(N*mc);
    double mean_d_wf_E = sum_d_wf_E/(N*mc);

    //myfile << "Energy squared = "<< sumKEsq/(mc) << endl;
    //myfile << "Variance = " << sumKEsq/(mc) - (mean_KE*mean_KE)<< endl;
    end=clock();

    cout << "Interaction and all are finished! Yay." << endl;
    vec mean_values = zeros(3);
    mean_values(0) = mean_KE;
    mean_values(1) = mean_d_wf;
    mean_values(2) = mean_d_wf_E;
    double var = sumKEsq/(mc*N) - mean_KE*mean_KE;
    myfile << "# Energy" << "     " << "Acceptance" << "   " << "CPU time" << "        " << "Solver" << endl;
    myfile << scientific << mean_KE << " " << scientific << accept/(mc*N) << " var: " << var << " " << scientific << ((double)end-(double)start)/CLOCKS_PER_SEC << "    " << 3 << "  # Interaction" << endl;
    return mean_values;
}

mat Interact::quantumF(const mat &R, double alpha_, mat rij){
    mat ngg = nablaphi(R,alpha_);
    mat nff = newnablaf(rij,R);
    return 2*ngg + 2*nff;
}

mat Interact::lapphi(const mat &R, double alpha_){
    mat lphi = zeros(N);
    mat newR = R;
    if(dim==3){
        newR.col(2) = R.col(2)*beta;
    }
    //double energy_r = 0;
    /*for(int i = 0; i < N; i++){
        for(int j = 0; j < dim; j++){
            energy_r += R(i,j)*R(i,j);
        }
    }
    lphi = -2*alpha_*((dim-1)+beta-2*alpha_*energy_r);
    return lphi;*/


    for(int k = 0; k < N; k++){
        if(dim == 3){
            lphi(k) = -2*alpha_*(dim-1 + beta-2*alpha_*dot(newR.row(k),newR.row(k))); // write more effecient, calculate 2alpha and beta^2 as own variables
        } else{
            lphi(k) = -2*alpha_*((dim-1) -2*alpha_*dot(newR.row(k),newR.row(k)));
        }
    }

    return lphi;
}

mat Interact::nablaphi(const mat &R, double alpha_){
    mat newR = R;
    if(dim==3){
        newR.col(2) = R.col(2)*beta;
    }

    return -2*alpha_*newR;

}

mat Interact::nablaphinablaF(mat R, mat distR, double alpha_){
    mat nphi = -2*alpha_*R;
    if(dim == 3){
        nphi = nphi(2)*beta;
    }
    double rkj = 0;
    mat sum = zeros(N);
    int k; int j;
    mat rk;
    mat nphik;
    for(k = 0; k < N; k++){
        nphik = nphi(k);
        rk = R(k);
        for(j = 0; j < N; j++){
            if(k != j){
                rkj = distR(k,j);
                sum(k) += dot(nphik,rk-R(j))*a/(rkj*rkj*(rkj - a));
            }
        }
    }
    return sum;
}

mat Interact::nablaf(mat R, mat distR){
    double rkj = 0;
    mat sum = zeros(N,dim);
    int k; int j;
    mat rk;
    for(k = 0; k < N; k++){
        rk = R.row(k);
        for(j = 0; j < N; j++){
            if(k != j){
                rkj = distR(k,j);
                sum.row(k) += (rk-R.row(j))*a/(rkj*rkj*(rkj - a));
            }
        }
    }
    return sum;
}

mat Interact::suma2(const mat &distanceR){
    mat summm = zeros(N);
    double rkj;
    for(int k = 0; k < N; k++){
        for(int j = 0; j < N; j++){
            if(k != j){
                rkj = distanceR(k,j);
                summm(k) = (dim - 1)*a/(rkj*rkj*(rkj - a)) + (a*a - 2*a*rkj)/(rkj*rkj*(rkj - a)*(rkj - a));
//                summm(k) = ((dim - 1)*a*rkj - (dim - 1)*a*a + a*a - 2*a*rkj)/(rkj*rkj*(rkj-a)*(rkj-a));
            }
        }
    }
    return summm;
}

mat Interact::newnablaf(const mat& init_distance, const mat &R){
    mat s = zeros(N,dim);
    double rkj = 0;
    double rkja = 0;
    mat vecrkj;
    mat vecrmi;
    //mat sumtot = zeros(N);
    for(int k = 0; k < N; k++){
        for(int j = 0; j < N; j++){
            if(k != j){
                rkj = init_distance(k,j);
                rkja = (rkj - a);
                vecrkj = R.row(k) - R.row(j);
                s.row(k) = vecrkj*a/(rkj*rkj*rkja);
            }
        }
    }
    return s;
}

mat Interact::nablafsquared(const mat& init_distance, const mat &R){
    mat s1 = zeros(N,dim);
    mat s2 = zeros(N,dim);
    double rkj = 0; double rmi = 0;
    double rkja = 0; double rmia = 0;

    mat vecrkj;
    mat vecrmi;
    mat sumtot = zeros(N);
    for(int k = 0; k < N; k++){
        for(int j = 0; j < N; j++){
            if(k != j){
                rkj = init_distance(k,j);
                rkja = (rkj - a);
                vecrkj = R.row(k) - R.row(j);
                s1.row(k) = vecrkj*a/(rkj*rkj*rkja);
            }
        }
    }
    for(int m = 0; m < N; m++){
        for(int i = 0; i < N; i++){
            if(m != i){
                rmi = init_distance(m,i);
                rmia = (rmi - a);
                vecrmi = R.row(m) - R.row(i);
                s2.row(m) = vecrmi*a/(rmi*rmi*rmia);
            }
        }
    }
    for(int p = 0; p < N; p++){
        sumtot(p) = dot(s1.row(p), s2.row(p));
        //REQUIRE(s1 == runsuma)
    }

    return sumtot;
}
/*
mat Interact::suma2(mat distanceR){
    int k; int j;
    mat suma = zeros(N);
    double rkj = 0;
    double rkja2 = 0;
    double a2 = a*a;
    for(k = 0; k < N; k++){
        for(j =0; j < N; j++){
            if(k != j){
                rkj = distanceR(k,j);
                rkja2 = (rkj - a)*(rkj - a);
                suma(k) -= a2/(rkj*rkj*rkja2);
            }
        }
    }

    return suma;
}
*/

double Interact::energy_interact(const mat& R, double alphanow){
    mat energy = zeros(N);
    double r2 = 0;
    double Vext = 0;
    mat rkj = distance_part(R);
    mat lphi = lapphi(R, alphanow);
    mat nphi = nablaphi(R, alphanow);
    mat nf = newnablaf(rkj,R);//nablaf(R,rkj);
    mat suma = suma2(rkj);
    mat nablaf2 = nablafsquared(rkj, R);
    double totsum = 0;
    double lphisum = 0;
    double nfsum = 0;
    double sumasum = 0;
    double nphisum = 0;
    double energysq = 0;
    double sumnf2 = 0;
    double sumnf2test = 0;
    for(int k = 0; k < N; k++){
        lphisum += lphi(k);
        nphisum += 2*dot(nphi.row(k),nf.row(k));
        nfsum += dot(nf.row(k),nf.row(k));
        sumasum += suma(k);
        sumnf2 += nablaf2(k);
        sumnf2test += dot(nf.row(k), nf.row(k));
        energy(k) = -0.5*(lphi(k) + 2*dot(nphi.row(k),nf.row(k)) + nablaf2(k) + suma(k));
        totsum += energy(k);
        energysq += energy(k)*energy(k);
    }

    //cout << "suma2 = " <<sumasum<< endl;
    //cout << "nablaf = " << nfsum << endl;
    //cout << "nphi = " << nphisum << endl;
    //cout << "lphi = " << lphisum << endl;
    //cout << "difference nf^2 = " << abs(sumnf2 - sumnf2test) << endl;

    double c = 0.5;
    for(int j = 0; j < N; j++) {
        r2 = 0;
        for(int q = 0; q < dim; q++) {
            r2 += R(j,q)*R(j,q);
        }
        Vext += c*r2; //calculate potential energy
    }
    //cout << "Vext = " << Vext << endl;
    return totsum + Vext /*+ lphi*/;
}

double Interact::best_alpha(){
    ofstream alphafile;ofstream alphafile2;
    alphafile.open("onlyg.dat");
    alphafile2.open("alphas_start_0.5.dat");
    double alpha_the_best;
    vec mean_values;
    vec alpha_ = zeros(1000);
    alpha_(0) = 0.5;
    double gamma = 0.02;
    double tol = 0.00001;
    for(int i=0;i<999;i++){
        mean_values = solve_interact(alphafile, alphafile2, alpha_(i));
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
    alphafile << alpha_the_best << endl;
    alphafile.close();
    return alpha_the_best;
}

vec Interact::dist_origo(const mat& R){
    vec d_origo = zeros(N);
    for(int i=0; i<N; i++){
        d_origo(i) = norm(R.row(i));
    }
    return d_origo;
}
