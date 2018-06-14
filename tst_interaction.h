#include <gtest/gtest.h>
#include <gmock/gmock-matchers.h>
#include "/home/hannahcb/project2ML/src/RBM/solver.h"
using namespace testing;
// test of the interaction
TEST(interaction, booltrue)
{
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 2;
    double omega = 1;
    double sig = 1;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = true;
    double spread = 0.5;
    Solver S(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);

    vec Rones = {0,0,-1,-1};
    double real_dist;
    vec R;
    double dist = 0;
    for(int i = 1; i<1000; i++){
        R = Rones*i;
        real_dist = 1/sqrt(2*i*i);
        dist = S.calc_interaction(R);
        EXPECT_NEAR(real_dist, dist, 1e-5) << "interaction not working";
    }
}

TEST(energy, zerobias){
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;
    double spread = 0.5;
    Solver S1(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
    vec a = zeros(numM);
    vec b = zeros(hidden);
    mat w = zeros(numM,hidden);
    vec X = zeros(numM);
    double energy = S1.E_L(a,b,w,X);
    double energy_al = S1.energy_analytic();
    ASSERT_FLOAT_EQ(energy,energy_al) << " energy not 2 with zero bias";
}

TEST(energy,onebias){
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;
    double spread = 0.5;
    Solver S1(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
    vec a = ones(numM);
    vec b = ones(hidden);
    mat w = ones(numM,hidden);
    vec X = ones(numM);
}

TEST(energy, Isimple){
    double rho = 0.01;
    int numpart = 3;
    int mc = 100000;
    int howmanyDs = 2;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;
    double spread = 0.5;
    Solver S1(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
    vec a = ones(numM);
    vec b = ones(hidden);
    mat w = ones(numM,hidden);
    vec X = ones(numM);
    vec xa = X-a;
    double sumxa = 0;
//    double un = 0;
    double sum;
    double totI = 0;
    for(int m = 0; m < numM; m++){
        double sumN = 0;
        sumxa = -xa(m)/(sig*sig);
        for(int n = 0; n < hidden; n++){
            sum = 0;
            for(int i = 0; i < numM; i++){
                sum += X(i)*w(i,n)/(sig*sig);
            }
//            un = -b(n) - sum;
            sumN += w(m,n)/(exp(-S1.u(b(n),X, w.col(n))) + 1);
        }
        sumxa += (1/(sig*sig))*sumN;
        totI += sumxa*sumxa;
        cout << totI << endl;
    }
//    double totI = sumxa + (1/(sig*sig))*sumN;
//    totI *= totI;
    double progI = S1.I(a,b,w,X);
    EXPECT_NEAR(totI,progI,1e-4) << "I not correct!";

}

TEST(energy, IIsimple){
    double rho = 0.01;
    int numpart = 3;
    int mc = 100000;
    int howmanyDs = 2;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;
    double spread = 0.5;
    Solver S1(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
    vec a = S1.init_a(spread);
    vec b = S1.init_b(spread);
    mat w = S1.init_w(spread);
    vec X = S1.init_X();
    double tempII;
    double eu = 0;
    double totII = 0;
    for(int i = 0; i< numM;i++){
        tempII = 0;
        for(int j = 0; j < hidden; j++){
            eu = exp(S1.u(b(j),X,w.col(j)));
            tempII +=(1.0/(sig*sig*sig*sig))*w(i,j)*w(i,j)*(eu/((1+eu)*(1+eu)));
        }
        totII += tempII;
    }

    double progII = S1.II(b,w,X);
    EXPECT_NEAR(totII,progII,1e-4) << "II not correct!";
}

TEST(energy,randinit){
    double rho = 0.01;
    int numpart = 3;
    int mc = 100000;
    int howmanyDs = 2;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;
    double spread = 0.5;
    Solver S1(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
    vec a = S1.init_a(spread);
    vec b = S1.init_b(spread);
    mat w = S1.init_w(spread);
    vec X = S1.init_X();
    double progE = S1.E_L(a,b,w,X);
    double X2 = 0;
    for(int m = 0; m < numM; m++){
        X2 += X(m)*X(m);
    }
    double E = 0.5*(-(S1.I(a,b,w,X) - numM/(sig*sig) + S1.II(b,w,X)) + omega*omega*X2);
    EXPECT_NEAR(E,progE,1e-5) << "energy not right";

}

// test of the wavefunction

TEST(wavefunc,zerobias){
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;
    double spread = 0.5;
    Solver S1(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
    vec a = zeros(numM);
    vec b = zeros(hidden);
    mat w = zeros(numM,hidden);
    vec X = zeros(numM);

    double wavefunczero = exp(0);
    double sum = 0;
    for(int j = 0; j < hidden; j++){
        sum += 1 + exp(0);
    }
    wavefunczero *= sum;
    double calcwavefunc = S1.wavefunc(a,b,w,X);
    ASSERT_FLOAT_EQ(wavefunczero,calcwavefunc) << "wavefunction not correct with zero bias";
}

TEST(wavefunc, onebias){
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;
    double spread = 0.5;
    Solver S1(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
    vec a = ones(numM);
    vec b = ones(hidden);
    mat w = ones(numM,hidden);
    vec X = ones(numM);
    double sum1 = 0;
    for(int i = 0; i < numM; i++){
        sum1 += (X(i) - a(i))*(X(i)-a(i))/(2*sig*sig);
    }
    double totsum = 1;
    for(int j = 0; j < hidden; j++){
        double sum3 = b(j);
        for(int k = 0; k < numM; k++){
            sum3 += X(k)*w(k,j)/(sig*sig);
        }
        totsum *= 1 + exp(sum3);
    }
    double wavefunc_an = exp(sum1)*totsum;
    double wavefunc_prog = S1.wavefunc(a,b,w,X);
    EXPECT_NEAR(wavefunc_an,wavefunc_prog, 1e-5) << "wavefunction not correct with ones in bias";


}

// test of u function
TEST(u,zeroinit){
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;
    double spread = 0.5;
    Solver S1(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
    vec a = zeros(numM);
    vec b = zeros(hidden);
    mat w = zeros(numM,hidden);
    vec X = zeros(numM);
    double progu;
    double u_an;
    double sum;
    for(int j = 0; j < hidden; j++){
        sum = 0;
        progu = S1.u(b(j),X,w.col(j));
        for(int i = 0; i < numM; i++){
            sum += X(i)*w(i,j)/(sig*sig);
        }
        u_an = b(j) + sum;
        EXPECT_NEAR(u_an,progu,1e-8) << "u not calculated correctly in zero case";
    }
}

TEST(u,oneinit){
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;    double spread = 0.5;
    Solver S1(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
    vec a = ones(numM);
    vec b = ones(hidden);
    mat w = ones(numM,hidden);
    vec X = ones(numM);
    double progu;
    double u_an;
    double sum;
    for(int j = 0; j < hidden; j++){
        sum = 0;
        progu = S1.u(b(j),X,w.col(j));
        for(int i = 0; i < numM; i++){
            sum += X(i)*w(i,j)/(sig*sig);
        }
        u_an = b(j) + sum;
        EXPECT_NEAR(u_an,progu,1e-8) << "u not calculated correctly in one case";
    }
}

TEST(u,randinit){
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 2;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 4;
    double dt = 0.01;
    bool interactionswitch = false;
    double spread = 0.5;
    Solver S1(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
    vec a = S1.init_a(spread);
    vec b = S1.init_b(spread);
    mat w = S1.init_w(spread);
    vec X = S1.init_X();
    double progu;
    double u_an;
    double sum;
    for(int j = 0; j < hidden; j++){
        sum = 0;
        progu = S1.u(b(j),X,w.col(j));
        for(int i = 0; i < numM; i++){
            sum += X(i)*w(i,j)/(sig*sig);
        }
        u_an = b(j) + sum;
        EXPECT_NEAR(u_an,progu,1e-8) << "u not calculated correctly in random case";
    }
    double E = S1.E_L(a,b,w,X);
    cout << "E "  << E << endl;
}


TEST(gibbs,wavefunc0){
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;
    double spread = 0.5;
    Solver S1(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
    vec a = zeros(numM);
    vec b = zeros(hidden);
    mat w = zeros(numM,hidden);
    vec X = zeros(numM);

    double wavefunczero = exp(0);
    double sum = 0;
    for(int j = 0; j < hidden; j++){
        sum += 1 + exp(0);
    }
    wavefunczero *= sum;
    wavefunczero = sqrt(wavefunczero);
    double calcwavefunc = S1.wavefunc_g(a,b,w,X);
    ASSERT_FLOAT_EQ(wavefunczero,calcwavefunc) << "wavefunction_Gibbs not correct with zero bias";

}

TEST(gibbs, onebias){
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;
    double spread = 0.5;
    Solver S1(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
    vec a = ones(numM);
    vec b = ones(hidden);
    mat w = ones(numM,hidden);
    vec X = ones(numM);
    double sum1 = 0;
    for(int i = 0; i < numM; i++){
        sum1 += (X(i) - a(i))*(X(i)-a(i))/(2*sig*sig);
    }
    double totsum = 1;
    for(int j = 0; j < hidden; j++){
        double sum3 = b(j);
        for(int k = 0; k < numM; k++){
            sum3 += X(k)*w(k,j)/(sig*sig);
        }
        totsum *= 1 + exp(sum3);
    }
    double wavefunc_an = sqrt(exp(sum1)*totsum);
    double wavefunc_prog = S1.wavefunc_g(a,b,w,X);
    EXPECT_NEAR(wavefunc_an,wavefunc_prog, 1e-5) << "wavefunction_Gibbs not correct with ones in bias";

}


TEST(gibbs, randbias){
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;
    double spread = 0.5;
    Solver S1(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
    vec a = S1.init_a(spread);
    vec b = S1.init_b(spread);
    mat w = S1.init_w(spread);
    vec X = S1.init_X();
    double sum1 = 0;
    for(int i = 0; i < numM; i++){
        sum1 -= (X(i) - a(i))*(X(i)-a(i))/(2*sig*sig);
    }
    double totsum = 1;
    for(int j = 0; j < hidden; j++){
        totsum *= 1 + exp(S1.u(b(j),X,w.col(j)));
//        double sum3 = b(j);
//        for(int k = 0; k < numM; k++){
//            sum3 += X(k)*w(k,j)/(sig*sig);
//        }
//        totsum *= 1 + exp(sum3);
    }

    double wavefunc_an = sqrt(exp(sum1)*totsum);
    double wavefunc_prog = S1.wavefunc_g(a,b,w,X);
    EXPECT_NEAR(wavefunc_an,wavefunc_prog, 1e-5) << "wavefunction_Gibbs not correct with random case";

}



TEST(gibbs, energy){
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 2;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;
    double spread = 0.5;
    Solver S1(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
    vec a = S1.init_a(spread);
    vec b = S1.init_b(spread);
    mat w = S1.init_w(spread);
    vec X = S1.init_X();
    double X2 = 0;
    double progELGibbs = S1.ELGibbs(a,b,w,X);
    for(int m = 0; m < numM; m++){
        X2 += X(m)*X(m);
    }
//    double Egib = 0.5*(-(0.25*S1.I(a,b,w,X)) - numM/(2*sig*sig) + 0.5*S1.II(b,w,X) + omega*omega*X2 );
    double E = 0.5*(-(0.25*S1.I(a,b,w,X) - numM/(2*sig*sig) + 0.5*S1.II(b,w,X)) + omega*omega*X2);
    EXPECT_NEAR(E,progELGibbs,1e-5) << "ELGibbs not working";
}


TEST(gibbs, energyzerobias){
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 2;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;
    double spread = 0.5;
    Solver S1(omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, interactionswitch, spread);
    vec a = zeros(numM);
    vec b = zeros(hidden);
    mat w = zeros(numM,hidden);
    vec X = zeros(numM);
    double X2 = 0;
    double progELGibbs = S1.ELGibbs(a,b,w,X);
    for(int m = 0; m < numM; m++){
        X2 += X(m)*X(m);
    }
//    double Egib = 0.5*(-(0.25*S1.I(a,b,w,X)) - numM/(2*sig*sig) + 0.5*S1.II(b,w,X) + omega*omega*X2 );
    double E = 0.5*(-(0.25*S1.I(a,b,w,X) - numM/(2*sig*sig) + 0.5*S1.II(b,w,X)) + omega*omega*X2);
    cout << E << endl;
    EXPECT_NEAR(E,progELGibbs,1e-5) << "ELGibbs not working, zerobias";
}
