#include <gtest/gtest.h>
#include <gmock/gmock-matchers.h>
#include "/home/hannahcb/project2ML/src/RBM/solver.h"
using namespace testing;

TEST(interaction, project2ml)
{
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 2;
    double hbar = 1;
    double mass = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = true;
    Solver S(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM,interactionswitch); // initialize Solver class

    vec Rones = {0,0,-1,-1};
    double real_dist;
    vec R;
    double dist = 0;
    for(int i = 1; i<1000; i++){
        R = Rones*i;
        real_dist = 1/sqrt(2*i*i);
        dist = S.calc_interaction(R);
        EXPECT_NEAR(real_dist, dist, 1e-5);
    }
}

TEST(energy, zerobias){
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 1;
    double hbar = 1;
    double mass = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;
    Solver S1(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM,interactionswitch);
    vec a = zeros(numM);
    vec b = zeros(hidden);
    mat w = zeros(numM,hidden);
    vec X = zeros(numM);
    double energy = S1.E_L(a,b,w,X);
    double energy_al = S1.energy_analytic();
    ASSERT_FLOAT_EQ(energy,energy_al);
}
TEST(energy,onebias){
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 1;
    double hbar = 1;
    double mass = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;
    Solver S1(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM,interactionswitch);
    vec a = ones(numM);
    vec b = ones(hidden);
    mat w = ones(numM,hidden);
    vec X = ones(numM);
}

TEST(wavefunc,zerobias){
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 1;
    double hbar = 1;
    double mass = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;
    Solver S1(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM,interactionswitch);
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
    ASSERT_FLOAT_EQ(wavefunczero,calcwavefunc);
}

TEST(wavefunc, onebias){
    double rho = 0.01;
    int numpart = 2;
    int mc = 100000;
    int howmanyDs = 1;
    double hbar = 1;
    double mass = 1;
    double omega = 1;
    double sig = 1;
    double numM = howmanyDs*numpart;
    int hidden = 2;
    double dt = 0.01;
    bool interactionswitch = false;
    Solver S1(hbar, mass, omega, rho, mc, numpart, howmanyDs, dt, sig, hidden, numM,interactionswitch);
    vec a = ones(numM);
    vec b = ones(hidden);
    mat w = ones(numM,hidden);
    vec X = ones(numM);
    double sum1 = 0;
    for(int i = 0; i < numM; i++){
        sum1 += (X(i) - a(i))*(X(i)-a(i))/(2*sig*sig);
    }
    double sum2 = 0;
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
    EXPECT_NEAR(wavefunc_an,wavefunc_prog, 1e-5);


}
