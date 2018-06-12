from numpy import *
import matplotlib.pyplot as plt

files_ = ['brute_N1_d1gam0.010000_H1_rho0.100000.dat','brute_N1_d1gam0.010000_H1_rho0.300000.dat','brute_N1_d1gam0.010000_H1_rho0.600000.dat','brute_N1_d1gam0.010000_H1_rho0.800000.dat','brute_N1_d1gam0.010000_H1_rho1.000000.dat']

dtype1 = [('E', 'f8'), ('V','f8'), ('T', 'f8'), ('A', 'f8')]
data_1 = loadtxt(files_[0],dtype=dtype1)
data_2 = loadtxt(files_[1], dtype=dtype1)
data_3 = loadtxt(files_[2], dtype=dtype1)
data_4 = loadtxt(files_[3], dtype=dtype1)
data_5 = loadtxt(files_[4], dtype=dtype1)

E_1 = data_1['E']
E_2 = data_2['E']
E_3 = data_3['E']
E_4 = data_4['E']
E_5 = data_5['E']

length = len(E_1)
gdc = linspace(0,length-1,length)

plt.plot(gdc,E_1)
plt.hold('on')
plt.plot(gdc,E_2)
plt.plot(gdc,E_3)
plt.plot(gdc,E_4)
plt.plot(gdc,E_5)

plt.xlabel('GDC')
plt.ylabel('Energy')
plt.legend([r'$\rho$ = 0.1',r'$\rho$ = 0.3',r'$\rho$ = 0.6',r'$\rho$ = 0.8',r'$\rho$ = 1.0',])
plt.show()
