from numpy import *
import matplotlib.pyplot as plt

files_ = ['brute_N2_d2gam0.100000_H2_rho1.000000.dat','imp_N2_d2gam1.000000_H2_dt0.100000.dat','gibbs_N2_d2gam0.100000_H2_sig0.700000.dat']

dtype1 = [('E', 'f8'), ('V','f8'), ('T', 'f8'), ('A', 'f8')]
dtype2 = [('E', 'f8'), ('V','f8'), ('T', 'f8')]
data_1 = loadtxt(files_[0],dtype=dtype1)
data_2 = loadtxt(files_[1], dtype=dtype1)
data_3 = loadtxt(files_[2], dtype=dtype2)
#data_4 = loadtxt(files_[3], dtype=dtype1)
#data_5 = loadtxt(files_[4], dtype=dtype1)
#data_6 = loadtxt(files_[5], dtype=dtype1)

E_1 = data_1['E']
E_2 = data_2['E']
E_3 = data_3['E']
#E_4 = data_4['E']
#E_5 = data_5['E']
#E_6 = data_6['E']


length = len(E_1)
gdc = linspace(0,length-1,length)

plt.plot(gdc,E_1)
plt.hold('on')
plt.plot(gdc,E_2)
plt.plot(gdc,E_3)
#plt.plot(gdc,E_4)
#plt.plot(gdc,E_5)
#plt.plot(gdc,E_6)
plt.xlabel('GDC')
plt.ylabel('Energy')
plt.legend(['Brute','Imp', 'Gibbs'])
plt.show()

