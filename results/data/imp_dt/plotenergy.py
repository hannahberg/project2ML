from numpy import *
import matplotlib.pyplot as plt

files_ = ['imp_N1_d1gam0.010000_H1_dt0.001000.dat','imp_N1_d1gam0.010000_H1_dt0.005000.dat','imp_N1_d1gam0.010000_H1_dt0.010000.dat','imp_N1_d1gam0.010000_H1_dt0.100000.dat','imp_N1_d1gam0.010000_H1_dt1.000000.dat']

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
plt.legend(['dt = 0.001','dt = 0.005','dt = 0.01','dt = 0.1','dt = 1.0'])
plt.show()
