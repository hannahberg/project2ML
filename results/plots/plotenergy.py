import numpy as np
import matplotlib.pyplot as plt
import glob
import pandas as pd

df = pd.read_csv('converge_N2_D2.dat',names=['energy'])
length = len(df.energy)
time = np.linspace(0,length-1,length)
plt.plot(df['energy'])
plt.xlabel('time')
plt.ylabel('energy')
plt.show()
