import numpy as np
import matplotlib.pyplot as plt

filename='output_force_xpol.dat'
n_part = 1189
n_lam = 100

lamda = np.genfromtxt(filename,usecols=(0),delimiter='\t',dtype=None)
col1 = np.genfromtxt(filename,usecols=(1),delimiter='\t',dtype=None)
col2 = np.genfromtxt(filename,usecols=(5),delimiter='\t',dtype=None)
col3 = np.genfromtxt(filename,usecols=(6),delimiter='\t',dtype=None)
col4 = np.genfromtxt(filename,usecols=(7),delimiter='\t',dtype=None)

lam = []
fx = []
fy = []
fz = []

sumfx = 0


for i in range(n_lam):
    sumfx = 0
    sumfy = 0
    sumfz = 0
    lam.append(lamda[n_part*i])
    for j in range(n_part):
        sumfx = sumfx + float(col2[n_part*i+j])
        sumfy = sumfy + float(col3[n_part*i+j])
        sumfz = sumfz + float(col4[n_part*i+j])
    fx.append(float(sumfx))
    fy.append(float(sumfy))
    fz.append(float(sumfz))


fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_title('Force', fontsize=20)
ax1.set_xlabel('lambda [nm]', fontsize=18)
ax1.set_ylabel('your y label...', fontsize=18)

ax1.plot(lam,fx, c='r', label='x')
ax1.plot(lam,fy, c='g', label='y')
ax1.plot(lam,fz, c='b', label='z')

leg = ax1.legend()

plt.show()


#print fx
#print fy
#print fz