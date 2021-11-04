import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.mlab import csv2rec

data = csv2rec('Au_Palik.dat', delimiter = '\t')



fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_title('Plot title...', fontsize=20)
ax1.set_xlabel('your x label..', fontsize=18)
ax1.set_ylabel('your y label...', fontsize=18)

ax1.plot(data['l'],data['n'], c='r', label='n')
ax1.plot(data['l'],data['k'], c='b', label='k')

l = data['l']
n = data['n']
k = data['k']
#epsr = []
#epsi = []

#epsr = n^2-k^2
#epsi = 2*n*k

print " l = ", l

#ax1.plot(data['l'],data['er'], c='r', label='$\epsilon^{r}$')
#ax1.plot(data['l'],data['ei'], c='b', label='$\epsilon^{i}$')

leg = ax1.legend()

plt.show()
