# -*- coding: utf-8 -*-
#
# Import a data file, an plot
#
# Nuno de Sousa
#
# March 2014

import numpy as np
import matplotlib.pyplot as plt
import math
import os

#function definition
def get_column(Matrix, col):
    result = []
    size = len(Matrix)
    for i in range (size):
        result.append(float(Matrix[i][col]))
    return result

def get_row(Matrix, row):
    result = []
    size = len(Matrix[0])
    for i in range (size):
        result.append(float(Matrix[row][i]))
    return result

def charge_Matrix(name, n_col):
    with open(name) as f:
        data = f.read()

    data = data.split('\n')

    Matrix = [[0 for x in xrange(n_col)] for x in xrange(len(data)-2)]

    for i in range (len(data)-1):
        if(i!=0):
            for j in range(n_col):
                Matrix[i-1][j] = data[i].split('\t')[j];
    return Matrix

def sec_ext(prefactor, xa, xb):
    sec = []
    for i in range(len(prefactor)):
        sec.append(prefactor[i]*(xa[i]+xb[i]))
    return sec

def sec_exta(prefactor, xa):
    sec = []
    for i in range(len(prefactor)):
        sec.append(prefactor[i]*(xa[i]))
    return sec

def sec_extb(prefactor, xb):
    sec = []
    for i in range(len(prefactor)):
        sec.append(prefactor[i]*(xb[i]))
    return sec

def sec_scat(prefactor, ra, ia, rb, ib):
    sec = []
    for i in range(len(prefactor)):
        sec.append(prefactor[i]*((ra[i]**2+ia[i]**2)+(rb[i]**2+ib[i]**2)))
    return sec

#this remove blank lines
os.system("grep -v '^$' outputprojectionaeo.dat > projectionaeo.dat")
os.system("grep -v '^$' outputprojectionbeo.dat > projectionbeo.dat")

a = charge_Matrix("output_projectionaeo_xpol.dat", 7)
b = charge_Matrix("output_projectionbeo_xpol.dat", 7)
dir = charge_Matrix("output_sections_xpol.dat", 5) #resultudo directo do DDA

#print ae11


lambdaae11 = get_column(a,0)
reae11 = get_column(a,3)
imae11 = get_column(a,4)
rebo11 = get_column(b,5)
imbo11 = get_column(b,6)

prefactor = [3*x**2/(2*math.pi) for x in lambdaae11]

sigma_ext = sec_ext(prefactor,reae11,rebo11)
sigma_scat = sec_scat(prefactor, reae11, imae11, rebo11, imbo11)
sigma_exta = sec_exta(prefactor, reae11)
sigma_extb = sec_exta(prefactor, rebo11)

fig = plt.figure()

ax1 = fig.add_subplot(111)

#axis configurations
#ax1.set_title('Plot title...', fontsize=20)
ax1.set_xlabel('$\lambda$ (nm)', fontsize=22)
ax1.set_ylabel('$\sigma_{ext} (nm^2)$', fontsize=22)
ax1.tick_params(labelsize=20)
#use the scientific notation for the y representation
ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')

#define the size of the legends
plt.rcParams.update({'font.size': 16, 'font.family': 'serif'})

ax1.plot(get_column(dir,0),get_column(dir,3), linewidth=2, c='r', label='$DDA_{ext}$')
#ax1.plot(get_column(a,0),sigma_scat, c='b', label='$\sigma_{scat}$')
ax1.plot(get_column(a,0),sigma_exta, linewidth=2, c='g', label='$\sigma_{ext}(a_1)$')
ax1.plot(get_column(a,0),sigma_extb, linewidth=2, c='black', label='$\sigma_{ext}(b_1)$')
ax1.plot(get_column(a,0),sigma_ext, linewidth=2, c='b', label='$\sigma_{ext}(a_1+b_1)$')


leg = ax1.legend()

fig.set_size_inches(13, 8.5)
fig.savefig('sections.png', dpi=100)

plt.show()

os.system("rm projectionaeo.dat projectionbeo.dat")


