# -*- coding: utf-8 -*-
#
# Import a data file, an plot
#
# Nuno de Sousa
#
# March 2014

import numpy as np
import matplotlib.pyplot as plt

#function definition
def get_column(Matrix, col):
    result = []
    size = len(Matrix)
    for i in range (size):
        result.append(Matrix[i][col])
    return result

def get_row(Matrix, row):
    result = []
    size = len(Matrix)
    for i in range (size):
        result.append(Matrix[row][i])
    return result


with open("Au_Palik.dat") as f:
    data = f.read()



data = data.split('\n')

Matrix = [[0 for x in xrange(3)] for x in xrange(len(data))]

print "size = ", len(data)

for i in range (len(data)):
    Matrix[i][0] = data[i].split('\t')[0];
    Matrix[i][1] = data[i].split('\t')[1];
    Matrix[i][2] = data[i].split('\t')[2];

n = get_column(Matrix,1)
k = get_column(Matrix,2)


epsr = []
epsi = []
l = []
for i in range(len(data)):
    l.append(float(Matrix[i][0])*1000)
    epsr.append(float(n[i])**2 - float(k[i])**2)
    epsi.append(2*float(n[i])*float(k[i]))
    print "lambda = ", l[i], " epsr =", epsr[i]


fig = plt.figure()

ax1 = fig.add_subplot(111)
#ax2 = fig.add_subplot(211)

ax1.set_title('Plot title...', fontsize=20)
ax1.set_xlabel('your x label..', fontsize=18)
ax1.set_ylabel('your y label...', fontsize=18)

#ax1.plot(l,get_column(Matrix,1), c='r', label='n')
#ax1.plot(l,get_column(Matrix,2), c='b', label='k')


with open("eps_gold.dat") as f1:
    data1 = f1.read()



data1 = data1.split('\n')

Matrix1 = [[0 for x in xrange(3)] for x in xrange(len(data1))]

print "size = ", len(data1)

for i in range (len(data1)):
    Matrix1[i][0] = data1[i].split('\t')[0];
    Matrix1[i][1] = data1[i].split('\t')[1];
    Matrix1[i][2] = data1[i].split('\t')[2];

l1 = get_column(Matrix1,0)
epsr1 = get_column(Matrix1,1)
epsi1 = get_column(Matrix1,2)



#for i in range(len(data)):
#l1.append(float(Matrix[i][0]))
    #epsr.append(float(n[i])*float(n[i]) - float(k[i])*float(k[i]))
    #epsi.append(2*float(n[i])*float(k[i]))




ax1.plot(l,epsr, c='black', label='$\epsilon^r$ Palik')
ax1.plot(l,epsi, c='g', label='$\epsilon^i$ Palik')
ax1.plot(l1,epsr1, c='blue', label='$\epsilon^r$ eq')
ax1.plot(l1,epsi1, c='red', label='$\epsilon^r$ eq')

leg = ax1.legend()

plt.show()


