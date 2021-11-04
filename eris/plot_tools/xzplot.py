import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.mlab import csv2rec
import cmath
import matplotlib.mlab as ml
import sys

print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)
print 'Test = ', str(sys.argv[1])



filename = 'Epartfields.dat'

edge = 15;

data = csv2rec('Epartfields.dat', delimiter = '\t')

#print data

#col1, col2, col3, col4, col5 = np.loadtxt(filename,skiprows=1, usecols=(5), unpack=True)
l = data['l']
i = data['i']
x = data['x']
y = data['y']
z = data['z']
rEx = data['rex']
iEx = data['iex']
rEy = data['rey']
iEy = data['iey']
rEz = data['rez']
iEz = data['iez']
#print "x = ", x
#print "rEy = ", rEy



npsec = 5

absE = [[0 for var in range(4)] for var in range(len(l))]

for i in range(len(l)):
    Ef = np.array([rEx[i]+1j*iEx[i],rEy[i]+1j*iEy[i] ,rEz[i]+1j*iEz[i]])
    Efconj = np.conjugate(Ef)
    absE[i][0] = round((x[i]/230.)*15.,0)
    absE[i][1] = round((y[i]/230.)*15.,0)
    absE[i][2] = round((z[i]/230.)*15.,0)
    absE[i][3] = np.inner(Ef,Efconj)
    #xcomp = complex(rEx,iEx)
    #ycomp = complex(rEy,iEy)
    #Ef=np.array[xcomp,ycomp]

#print "absE = ", absE


#f = open('workfile.dat', 'w')
#for i in range(len(l)):
#    f.write("%s\t%s\t%s\t%s\n" %(absE[i][0], absE[i][1], absE[i][2], absE[i][3]))
#f.closed

sqtable = [[0 for var in range(4)] for var in range(edge**3)]

cum = 0

for i in range(edge):
    for j in range(edge):
        for k in range(edge):
            sqtable[cum][0] = 2*i-14
            sqtable[cum][1] = 2*j-14
            sqtable[cum][2] = 2*k-14
            cum = cum + 1


#print "sq = ", sqtable

#print "leng = ", len(absE)

#for i in range(len(sqtable)):
#    print [sq]

for i in range(len(sqtable)):
    for j in range(0,len(absE)):
        #print "i = ", i , "   j = ", j
        #print "abs = ", absE[j][0], "  sqtable = ", sqtable[i][0]
        if(absE[j][0] == sqtable[i][0]):
            if(absE[j][1] == sqtable[i][1]):
                if(absE[j][2] == sqtable[i][2]):
                    sqtable[i][3] = absE[j][3]

#print "sq = ", len(sqtable)


#f = open('workfile2.dat', 'w')
#for i in range(len(sqtable)):
#    f.write("%s\t%s\t%s\t%s\n" %(sqtable[i][0], sqtable[i][1], sqtable[i][2], sqtable[i][3]))
#f.closed


#pass to a square table



Erep_a = []
Erep_b = []
Erep_int = []


#xz data
for i in range(len(sqtable)):
    if(sqtable[i][1] == 0.):
        Erep_a.append(sqtable[i][0])
        Erep_b.append(sqtable[i][2])
        Erep_int.append(sqtable[i][3])

#f = open('workfile3.dat', 'w')
#for i in range(len(Erep_int)):
#    f.write("%s\t%s\t%s\t%s\n" %(i, Erep_a[i], Erep_b[i], Erep_int[i]))
#f.closed


Erep_a = np.r_[Erep_a,min(Erep_a),max(Erep_a)]
Erep_b = np.r_[Erep_b,min(Erep_b),max(Erep_b)]
Erep_int = np.r_[Erep_int,min(Erep_int),max(Erep_int)]
xi = np.linspace(min(Erep_a), max(Erep_a), max(Erep_a))
yi = np.linspace(min(Erep_b), max(Erep_b), max(Erep_b))
zi = ml.griddata(Erep_a, Erep_b, Erep_int, xi, yi)
#print "zi = ",zi

#print "xi = ", xi

alternative_matrix = [[0 for var in range(edge)] for var in range(edge)]
beta = 0
alpha = 0

#print "leng Erep = ", len(Erep_int), " ", len(Erep_a)

for i in range(len(Erep_int)-2):
    #print "i = ", i, "alpha = ", alpha, " beta = ", beta
    alternative_matrix[beta][alpha] = Erep_int[i].real
    alpha = alpha + 1
    #print alternative_matrix
    if(alpha == edge):
        alpha = 0
        beta = beta + 1

#print alternative_matrix



plt.imshow(alternative_matrix,interpolation='none')
plt.xlabel('Z')
plt.ylabel('X')
plt.grid(True)
plt.title("$\lambda$ = %s (nm)"%(l[1]))
#plt.show()

#plt.contour(xi, yi, zi, 30, linewidths = 0.5, colors = 'k')
#plt.pcolormesh(xi, yi, zi, cmap = plt.get_cmap('hot'))

#plt.colorbar()
#plt.scatter(Erep_a, Erep_b, marker = 'o', c = 'b', s = 10, zorder = 10)
#plt.xlim(min(Erep_a), max(Erep_a))
#plt.ylim(min(Erep_b), max(Erep_b))

#plt.ylabel(r'$\ln\left(\frac{x_a-x_b}{x_a-x_c}\right)$')
#plt.xlabel(r'x dimension', fontsize=20)
#plt.ylabel(r'y dimension', fontsize=20)
#plt.tick_params(labelsize=20)
#plt.xticks([0.4,0.14,0.2,0.2], fontsize = 50)

#Modification of the name

id_name = str(sys.argv[1])

str = "_";
seq = ("plotxz", id_name); # This is sequence of strings.
print str.join( seq );

seq = str.join(seq);

str2= ".";
seq2= (seq,"png");
final_name = str2.join(seq2);

plt.savefig(final_name)
#plt.show()




#3D plot instructions

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(data['x'],data['y'], data['z'])

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z Label')

#plt.show()