# -*- coding: utf-8 -*-
#
# Import a data file, an plot
#
# Nuno de Sousa
#
# March 2014

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import csv2rec
from mpl_toolkits.mplot3d import Axes3D

data = csv2rec('power1.dat', delimiter='\t')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist((data['x']+data['y']+data['z'])/3,bins=100, range=(1.0,1.1))

plt.show()
