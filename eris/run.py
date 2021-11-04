#!/usr/bin/python

import os

number = 12
nodes = 4
total = number*nodes

os.system("qsub -q gsnm -l nodes=4:ppn=12 script")

#string1 = "qsub -q gsnm -l nodes="
#string2 = ":ppn="
#string3 = " script "
#string = string1 + str(nodes) + string2 + str(number) + string3 + str(total)

#os.system(string)


#os.system("qsub -q gsnm -l nodes=gsnm13.icmm.csic.es:ppn=12 script")
#os.system("qsub -q gsnm -l nodes=1:ppn=12 script")
#os.system("qsub -q gsnm -l nodes=2:ppn=12 script")
#two nodes with 12 cores each