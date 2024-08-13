#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.stats as sc
import sys
import argparse

parser = argparse.ArgumentParser(description = 'Use to two bar graphs on the same axes')

# Arguments and help options

parser.add_argument("-d1", "--data1", dest = "data1", help = "Data 1")
parser.add_argument("-d2", "--data2", dest = "data2", help = "Data 2")
parser.add_argument("-yl", "--ylab", dest = "ylab", help = "Y-label")
parser.add_argument("-t", "--title", dest = "title", help = "Title")
parser.add_argument("-o", "--output", dest = "output", help = "Output file name")
parser.add_argument("-n", "--norm", dest = "norm", default=0, help="Optional: Provide a normalization constant")
# Print --help message if no argument flags are provided by user
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

data1=open(args.data1, 'r')
data1_lines = data1.readlines()
data1.close()
#Save data as numpy arrays to be used for plotting
dat1=[]
for line in data1_lines:
    col = line.split()
    if len(col) != 1:
        dat1.append(float(col[1]))
    else:
        dat1.append(float(line))
dat1=np.asarray(dat1)

data2=open(args.data2, 'r')
data2_lines = data2.readlines()
data2.close()
#Save data as numpy arrays to be used for plotting
dat2=[]
for line in data2_lines:
    col = line.split()
    if len(col) != 1:
        dat2.append(float(col[1]))
    else:
        dat2.append(float(line))
dat2=np.asarray(dat2)

#Normalize data
if args.norm != 0:
    dat1=dat1/float(args.norm)
    dat2=dat2/float(args.norm)
else:
    pass

#Plot
pos = [2,6,10,14,18,22]
pos=np.asarray(pos)
labels = ['L/MDH',  'Terpene Synthase',  'Kinase',  'KaiB',  'a/b Hydrolase',  'Teleost LDH']
c=['deepskyblue','orange']
plt.figure()
plt.bar(pos - 0.55, dat1, width=1,color=c[0])
plt.bar(pos + 0.55, dat2, width=1,color=c[1])
plt.title(args.title)
plt.ylabel(args.ylab)
plt.xticks(pos,labels,rotation=90)
plt.legend(['LG+FO+G4','Poisson+FO+G4'],loc='upper left',bbox_to_anchor=(1,1))
plt.show()
plt.savefig(args.output,bbox_inches='tight',pad_inches=0.2)

