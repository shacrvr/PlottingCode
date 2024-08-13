#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.stats as sc
import sys
import argparse

parser = argparse.ArgumentParser(description = 'Use to quickly plot a bar graph')

# Arguments and help options

parser.add_argument("-d", "--data", dest = "data", help = "Data")
parser.add_argument("-yl", "--ylab", dest = "ylab", help = "Y-label")
parser.add_argument("-t", "--title", dest = "title", help = "Title")
parser.add_argument("-o", "--output", dest = "output", help = "Output file name")
parser.add_argument("-s", "--scatter", action="store_true",help="Optional: User scatter instead of hist")
parser.add_argument("-n", "--norm", dest = "norm", default=0, help = "Optional: Normalize data")
# Print --help message if no argument flags are provided by user
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

data=open(args.data, 'r')
data_lines = data.readlines()
data.close()
#Save data as numpy arrays to be used for plotting
dat=[]
for line in data_lines:
    col = line.split()
    if len(col) != 1:
        dat.append(float(col[1]))
    else:
        dat.append(float(line))
dat=np.asarray(dat)

#To plot diff
best = np.min(dat)
dat = dat - best

#Normalize data
if args.norm != 0:
    dat=dat/float(args.norm)
else:
    pass

#Plot
#pos = [2,6,10,14,18,22,26]
#labels = ['BAli-Phy','CLUSTALO','FFTNS1','LINSI','MUSCLE','PROBCONS','TCOFFEE']
#c = ['grey','lightcoral','peru','darkorange','olive','deepskyblue','orchid']
#m = ['o','v','^','p','D','*','P']
pos = [2,6,10,14,18,22]
labels = ['CLUSTALO','FFTNS1','LINSI','MUSCLE','PROBCONS','TCOFFEE']
c = ['lightcoral','peru','darkorange','olive','deepskyblue','orchid']
m = ['v','^','p','D','*','P']
plt.figure()
if args.scatter is True:
    for i in range(len(pos)):
        plt.scatter(pos[i], dat[i], color=c[i], marker=m[i], s=100)
else:
    plt.bar(pos, dat, color=c, width=2)
plt.title(args.title)
plt.ylabel(args.ylab)
plt.xticks(pos,labels)
plt.show()
plt.savefig(args.output,bbox_inches='tight',pad_inches=0.2)

if args.norm != 0:
    print(best/float(args.norm))
else:
    print(best)
