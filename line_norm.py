#!/usr/bin/env python3

# Use this script to generate the average lnL vs Fraction Incorrect plot
# CVRSharma 2024-04-04

import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.stats as sc
import sys
import argparse

parser = argparse.ArgumentParser(description = 'Use to quickly plot a line')

# Arguments and help options

parser.add_argument("-x", "--xdata", dest = "xdata", help = "Data for x-axis")
parser.add_argument("-y", "--ydata", dest = "ydata", help = "Data for y-axis")
parser.add_argument("-xl", "--xlab", dest = "xlab", help = "X-label")
parser.add_argument("-yl", "--ylab", dest = "ylab", help = "Y-label")
parser.add_argument("-t", "--title", dest = "title", help = "Title")
parser.add_argument("-o", "--output", dest = "output", help = "Output file name")
parser.add_argument("-n", "--norm", default=1, type=int,dest="norm", help= "Optional: Provide a normalization constant")
parser.add_argument("--flipit", action="store_true",dest="flipit", help= "Optional: Plot 1-y instead of y")

# Print --help message if no argument flags are provided by user
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

pp=open(args.xdata, 'r')
fc=open(args.ydata, 'r')

fc_lines = fc.readlines()
fc.close()

pp_lines = pp.readlines()
pp.close()

#Save data as numpy arrays to be used for plotting
frac=[]
for line in fc_lines:
    data = line.split()
    if len(data) != 1:
        frac.append(float(data[1]))
    else:
        frac.append(float(line))
frac=np.asarray(frac)

prob=[]
for line in pp_lines:
    data = line.split()
    if len(data) != 1:
        prob.append(float(data[1]))
    else:
        prob.append(float(line))
prob=np.asarray(prob)

#Apply normalization
if args.norm == 1:
    pass
else:
    frac=frac/args.norm
    prob=prob/args.norm
#Flip-it
if args.flipit is True:
    frac = 1.0 - frac
else:
    pass

#Calculate line of best fit
slope,yinter,r2,p,std_err = sc.linregress(prob,frac)
fit=slope*prob + yinter

"""
#Plot
plt.figure
labels = ['Fit','BAli-Phy',  'CLUSTALO',  'FFTNS1',  'FSA',  'LINSI',  'MUSCLE',  'PAGAN2',  'PRANK',  'PROBCONS',  'TCOFFEE']
c = ['grey','lightcoral','peru','darkorange','olive','deepskyblue','orchid','firebrick','lightsteelblue','yellowgreen']
m = ['o','v','^','p','D','*','P','>','X','<']
#plt.scatter(prob,frac,color='deepskyblue')
#plt.plot(prob,fit,color='orange')
for i in range(len(c)):
    plt.scatter(prob[i], frac[i], color=c[i], marker=m[i], s=100)
plt.plot(prob,fit,color='black')
plt.title(args.title)
plt.title(args.title)
plt.xlabel(args.xlab)
plt.ylabel(args.ylab)
#plt.xlim(0,1)
#plt.ylim(0,1)
plt.text(np.min(prob),0.8639, "$R^2$ = "+str(round(r2,3))+"\nslope = "+str(round(slope,3)), fontsize = 12,horizontalalignment='left',verticalalignment='bottom')
plt.legend(labels,loc='upper left', bbox_to_anchor=(1,1))
plt.show()
plt.savefig(args.output,bbox_inches="tight")

print("slope = "+str(slope))
"""
print("p-value = "+str(p))