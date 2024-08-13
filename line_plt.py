#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.stats as sc
from scipy.optimize import curve_fit
import sys
import argparse

parser = argparse.ArgumentParser(description = 'Use to quickly plot a curve and perform linear, exp, or log fits of data')

# Arguments and help options

parser.add_argument("-x", "--xdata", dest = "xdata", help = "Data for x-axis")
parser.add_argument("-y", "--ydata", dest = "ydata", help = "Data for y-axis")
parser.add_argument("-xl", "--xlab", dest = "xlab", help = "X-label")
parser.add_argument("-yl", "--ylab", dest = "ylab", help = "Y-label")
parser.add_argument("-t", "--title", dest = "title", help = "Title")
parser.add_argument("-o", "--output", dest = "output", help = "Output file name")
parser.add_argument("--cum", dest = "cum",action="store_true", help = "Optional: Plot using the cumulative y-variable")
parser.add_argument("--exp", dest = "exp",action="store_true", help = "Optional: fit using an exponential function")
parser.add_argument("--log", dest = "log",action="store_true", help = "Optional: fit using an logarithmic function")


# Print --help message if no argument flags are provided by user
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

x=open(args.xdata, 'r')
y=open(args.ydata, 'r')

y_lines = y.readlines()
y.close()

x_lines = x.readlines()
x.close()

#Save data as numpy arrays to be used for plotting
yy=[]
for line in y_lines:
    data = line.split()
    if len(data) != 1:
        yy.append(float(data[1]))
    else:
        yy.append(float(line))
yy=np.asarray(yy)

xx=[]
for line in x_lines:
    data = line.split()
    if len(data) != 1:
        xx.append(float(data[1]))
    else:
        xx.append(float(line))
xx=np.asarray(xx)


# Check if cum = true
fit_check=0
if args.cum is True:
    num_bins=len(yy)
    c,b = np.histogram(yy,num_bins)
    pdf = c/sum(c)
    cdf = np.cumsum(pdf)
    yy =cdf
elif args.exp is True:
    fit_check=2
    def func(x,a,b,c):
        return a * np.exp(b * x) + c
    popt, pcov = curve_fit(func, xx, yy)
    #c = popt[0]-1.0
elif args.log is True:
    fit_check=2
    def func(x, a, b):
        return a * np.log(x) + b
    popt, pcov = curve_fit(func, xx, yy)
    #c = 0
else:
    #Calculate line of best fit
    fit_check=1
    slope,yinter,r2,p,std_err = sc.linregress(xx,yy)
    fit=slope*xx + yinter

#Plot
plt.figure
plt.scatter(xx,yy,color='orange')
if fit_check == 1:
    plt.plot(xx,fit,color='purple')
elif fit_check == 2:
    plt.plot(xx, func(xx, *popt), 'g--',
         label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
    plt.legend()
    #leg = np.append(popt,c)
    #plt.plot(xx, func(xx, *popt) - c, 'g--',
    #     label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(leg))
    #plt.legend()
plt.title(args.title)
plt.xlabel(args.xlab)
plt.ylabel(args.ylab)
#plt.xlim(0,1)
#plt.ylim(0,1)
if fit_check==1:
    plt.legend(['R^2 = {0}'.format(round(r2,3)),'slope = {0}'.format(round(slope,3))],loc='upper left')
plt.show()
plt.savefig(args.output,bbox_inches="tight")

if fit_check == 1:
    print("slope = "+str(slope))
elif fit_check == 2:
    check=np.linalg.cond(pcov)
    print("Condition # of Covariance Matrix = "+str(check))
    print("")
    print(np.diag(pcov))
