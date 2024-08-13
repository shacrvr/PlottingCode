#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.stats as sc
import sys
import argparse

parser = argparse.ArgumentParser(description = 'Use to quickly plot a histogram')

# Arguments and help options

parser.add_argument("-n", "--nbins", dest = "nbins",type=int, help = "Number of bins")
parser.add_argument("-d", "--data", dest = "data", help = "Data")
parser.add_argument("-c", "--colour", dest = "colour", default="blue", help = "Bin Colour; default is blue")
parser.add_argument("-xl", "--xlab", dest = "xlab", help = "X-label")
parser.add_argument("-t", "--title", dest = "title", help = "Title")
parser.add_argument("-o", "--output", dest = "output", help = "Output file name")
parser.add_argument("--cum", dest = "cum",action="store_true", help = "Plot as a cumulative histogram")

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
    entry = col[1:]
    for i in entry:
        # Special try/except statement to handle when entry is empty i.e all worong/correct reconstruction (see Fish LDH dataset)
        try:
            ii = float(i.strip("[]").strip(","))
            dat.append(ii)
        except ValueError:
            print(col[0])

#Plot
plt.figure
num_bins=args.nbins
if args.cum is True:
    c,b = np.histogram(dat,num_bins)
    pdf = c/sum(c)
    cdf = np.cumsum(pdf)
    plt.plot(b[1:],cdf,color = args.colour)
    plt.ylabel("CDF")

else:
    plt.hist(dat,num_bins, color = args.colour)
    plt.ylabel('Count')

plt.title(args.title)
plt.xlabel(args.xlab)
#plt.xlim(0,1)
#plt.ylim(0,1)
plt.legend(['Mean = {0}'.format(round(np.mean(dat),3))],loc='upper left')
plt.show()
plt.savefig(args.output)

print("mean = "+str(np.mean(dat)))
