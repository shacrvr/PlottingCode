#!/usr/bin/env python3

"""
2024-07-01

Script will output the min-max scaled to alignment tool and protein family file values to be used to color
a PDB structure in PyMol

Note: Code is formatted to handle output from Mike Sennett's ESR scripts i.e. the naming convention of the recon msa 

CVRSharma
Theobald Lab, Brandeis University
Date created : 2024-07-01
Date modified : 
"""

# Import modules
import Bio
from Bio import AlignIO
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.stats as sc
import glob

 # Read in input arguments
seq_name = sys.argv[1]

""" Define helper functions """

# Read in alignment
format="fasta"
def aln_read(aln):
    aln_obj = AlignIO.read(open(aln),format)
    return aln_obj

# Determine number of taxa in alignment
def count_taxa(file):
    with open(file, "r") as handle:
        records = list(Bio.SeqIO.parse(handle,format))
        n_taxa = len(records)

    return n_taxa

# Pull specific seqeunce by name from alignment object
def seq_grabber(name,msa,N):
    for i in range(N):
        if msa[i].id == name:
            seq=msa[i]
        else:
            pass
    return seq

# Determine total # of mistakes between two seqeunces assuming equal length and identical gaps, outputting  binary string of sequence mistakes
# excluding gaps in sequence
def seq_compare(seq1,seq2):
    mistake_list=[]
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            if seq1[i] != "-":
                mistake_list.append(float(0))
            else:
                pass
        else:
            mistake_list.append(float(1))
    mistake_arr=np.asarray(mistake_list)
    return mistake_arr

# Determine total # of mistakes per column between two MSAs assuming equal length and identical gaps
def msa_compare(msa1,msa2,N):
    n=len(msa1[0])
    mistakes=np.zeros(n)
    for i in range(N):
        seq1 = msa1[i]
        seq2 = msa2[i]
        # Ensure names match
        if seq1.id == seq2.id:
            pass
        else:
            for m in range(N):
                if msa2[m].id == seq1.id:
                    seq2 = msa2[m]
                else:
                    pass
        for j in range(len(seq1)):
            if seq1[j] == seq2[j]:
                if seq1[j] == "-":
                    pass
                else:
                    pass
            else:
                mistakes[j]+=1
    return mistakes, np.sum(mistakes)

# Prep regex file inputs to pass into python code
def translator(input):
    query=glob.glob(input)
    processed=str(query).strip("[]").strip("''")
    return processed

#Min-max scale data between 0-100
def minmax(file):
    scaled_file=[]
    min=np.min(file)
    max=np.max(file)
    diff= max - min
    for i in file:
        zi = ((i-min)/diff)*100
        scaled_file.append(zi)
    return scaled_file,min,max

def scaler(file,min,max,diff):
    scaled=[]
    for i in file:
        zi = ((i-min)/diff)*100
        scaled.append(zi)
    return scaled 

def file_writer(outpath,file):
    stdout_fileno = sys.stdout
    sys.stdout = open(outpath, "w")
    for i in file:
        print(i)
    sys.stdout.close()
    sys.stdout = stdout_fileno

""" Perform calculations"""

#BAliPhy
aln_ex_ba = translator("BAliPhy/esr/*extant_0.a2m")
msa_ex_ba = aln_read(aln_ex_ba)
aln_re_ba = translator("BAliPhy/esr/*recon_0.a2m")
msa_re_ba = aln_read(aln_re_ba)

#CLUSTALO
aln_ex_co = translator("CLUSTALO/esr/*extant_0.a2m")
msa_ex_co = aln_read(aln_ex_co)
aln_re_co = translator("CLUSTALO/esr/*recon_0.a2m")
msa_re_co = aln_read(aln_re_co)

#FFT-NS-1
aln_ex_ff = translator("FFTNS1/esr/*extant_0.a2m")
msa_ex_ff = aln_read(aln_ex_ff)
aln_re_ff = translator("FFTNS1/esr/*recon_0.a2m")
msa_re_ff = aln_read(aln_re_ff)

#FSA
aln_ex_fs = translator("FSA/esr/*extant_0.a2m")
msa_ex_fs = aln_read(aln_ex_fs)
aln_re_fs = translator("FSA/esr/*recon_0.a2m")
msa_re_fs = aln_read(aln_re_fs)

#L-INS-i
aln_ex_li = translator("LINSI/esr/*extant_0.a2m")
msa_ex_li = aln_read(aln_ex_li)
aln_re_li = translator("LINSI/esr/*recon_0.a2m")
msa_re_li = aln_read(aln_re_li)

#MUSCLE
aln_ex_mu = translator("MUSCLE/esr/*extant_0.a2m")
msa_ex_mu = aln_read(aln_ex_mu)
aln_re_mu = translator("MUSCLE/esr/*recon_0.a2m")
msa_re_mu = aln_read(aln_re_mu)

#PAGAN2
aln_ex_p2 = translator("PAGAN2/esr/*extant_0.a2m")
msa_ex_p2 = aln_read(aln_ex_p2)
aln_re_p2 = translator("PAGAN2/esr/*recon_0.a2m")
msa_re_p2 = aln_read(aln_re_p2)

#PRANK
aln_ex_pk = translator("PRANK/esr/*extant_0.a2m")
msa_ex_pk = aln_read(aln_ex_pk)
aln_re_pk = translator("PRANK/esr/*recon_0.a2m")
msa_re_pk = aln_read(aln_re_pk)

#PROBCONS
aln_ex_pb = translator("PROBCONS/esr/*extant_0.a2m")
msa_ex_pb = aln_read(aln_ex_pb)
aln_re_pb = translator("PROBCONS/esr/*recon_0.a2m")
msa_re_pb = aln_read(aln_re_pb)

#TCOFFEE
aln_ex_tf = translator("TCOFFEE/esr/*extant_0.a2m")
msa_ex_tf = aln_read(aln_ex_tf)
aln_re_tf = translator("TCOFFEE/esr/*recon_0.a2m")
msa_re_tf = aln_read(aln_re_tf)

# Determine number of taxa
N = count_taxa(aln_ex_ba)

# Generate mistake vectors and total mistakes for each alignment tool
ba_mis,ba_sum_mis = msa_compare(msa_ex_ba,msa_re_ba,N)
co_mis,co_sum_mis = msa_compare(msa_ex_co,msa_re_co,N)
ff_mis,ff_sum_mis = msa_compare(msa_ex_ff,msa_re_ff,N)
fs_mis,fs_sum_mis = msa_compare(msa_ex_fs,msa_re_fs,N)
li_mis,li_sum_mis = msa_compare(msa_ex_li,msa_re_li,N)
mu_mis,mu_sum_mis = msa_compare(msa_ex_mu,msa_re_mu,N)
p2_mis,p2_sum_mis = msa_compare(msa_ex_p2,msa_re_p2,N)
pk_mis,pk_sum_mis = msa_compare(msa_ex_pk,msa_re_pk,N)
pb_mis,pb_sum_mis = msa_compare(msa_ex_pb,msa_re_pb,N)
tf_mis,tf_sum_mis = msa_compare(msa_ex_tf,msa_re_tf,N)

# Sum total mistakes across alignments
grand_total_mis = ba_sum_mis + co_sum_mis + ff_sum_mis + fs_sum_mis + li_sum_mis + mu_sum_mis + p2_sum_mis + pk_sum_mis + pb_sum_mis + tf_sum_mis

# Normalize mistake vectors by sum total of all mistakes
norm_mis_ba = ba_mis/grand_total_mis
norm_mis_co = co_mis/grand_total_mis
norm_mis_ff = ff_mis/grand_total_mis
norm_mis_fs = fs_mis/grand_total_mis
norm_mis_li = li_mis/grand_total_mis
norm_mis_mu = mu_mis/grand_total_mis
norm_mis_p2 = p2_mis/grand_total_mis
norm_mis_pk = pk_mis/grand_total_mis
norm_mis_pb = pb_mis/grand_total_mis
norm_mis_tf = tf_mis/grand_total_mis

# Binary encode sequence as corrrect (0) or incorrect (1)
ba_ex=seq_grabber(seq_name,msa_ex_ba,N)
ba_re=seq_grabber(seq_name+"_recon",msa_re_ba,N)
ba_seq_mis=seq_compare(ba_ex,ba_re)
co_ex=seq_grabber(seq_name,msa_ex_co,N)
co_re=seq_grabber(seq_name+"_recon",msa_re_co,N)
co_seq_mis=seq_compare(co_ex,co_re)
ff_ex=seq_grabber(seq_name,msa_ex_ff,N)
ff_re=seq_grabber(seq_name+"_recon",msa_re_ff,N)
ff_seq_mis=seq_compare(ff_ex,ff_re)
fs_ex=seq_grabber(seq_name,msa_ex_fs,N)
fs_re=seq_grabber(seq_name+"_recon",msa_re_fs,N)
fs_seq_mis=seq_compare(fs_ex,fs_re)
li_ex=seq_grabber(seq_name,msa_ex_li,N)
li_re=seq_grabber(seq_name+"_recon",msa_re_li,N)
li_seq_mis=seq_compare(li_ex,li_re)
mu_ex=seq_grabber(seq_name,msa_ex_mu,N)
mu_re=seq_grabber(seq_name+"_recon",msa_re_mu,N)
mu_seq_mis=seq_compare(mu_ex,mu_re)
p2_ex=seq_grabber(seq_name,msa_ex_p2,N)
p2_re=seq_grabber(seq_name+"_recon",msa_re_p2,N)
p2_seq_mis=seq_compare(p2_ex,p2_re)
pk_ex=seq_grabber(seq_name,msa_ex_pk,N)
pk_re=seq_grabber(seq_name+"_recon",msa_re_pk,N)
pk_seq_mis=seq_compare(pk_ex,pk_re)
pb_ex=seq_grabber(seq_name,msa_ex_pb,N)
pb_re=seq_grabber(seq_name+"_recon",msa_re_pb,N)
pb_seq_mis=seq_compare(pb_ex,pb_re)
tf_ex=seq_grabber(seq_name,msa_ex_tf,N)
tf_re=seq_grabber(seq_name+"_recon",msa_re_tf,N)
tf_seq_mis=seq_compare(tf_ex,tf_re)

# Only keep values associated with sequence mistakes for output
ba_tool_norm=[]
ba_all_norm=[]
for i in range(len(ba_ex)):
    if ba_ex[i] == "-":
        pass
    else:
        ba_tool_norm.append(ba_mis[i])
        ba_all_norm.append(norm_mis_ba[i])
ba_tool_norm=ba_tool_norm/ba_sum_mis
ba_tool_arr=np.asarray(ba_tool_norm)
ba_all_arr=np.asarray(ba_all_norm)
ba_all_col = np.multiply(ba_seq_mis,ba_all_arr)
ba_tool_col = np.multiply(ba_seq_mis,ba_tool_arr)

co_tool_norm=[]
co_all_norm=[]
for i in range(len(co_ex)):
    if co_ex[i] == "-":
        pass
    else:
        co_tool_norm.append(co_mis[i])
        co_all_norm.append(norm_mis_co[i])
co_tool_norm=co_tool_norm/co_sum_mis
co_tool_arr=np.asarray(co_tool_norm)
co_all_arr=np.asarray(co_all_norm)
co_all_col = np.multiply(co_seq_mis,co_all_arr)
co_tool_col = np.multiply(co_seq_mis,co_tool_arr)

ff_tool_norm=[]
ff_all_norm=[]
for i in range(len(ff_ex)):
    if ff_ex[i] == "-":
        pass
    else:
        ff_tool_norm.append(ff_mis[i])
        ff_all_norm.append(norm_mis_ff[i])
ff_tool_norm=ff_tool_norm/ff_sum_mis
ff_tool_arr=np.asarray(ff_tool_norm)
ff_all_arr=np.asarray(ff_all_norm)
ff_all_col = np.multiply(ff_seq_mis,ff_all_arr)
ff_tool_col = np.multiply(ff_seq_mis,ff_tool_arr)

fs_tool_norm=[]
fs_all_norm=[]
for i in range(len(fs_ex)):
    if fs_ex[i] == "-":
        pass
    else:
        fs_tool_norm.append(fs_mis[i])
        fs_all_norm.append(norm_mis_fs[i])
fs_tool_norm=fs_tool_norm/fs_sum_mis
fs_tool_arr=np.asarray(fs_tool_norm)
fs_all_arr=np.asarray(fs_all_norm)
fs_all_col = np.multiply(fs_seq_mis,fs_all_arr)
fs_tool_col = np.multiply(fs_seq_mis,fs_tool_arr)

li_tool_norm=[]
li_all_norm=[]
for i in range(len(li_ex)):
    if li_ex[i] == "-":
        pass
    else:
        li_tool_norm.append(li_mis[i])
        li_all_norm.append(norm_mis_li[i])
li_tool_norm=li_tool_norm/li_sum_mis
li_tool_arr=np.asarray(li_tool_norm)
li_all_arr=np.asarray(li_all_norm)
li_all_col = np.multiply(li_seq_mis,li_all_arr)
li_tool_col = np.multiply(li_seq_mis,li_tool_arr)

mu_tool_norm=[]
mu_all_norm=[]
for i in range(len(mu_ex)):
    if mu_ex[i] == "-":
        pass
    else:
        mu_tool_norm.append(mu_mis[i])
        mu_all_norm.append(norm_mis_mu[i])
mu_tool_norm=mu_tool_norm/mu_sum_mis
mu_tool_arr=np.asarray(mu_tool_norm)
mu_all_arr=np.asarray(mu_all_norm)
mu_all_col = np.multiply(mu_seq_mis,mu_all_arr)
mu_tool_col = np.multiply(mu_seq_mis,mu_tool_arr)

p2_tool_norm=[]
p2_all_norm=[]
for i in range(len(p2_ex)):
    if p2_ex[i] == "-":
        pass
    else:
        p2_tool_norm.append(p2_mis[i])
        p2_all_norm.append(norm_mis_p2[i])
p2_tool_norm=p2_tool_norm/p2_sum_mis
p2_tool_arr=np.asarray(p2_tool_norm)
p2_all_arr=np.asarray(p2_all_norm)
p2_all_col = np.multiply(p2_seq_mis,p2_all_arr)
p2_tool_col = np.multiply(p2_seq_mis,p2_tool_arr)

pk_tool_norm=[]
pk_all_norm=[]
for i in range(len(pk_ex)):
    if pk_ex[i] == "-":
        pass
    else:
        pk_tool_norm.append(pk_mis[i])
        pk_all_norm.append(norm_mis_pk[i])
pk_tool_norm=pk_tool_norm/pk_sum_mis
pk_tool_arr=np.asarray(pk_tool_norm)
pk_all_arr=np.asarray(pk_all_norm)
pk_all_col = np.multiply(pk_seq_mis,pk_all_arr)
pk_tool_col = np.multiply(pk_seq_mis,pk_tool_arr)

pb_tool_norm=[]
pb_all_norm=[]
for i in range(len(pb_ex)):
    if pb_ex[i] == "-":
        pass
    else:
        pb_tool_norm.append(pb_mis[i])
        pb_all_norm.append(norm_mis_pb[i])
pb_tool_norm=pb_tool_norm/pb_sum_mis
pb_tool_arr=np.asarray(pb_tool_norm)
pb_all_arr=np.asarray(pb_all_norm)
pb_all_col = np.multiply(pb_seq_mis,pb_all_arr)
pb_tool_col = np.multiply(pb_seq_mis,pb_tool_arr)

tf_tool_norm=[]
tf_all_norm=[]
for i in range(len(tf_ex)):
    if tf_ex[i] == "-":
        pass
    else:
        tf_tool_norm.append(tf_mis[i])
        tf_all_norm.append(norm_mis_tf[i])
tf_tool_norm=tf_tool_norm/tf_sum_mis
tf_tool_arr=np.asarray(tf_tool_norm)
tf_all_arr=np.asarray(tf_all_norm)
tf_all_col = np.multiply(tf_seq_mis,tf_all_arr)
tf_tool_col = np.multiply(tf_seq_mis,tf_tool_arr)

""" Determine min and max for all data sets  """ 
scaled_ba,ba_min,ba_max=minmax(ba_tool_col)
all_scaled_ba,all_ba_min,all_ba_max=minmax(ba_all_col)
tool_min=ba_min
tool_max=ba_max
all_min=all_ba_min
all_max=all_ba_max

scaled_co,co_min,co_max=minmax(co_tool_col)
all_scaled_co,all_co_min,all_co_max=minmax(co_all_col)
if co_min < tool_min:
    tool_min = co_min
else:
    pass
if co_max > tool_max:
    tool_max = co_max
else:
    pass
if all_co_min < all_min:
    all_min = all_co_min
else:
    pass
if co_max > all_max:
    all_max = all_co_max
else:
    pass

scaled_ff,ff_min,ff_max=minmax(ff_tool_col)
all_scaled_ff,all_ff_min,all_ff_max=minmax(ff_all_col)
if ff_min < tool_min:
    tool_min = ff_min
else:
    pass
if ff_max > tool_max:
    tool_max = ff_max
else:
    pass
if all_ff_min < all_min:
    all_min = all_ff_min
else:
    pass
if ff_max > all_max:
    all_max = all_ff_max
else:
    pass

scaled_fs,fs_min,fs_max=minmax(fs_tool_col)
all_scaled_fs,all_fs_min,all_fs_max=minmax(fs_all_col)
if fs_min < tool_min:
    tool_min = fs_min
else:
    pass
if fs_max > tool_max:
    tool_max = fs_max
else:
    pass
if all_fs_min < all_min:
    all_min = all_fs_min
else:
    pass
if fs_max > all_max:
    all_max = all_fs_max
else:
    pass

scaled_li,li_min,li_max=minmax(li_tool_col)
all_scaled_li,all_li_min,all_li_max=minmax(li_all_col)
if li_min < tool_min:
    tool_min = li_min
else:
    pass
if li_max > tool_max:
    tool_max = li_max
else:
    pass
if all_li_min < all_min:
    all_min = all_li_min
else:
    pass
if li_max > all_max:
    all_max = all_li_max
else:
    pass

scaled_mu,mu_min,mu_max=minmax(mu_tool_col)
all_scaled_mu,all_mu_min,all_mu_max=minmax(mu_all_col)
if mu_min < tool_min:
    tool_min = mu_min
else:
    pass
if mu_max > tool_max:
    tool_max = mu_max
else:
    pass
if all_mu_min < all_min:
    all_min = all_mu_min
else:
    pass
if mu_max > all_max:
    all_max = all_mu_max
else:
    pass

scaled_p2,p2_min,p2_max=minmax(p2_tool_col)
all_scaled_p2,all_p2_min,all_p2_max=minmax(p2_all_col)
if p2_min < tool_min:
    tool_min = p2_min
else:
    pass
if p2_max > tool_max:
    tool_max = p2_max
else:
    pass
if all_p2_min < all_min:
    all_min = all_p2_min
else:
    pass
if p2_max > all_max:
    all_max = all_p2_max
else:
    pass

scaled_pk,pk_min,pk_max=minmax(pk_tool_col)
all_scaled_pk,all_pk_min,all_pk_max=minmax(pk_all_col)
if pk_min < tool_min:
    tool_min = pk_min
else:
    pass
if pk_max > tool_max:
    tool_max = pk_max
else:
    pass
if all_pk_min < all_min:
    all_min = all_pk_min
else:
    pass
if pk_max > all_max:
    all_max = all_pk_max
else:
    pass

scaled_pb,pb_min,pb_max=minmax(pb_tool_col)
all_scaled_pb,all_pb_min,all_pb_max=minmax(pb_all_col)
if pb_min < tool_min:
    tool_min = pb_min
else:
    pass
if pb_max > tool_max:
    tool_max = pb_max
else:
    pass
if all_pb_min < all_min:
    all_min = all_pb_min
else:
    pass
if pb_max > all_max:
    all_max = all_pb_max
else:
    pass

scaled_tf,tf_min,tf_max=minmax(tf_tool_col)
all_scaled_tf,all_tf_min,all_tf_max=minmax(tf_all_col)
if tf_min < tool_min:
    tool_min = tf_min
else:
    pass
if tf_max > tool_max:
    tool_max = tf_max
else:
    pass
if all_tf_min < all_min:
    all_min = all_tf_min
else:
    pass
if tf_max > all_max:
    all_max = all_tf_max
else:
    pass

tool_diff = tool_max - tool_min
all_diff = all_max - all_min

ba_tool_scaled=scaler(ba_tool_col,tool_min,tool_max,tool_diff)
ba_all_scaled=scaler(ba_all_col,all_min,all_max,all_diff)

co_tool_scaled=scaler(co_tool_col,tool_min,tool_max,tool_diff)
co_all_scaled=scaler(co_all_col,all_min,all_max,all_diff)

ff_tool_scaled=scaler(ff_tool_col,tool_min,tool_max,tool_diff)
ff_all_scaled=scaler(ff_all_col,all_min,all_max,all_diff)

fs_tool_scaled=scaler(fs_tool_col,tool_min,tool_max,tool_diff)
fs_all_scaled=scaler(fs_all_col,all_min,all_max,all_diff)

li_tool_scaled=scaler(li_tool_col,tool_min,tool_max,tool_diff)
li_all_scaled=scaler(li_all_col,all_min,all_max,all_diff)

mu_tool_scaled=scaler(mu_tool_col,tool_min,tool_max,tool_diff)
mu_all_scaled=scaler(mu_all_col,all_min,all_max,all_diff)

p2_tool_scaled=scaler(p2_tool_col,tool_min,tool_max,tool_diff)
p2_all_scaled=scaler(p2_all_col,all_min,all_max,all_diff)

pk_tool_scaled=scaler(pk_tool_col,tool_min,tool_max,tool_diff)
pk_all_scaled=scaler(pk_all_col,all_min,all_max,all_diff)

pb_tool_scaled=scaler(pb_tool_col,tool_min,tool_max,tool_diff)
pb_all_scaled=scaler(pb_all_col,all_min,all_max,all_diff)

tf_tool_scaled=scaler(tf_tool_col,tool_min,tool_max,tool_diff)
tf_all_scaled=scaler(tf_all_col,all_min,all_max,all_diff)

""" Write outputs to file """
# BAliPhy
file_writer("BAliPhy/ba_tool_scaled-d.txt",ba_tool_scaled)
file_writer("BAliPhy/ba_all_scaled-d.txt",ba_all_scaled)
file_writer("BAliPhy/self_scaled_ba-d.txt",scaled_ba)
file_writer("BAliPhy/self_all_scaled_ba-d.txt",all_scaled_ba)

# CLUSTALO
file_writer("CLUSTALO/co_tool_scaled-d.txt",co_tool_scaled)
file_writer("CLUSTALO/co_all_scaled-d.txt",co_all_scaled)
file_writer("CLUSTALO/self_scaled_co-d.txt",scaled_co)
file_writer("CLUSTALO/self_all_scaled_co-d.txt",all_scaled_co)

# FFT-NS-1
file_writer("FFTNS1/ff_tool_scaled-d.txt",ff_tool_scaled)
file_writer("FFTNS1/ff_all_scaled-d.txt",ff_all_scaled)
file_writer("FFTNS1/self_scaled_ff-d.txt",scaled_ff)
file_writer("FFTNS1/self_all_scaled_ff-d.txt",all_scaled_ff)

# FSA
file_writer("FSA/fs_tool_scaled-d.txt",fs_tool_scaled)
file_writer("FSA/fs_all_scaled-d.txt",fs_all_scaled)
file_writer("FSA/self_scaled_fs-d.txt",scaled_fs)
file_writer("FSA/self_all_scaled_fs-d.txt",all_scaled_fs)

# L-INS-i
file_writer("LINSI/li_tool_scaled-d.txt",li_tool_scaled)
file_writer("LINSI/li_all_scaled-d.txt",li_all_scaled)
file_writer("LINSI/self_scaled_li-d.txt",scaled_li)
file_writer("LINSI/self_all_scaled_li-d.txt",all_scaled_li)

# MUSCLE
file_writer("MUSCLE/mu_tool_scaled-d.txt",mu_tool_scaled)
file_writer("MUSCLE/mu_all_scaled-d.txt",mu_all_scaled)
file_writer("MUSCLE/self_scaled_mu-d.txt",scaled_mu)
file_writer("MUSCLE/self_all_scaled_mu-d.txt",all_scaled_mu)

# PAGAN 2
file_writer("PAGAN2/p2_tool_scaled-d.txt",p2_tool_scaled)
file_writer("PAGAN2/p2_all_scaled-d.txt",p2_all_scaled)
file_writer("PAGAN2/self_scaled_p2-d.txt",scaled_p2)
file_writer("PAGAN2/self_all_scaled_p2-d.txt",all_scaled_p2)

# PRANK
file_writer("PRANK/pk_tool_scaled-d.txt",pk_tool_scaled)
file_writer("PRANK/pk_all_scaled-d.txt",pk_all_scaled)
file_writer("PRANK/self_scaled_pk-d.txt",scaled_pk)
file_writer("PRANK/self_all_scaled_pk-d.txt",all_scaled_pk)

# PROBCONS
file_writer("PROBCONS/pb_tool_scaled-d.txt",pb_tool_scaled)
file_writer("PROBCONS/pb_all_scaled-d.txt",pb_all_scaled)
file_writer("PROBCONS/self_scaled_pb-d.txt",scaled_pb)
file_writer("PROBCONS/self_all_scaled_pb-d.txt",all_scaled_pb)

# PROBCONS
file_writer("TCOFFEE/tf_tool_scaled-d.txt",tf_tool_scaled)
file_writer("TCOFFEE/tf_all_scaled-d.txt",tf_all_scaled)
file_writer("TCOFFEE/self_scaled_tf-d.txt",scaled_tf)
file_writer("TCOFFEE/self_all_scaled_tf-d.txt",all_scaled_tf)