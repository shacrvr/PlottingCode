#!/usr/bin/env python3

"""
2024-06-28

Hard-coded to generate a SASA vs Proption of total mistakes of a set plot for a protein family

Plots the SASA of sequence reconstruction mistakes against the proportion
of mistakes for that residue containing alignment column

Proportion of Mistakes is calculated using total number of mistakes for an alignment set for a 
specific protein family

Script uses the biopython module to calculate the solvent-accesible surface-area (SASA)
of a provided PDB structure file.
Biopython's implemntation of the Shrake & Rupley algorithm is used using default parameters

Note: Code is formatted to handle output from Mike Sennett's ESR scripts i.e. the naming convention of the recon msa (line 172)

CVRSharma
Theobald Lab, Brandeis University
Date created : 2024-06-28
Date modified : 
"""

# Import modules
import Bio
from Bio.PDB import *
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.PDBList import PDBList
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio import AlignIO
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.stats as sc
import glob

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

# Calculate per residue SASA for a PDB file
def sasacal(pdbin,chname):
    downloader=PDBList()
    pdb_file=pdbin
    pdb=pdb_file.strip(".pdb")
    #pdb_file=downloader.retrieve_pdb_file(pdb_code=pdb)
    parser=PDBParser(QUIET=1)

    structure = parser.get_structure(pdb,pdb_file)

    chainXname=chname
    mod_ch1=Model(1)

    chain1 = Chain(chainXname)

    ch1_struct=Structure("ch1")

    num_count=0
    for resi_ind, resi in enumerate(structure[0][chainXname].get_residues()):

        res_id = resi.get_full_id()
        resi_name=resi.get_resname()
        residue = Residue(res_id,resi_name,' ')
        if resi_name !='HOH':
            for at in structure[0][chainXname][res_id[3]].get_atoms():
                residue.add(at)
            chain1.add(residue)
    mod_ch1.add(chain1)
    ch1_struct.add(mod_ch1)

    sr = ShrakeRupley()
    sr.compute(ch1_struct[1], level="R")
    out=[]
    for res in ch1_struct[1][chainXname]:
        #out.append((res.get_resname(),round(res.sasa,2)))
        out.append(round(res.sasa,3))
        #out.append(res.sasa)

    out_arr=np.asarray(out)
    return out_arr

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

""" Perform SASA calculations"""
 # Read in input arguments
seq_name = sys.argv[1]
prot_fam = sys.argv[2]

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

# Get AlphaFold rank 1 pdb model
pdb = translator("../SASA_scripts_structures/"+str(seq_name)+"*/*rank_001*.pdb")

# Perform SASA calculation
sasa=sasacal(pdb,"A")

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

# Only keep values associated with sequence mistakes for plots
ba_sasa=np.multiply(ba_seq_mis,sasa)
ba_msa_mis=[]
for i in range(len(ba_ex)):
    if ba_ex[i] == "-":
        pass
    else:
        ba_msa_mis.append(norm_mis_ba[i])
ba_arr_mis=np.asarray(ba_msa_mis)
ba_col_mis = np.multiply(ba_seq_mis,ba_arr_mis)

co_sasa=np.multiply(co_seq_mis,sasa)
co_msa_mis=[]
for i in range(len(co_ex)):
    if co_ex[i] == "-":
        pass
    else:
        co_msa_mis.append(norm_mis_co[i])
co_arr_mis=np.asarray(co_msa_mis)
co_col_mis = np.multiply(co_seq_mis,co_arr_mis)

ff_sasa=np.multiply(ff_seq_mis,sasa)
ff_msa_mis=[]
for i in range(len(ff_ex)):
    if ff_ex[i] == "-":
        pass
    else:
        ff_msa_mis.append(norm_mis_ff[i])
ff_arr_mis=np.asarray(ff_msa_mis)
ff_col_mis = np.multiply(ff_seq_mis,ff_arr_mis)

fs_sasa=np.multiply(fs_seq_mis,sasa)
fs_msa_mis=[]
for i in range(len(fs_ex)):
    if fs_ex[i] == "-":
        pass
    else:
        fs_msa_mis.append(norm_mis_fs[i])
fs_arr_mis=np.asarray(fs_msa_mis)
fs_col_mis = np.multiply(fs_seq_mis,fs_arr_mis)

li_sasa=np.multiply(li_seq_mis,sasa)
li_msa_mis=[]
for i in range(len(li_ex)):
    if li_ex[i] == "-":
        pass
    else:
        li_msa_mis.append(norm_mis_li[i])
li_arr_mis=np.asarray(li_msa_mis)
li_col_mis = np.multiply(li_seq_mis,li_arr_mis)

mu_sasa=np.multiply(mu_seq_mis,sasa)
mu_msa_mis=[]
for i in range(len(mu_ex)):
    if mu_ex[i] == "-":
        pass
    else:
        mu_msa_mis.append(norm_mis_mu[i])
mu_arr_mis=np.asarray(mu_msa_mis)
mu_col_mis = np.multiply(mu_seq_mis,mu_arr_mis)

p2_sasa=np.multiply(p2_seq_mis,sasa)
p2_msa_mis=[]
for i in range(len(p2_ex)):
    if p2_ex[i] == "-":
        pass
    else:
        p2_msa_mis.append(norm_mis_p2[i])
p2_arr_mis=np.asarray(p2_msa_mis)
p2_col_mis = np.multiply(p2_seq_mis,p2_arr_mis)

pk_sasa=np.multiply(pk_seq_mis,sasa)
pk_msa_mis=[]
for i in range(len(pk_ex)):
    if pk_ex[i] == "-":
        pass
    else:
        pk_msa_mis.append(norm_mis_pk[i])
pk_arr_mis=np.asarray(pk_msa_mis)
pk_col_mis = np.multiply(pk_seq_mis,pk_arr_mis)

pb_sasa=np.multiply(pb_seq_mis,sasa)
pb_msa_mis=[]
for i in range(len(pb_ex)):
    if pb_ex[i] == "-":
        pass
    else:
        pb_msa_mis.append(norm_mis_pb[i])
pb_arr_mis=np.asarray(pb_msa_mis)
pb_col_mis = np.multiply(pb_seq_mis,pb_arr_mis)

tf_sasa=np.multiply(tf_seq_mis,sasa)
tf_msa_mis=[]
for i in range(len(tf_ex)):
    if tf_ex[i] == "-":
        pass
    else:
        tf_msa_mis.append(norm_mis_tf[i])
tf_arr_mis=np.asarray(tf_msa_mis)
tf_col_mis = np.multiply(tf_seq_mis,tf_arr_mis)

""" Generate Plot of SASA vs Propotion of Total Mistakes for Incorrect residues """
plt.figure
labels = ['BAli-Phy',  'CLUSTALO',  'FFTNS1',  'FSA',  'LINSI',  'MUSCLE',  'PAGAN2',  'PRANK',  'PROBCONS',  'TCOFFEE']
c = ['grey','lightcoral','peru','darkorange','olive','deepskyblue','orchid','firebrick','lightsteelblue','yellowgreen']
m = ['o','v','^','p','D','*','P','>','X','<']
plt.scatter(ba_sasa,ba_col_mis,color=c[0],marker=m[0])
plt.scatter(co_sasa,co_col_mis,color=c[1],marker=m[1])
plt.scatter(ff_sasa,ff_col_mis,color=c[2],marker=m[2])
plt.scatter(fs_sasa,fs_col_mis,color=c[3],marker=m[3])
plt.scatter(li_sasa,li_col_mis,color=c[4],marker=m[4])
plt.scatter(mu_sasa,mu_col_mis,color=c[5],marker=m[5])
plt.scatter(p2_sasa,p2_col_mis,color=c[6],marker=m[6])
plt.scatter(pk_sasa,pk_col_mis,color=c[7],marker=m[7])
plt.scatter(pb_sasa,pb_col_mis,color=c[8],marker=m[8])
plt.scatter(tf_sasa,tf_col_mis,color=c[9],marker=m[9])
plt.title(prot_fam)
plt.xlabel("SASA")
plt.ylabel("Proportion of Total Mistakes")
plt.legend(labels, loc='upper left', bbox_to_anchor=(1,1))
plt.show()
plt.savefig("SASAvAllMis" ,bbox_inches="tight")