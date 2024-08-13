#!/usr/bin/env python3

"""
2024-06-18

Script uses the biopython module to calculate the solvent-accesible surface-area (SASA)
of a provided PDB structure file.
Biopython's implemntation of the Shrake & Rupley algorithm is used using default parameters

Note: Code is formatted to handle output from Mike Sennett's ESR scripts i.e. the naming convention of the recon msa (line 172)

CVRSharma
Theobald Lab, Brandeis University
Date created : 2024-06-18
Date modified : 2024-06-24

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
import argparse
import matplotlib.pyplot as plt
import scipy.stats as sc

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
    return mistakes

""" Input arguments and help options """
parser = argparse.ArgumentParser(description = 'Calculate SASA values for residues in a PDB file')

parser.add_argument("-s", "--seq",  dest = "seq", help = "Name of sequence to pull from input MSA 1.")
parser.add_argument("-p", "--pdb",  dest = "pdb", help = "PDB File associated with sequence.")
parser.add_argument("-c", "--chain", dest = "chain", default="A", help = "Chain in file to use. Defaults to A.")
parser.add_argument("-m1", "--msa1", dest="msa1", help = "Input MSA 1 (Extant)")
parser.add_argument("-m2", "--msa2", dest="msa2", help = "Input MSA 2 (Recon).")
parser.add_argument("-t", "--title", dest="title", help = "Alignment tool being plotted.")

# Print --help message if no argument flags are provided by user
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

""" Perform SASA calculations"""
 # Read in input arguments
seq_name=args.seq
pdb_file=args.pdb
ch_name=args.chain
msa1=aln_read(args.msa1)
msa2=aln_read(args.msa2)

# Determine number of Taxa in alignment
N = count_taxa(args.msa1)

# Grab extant sequence
extant_seq=seq_grabber(seq_name,msa1,N)

# Grab recon sequence
recon_seq=seq_grabber(seq_name+"_recon",msa2,N)

# Compare extant to recon
seq_mis=seq_compare(extant_seq,recon_seq)

# Compare extant msa to recon msa
msa_mis=msa_compare(msa1,msa2,N)
# Determine proportion of mistakes per column in alignment by nomalizing # of mistakes per column by total mistakes in alignment
norm=np.sum(msa_mis)
norm_msa_mis=msa_mis/norm

# Perform SASA calculation
raw_sasa=sasacal(pdb_file,ch_name)

# Get proportion of mistakes for columns which contain sequence residue i
raw_seq_prop_mis=[]
for i in range(len(norm_msa_mis)):
    if extant_seq[i] == "-":
        pass
    else:
        raw_seq_prop_mis.append(norm_msa_mis[i])

# Only keep values associated with sequence mistakes
sasa=[]
seq_prop_mis=[]
for k in range(len(seq_mis)):
    if seq_mis[k] == 0:
        pass
    else:
        sasa.append(raw_sasa[k])
        seq_prop_mis.append(raw_seq_prop_mis[k])

#print(sasa)
#print(seq_prop_mis)

""" Plot SASA value of mistake against proportion of mistakes of mistake containing column """
# Calculate line of best fit
slope,yinter,r2,p,std_err = sc.linregress(sasa,seq_prop_mis)

#print(slope)

#fit=slope*sasa + yinter

#Plot
plt.figure
#labels = ['Fit', seq_name]
labels = [seq_name]
#plt.plot(sasa,fit,color='orange')
plt.scatter(sasa,seq_prop_mis,color="deepskyblue")
plt.title(args.title)
plt.xlabel("SASA")
plt.ylabel("Proportion of Mistakes in Column of Residue")
plt.text(np.max(sasa)*.8,np.max(seq_prop_mis)*.9, "$R^2$ = "+str(round(r2,3))+"\np-value = "+str(round(p,3)), fontsize = 12)
plt.legend(labels,loc='upper left', bbox_to_anchor=(1,1))
plt.show()
plt.savefig("SASAvMis",bbox_inches="tight")

