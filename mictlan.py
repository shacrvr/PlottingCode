#!/usr/bin/env python3

"""
2024-06-27

Script uses the biopython module to calculate the solvent-accesible surface-area (SASA)
of a provided PDB structure file.
Biopython's implemntation of the Shrake & Rupley algorithm is used using default parameters
Script will additionally produce the total number of mistakes between an extant and its reconstruction

Note: Code is formatted to handle output from Mike Sennett's ESR scripts i.e. the naming convention of the recon msa (line 172)

CVRSharma
Theobald Lab, Brandeis University
Date created : 2024-06-27
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
import argparse

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

# Count non-gap charcters in sequence
def gapcount(seq):
    count=0
    for i in range(len(seq)):
        if seq[i] == "-":
            pass
        else:
            count+=1
    return count


""" Input arguments and help options """
parser = argparse.ArgumentParser(description = 'Calculate sum of SASA values for mistakes and total mistakes // Output : sum SASA, sum Mistakes, FC, len(seq)')

parser.add_argument("-s", "--seq",  dest = "seq", help = "Name of sequence to pull from input MSA 1.")
parser.add_argument("-p", "--pdb",  dest = "pdb", help = "PDB File associated with sequence.")
parser.add_argument("-c", "--chain", dest = "chain", default="A", help = "Chain in file to use. Defaults to A.")
parser.add_argument("-m1", "--msa1", dest="msa1", help = "Input MSA 1 (Extant)")
parser.add_argument("-m2", "--msa2", dest="msa2", help = "Input MSA 2 (Recon).")

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

# Perform SASA calculation
sasa=sasacal(pdb_file,ch_name)

# Array multiplication to get SASA of mistakes only
sasa_mis = sasa * seq_mis

# Perform sums
tot_sasa=np.sum(sasa_mis)
tot_mis=np.sum(seq_mis)
length=gapcount(extant_seq)

# Calculate Fraction Correct
correct = length - tot_mis
fc = correct/length

# Print sums to screen
print("  SASA  Mistakes  FC  Length")
print("================================")
print(round(tot_sasa,3)," "+str(tot_mis),"  "+str(round(fc,3))," "+str(length))









