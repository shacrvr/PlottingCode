#!/usr/bin/env python3

"""
Code will calculate the Levenshtein distance between two alignments
Adapted from psuedocode from : Hjelmqvist, Sten (26 March 2012), Fast, memory efficient Levenshtein algorithm

CVRS
2023-11-29
Theobald Lab, Brandeis University

"""

# Import dependencies
from tqdm import tqdm
import numpy as np
import argparse
import sys
import Bio
from Bio import AlignIO

parser = argparse.ArgumentParser(description = 'LEViathan: Levenshtein Distance Calculator for Alignments')

# Arguments and help options

parser.add_argument("-a", "--aln1", dest = "aln1", help = "Alignment 1")
parser.add_argument("-b", "--aln2", dest = "aln2", help = "Alignment 2")
parser.add_argument("-f", "--format", dest = "format", default = "fasta", help = "Optional: Alignemnt file format. Default is FASTA")

# Print --help message if no argument flags are provided by user
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

""" Define helper functions """

# Read in alignment
def aln_read(aln):
    aln_obj = AlignIO.read(open(aln), args.format)
    return aln_obj

# Determine number of taxa in alignment
def count_taxa(file):
    with open(file, "r") as handle:
        records = list(Bio.SeqIO.parse(handle, args.format))
        n_taxa = len(records)

    return n_taxa
# Calculate Ld for two strings
def lev(string1,string2):
    # Initilize vectors
    n = len(string1)+1
    m = len(string2)+1
    v0 = np.zeros(m)
    v1 = v0.copy()

    # Populate vector v1
    for i in range(1,n):
        # Initialize first element of vector v1
        v1[0] = i
        # Step through
        for j in range(1,m):
            # calculate cost
            cost = 0 if string1[i-1] == string2[j-1] else 1
            v1[j] = min(
                v0[j] + 1, #Deletion
                v1[j-1]+1, #Insertion
                v0[j-1]+cost, #Substitution 
            )
        # Swap v0 with v1 and conitnue iteration
        v0=v1.copy()
    
    return v1[m-1]

""" Calculate average Levenshtein distance between the alignments """

# Read in alignments
aln1 = aln_read(args.aln1)
aln2 = aln_read(args.aln2)

N = count_taxa(args.aln1)
# Comapre alignments sequence by sequence
tot_Ld = 0
for i in tqdm(range(N)):
    record1 = aln1[i]
    record2 = aln2[i]
    # Check that seq IDs match
    if record1.id == record2.id:
        Ld = lev(record1, record2)
    else:
        for j in range(N):
            if aln2[j].id == record1.id:
                Ld = lev(record1,aln2[j])
                break
            else:
                pass
    tot_Ld+=Ld

print("=================================================================================")
print("Average Levenshtein Distance between Alignments = "+str(tot_Ld/N))
print("Total Levenshtein Distance between Alignments = "+str(tot_Ld))
print("Number of Sequence Pairs Examined : "+str(N)+"\n")








    
    
