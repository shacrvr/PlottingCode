#!/usr/bin/env python3

# Script will sample from ancestor distribution for a given node provided an IQTREE2 state file and a binary state file for gapping

import sys
import argparse
import numpy as np
import math
from numpy.random import Generator, PCG64, default_rng
from datetime import datetime

parser = argparse.ArgumentParser(description = 'Ancestor Distribution Sampler')

# Arguments and help options

parser.add_argument("-f", "--state_file", dest = "state_file", help = "IQTREE2 formatted state file")
parser.add_argument("-n", "--node", dest = "node", help = "Node to sample from")
parser.add_argument("-b", "--binary_state", dest = "binary_file", help = "Binary state file for gapping")
parser.add_argument("-s", "--samples", dest = "samples", help = "Number of sequences to sample from distribution. Uses SMP gap sequence if -g flag not provided")
parser.add_argument("-g", "--gap_sampler", action='store_true', help = "Optional: provide flag, -g, to sample in/dels in addition to probabalities")
parser.add_argument("-d", "--debug", action='store_true', help = "Optional: provide flag, -d, to add list of seeds to output")

# Print --help message if no argument flags are provided by user
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

nod = str(args.node) # Node to reconstruct and sample

# Import state file
statefile = ()
statefile = open(args.state_file, 'r')
state_lines = statefile.readlines()
statefile.close()

# Import binary state file
bin_state = ()
bin_state = open(args.binary_file, 'r')
bin_lines = bin_state.readlines()
bin_state.close()

# Save posterior probability distribution for node of interest
node_state = []
for line in state_lines:
    if line.startswith("Node"+nod+'\t'):
        line = line.replace('\n','')
        node_state.append(line.split('\t'))
    else:
        pass
state_vector = []
state_vector = np.reshape(node_state,(-1,len(node_state[0])))
pp_vector = np.delete(state_vector, [0,1,2], 1)

# Save gap probability distribution for node of interest
gap_state = []
for line in bin_lines:
    if line.startswith("Node"+nod+'\t'):
        line = line.replace('\n', '')
        gap_state.append(line.split('\t'))
    else:
        pass
bin_vector = []
bin_vector = np.reshape(gap_state,(-1,len(gap_state[0])))
gap_vector = np.delete(bin_vector, [0,1,2], 1)

# Dictionary for amino acids ordered as they appear in the IQTREE2 state file (alphabetical by name)
AA_dict={"A":0, "R":1, "N":2, "D":3, "C":4, "Q":5, "E":6, "G":7, "H":8, "I":9, "L":10, "K":11, "M":12, "F":13, "P":14, "S":15, "T":16, "W":17, "Y":18, "V":19}

# Generate SMP sequence with gaps
smp_list = ['>node'+nod+'_smp']
smp_ave_prob = 0
smp_log_prob = 0
smp_gap_prob = 0
smp_total_prob = 0
res_count = 0
for i in range(len(gap_vector)):
    if gap_vector[i][0] > gap_vector[i][1]:
        smp_list.append('-')
        smp_gap_prob+=np.log(np.float(gap_vector[i][0]))
    else:
        for j in range(len(pp_vector[i])):
            position = j
            hi_prob = max(pp_vector[i])
            res_count+=1
            smp_ave_prob+=np.float(hi_prob)
            if pp_vector[i][position] == hi_prob:
                for key, val in AA_dict.items():
                    if val == position:
                        smp_list.append(key)
                        smp_log_prob+=np.log(np.float(pp_vector[i][position]))
                        smp_total_prob+=np.log(np.float(gap_vector[i][1]))
newline = 69
prob = smp_ave_prob/np.float(res_count)
smp_sum = smp_log_prob + smp_gap_prob
smp_total_prob+=smp_sum
for n in range(len(smp_list)):
    if n == 0:
        print(smp_list[n])
    elif n == newline:
         print(smp_list[n]+'\r')
         newline+=69
    else:
        print(smp_list[n],end="")
print("")
print('smp ave prob:', "%11.8f " % prob)
print('smp log prob:', "%11.8f " % smp_log_prob)
print('smp total prob:', "%11.8f " % smp_total_prob)
print("")
print("")

# Output smp to smp output file 
smp_name = 'node'+nod+'_smp.fasta'
stdout_fileno = sys.stdout
sys.stdout = open(smp_name, 'w')
newline = 69
for n in range(len(smp_list)):
    if n == 0:
        print(smp_list[n])
    elif n == newline:
         print(smp_list[n]+'\r')
         newline+=69
    else:
        print(smp_list[n],end="")
print("")
print('smp ave prob:', "%11.8f " % prob)
print('smp log prob:', "%11.8f " % smp_log_prob)
sys.stdout.close()
sys.stdout = stdout_fileno
# If sampling flag provided sample from ancestral distribution
NoneType = type(None)
if type(args.samples) != NoneType and args.gap_sampler is False:
    print('Sample Seqeunces from ancestral distribution')
    print("")
    print("")
    samp = int(args.samples) # Number of sequences to sample
    count = 0
    seed_list=['seeds used in sampling']
    while count < samp:
        count+=1
        seq = ['>node'+nod+'_'+str(count)]
        ave_prob = 0
        log_prob = 0
        res_count = 0
        # Check if site should be a gap, and, if it is, skip that position
        for i in range(len(gap_vector)):
            if gap_vector[i][0] > gap_vector[i][1]:
                seq.append('-')
            else:
                # Generate cumulative probability distribution from site probability distribution
                cum_dis = []
                for j in range(len(pp_vector[i])):
                    if j == 0:
                        cum_dis.append(np.float(pp_vector[i][j]))
                    else:
                        sum_prob = np.float(pp_vector[i][j]) + cum_dis[j-1]
                        cum_dis.append(sum_prob)
                #Generate random seed for sampling 
                dt = datetime.now()
                seed = int(round(dt.timestamp())) + np.random.default_rng().integers(low=1,high=10000) # This generates an integer seed using the system date-time in seconds
                seed_list.append(seed)
                rng = Generator(PCG64(seed)) # This calls a new instance of the generator using the seed
                roll = rng.random() # This will generate a random float using the generator and the above seed
                check = 0
                for c in range(len(cum_dis)):
                    if check == 1:
                        break
                    elif c == 0 and roll <= cum_dis[c]:
                        for key,val in AA_dict.items():
                            if val == c:
                                ave_prob+=np.float(pp_vector[i][c])
                                log_prob+=np.log(np.float(pp_vector[i][c]))
                                res_count+=1
                                seq.append(key)
                                check+=1
                                break
                    elif roll > cum_dis[c-1] and roll <= cum_dis[c] and c != 0:
                        for key,val in AA_dict.items():
                            if val == c:
                                ave_prob+=np.float(pp_vector[i][c])
                                log_prob+=np.log(np.float(pp_vector[i][c]))
                                res_count+=1
                                seq.append(key)
                                check+=1
                                break
                    else:
                        pass
        newline = 69
        prob = ave_prob/np.float(res_count)
        for n in range(len(seq)):
            if n == 0:
                print(seq[n])
            elif n == newline:
                print(seq[n]+'\r')
                newline+=69
            else:
                print(seq[n],end="")
        print("")
        print('ave prob:', "%11.8f " % prob)
        print('log prob:', "%11.8f " % log_prob)
        print("")
        print("")
    if args.debug is True:
        stdout_fileno = sys.stdout
        sys.stdout = open('sampling_seeds.txt', 'w')
        print('Seeds Used for Sampling')
        for seed in seed_list:
            print(seed)
        sys.stdout.close()
        sys.stdout = stdout_fileno
else:
    pass
# If -g flag provided sample in/dels
if args.gap_sampler is True and type(args.samples) != NoneType:
    print('In/Del Samples')
    print("")
    print("")
    samp = int(args.samples) # Number of sequences to sample
    countgaps = 0
    seed_listgap = ["Seeds used for In/Del Sampling"]
    while countgaps < samp:
        countgaps+=1
        gaps = ['>node'+nod+'_gap_s'+str(countgaps)]
        ave_prob = 0 # arithmetic average probablity of sequnce i.e. posterior probability of sequence
        log_prob = 0 # log probability of sequence
        gap_prob = 0 # log probability of gaps
        total_prob = 0 # log sum of the probablity of all residues given the probablity of there not being a gap and the probability of all gaps
        res_count = 0
        for k in range(len(gap_vector)):
            #Generate cumulative distribution for sampling in/dels
            cum_gap = []
            for j in range(len(gap_vector[k])):
                if j == 0:
                    cum_gap.append(np.float(gap_vector[k][j]))
                else:
                    sum_probgap = np.float(gap_vector[k][j]) + cum_gap[j-1]
                    cum_gap.append(sum_probgap)
            #Generate random seed for sampling
            dtgap = datetime.now()
            seedgap = int(round(dtgap.timestamp())) + np.random.default_rng().integers(low=1,high=10000) # This generates an integer seed using the system date-time in seconds
            seed_listgap.append(seedgap)
            rng = Generator(PCG64(seedgap)) # This calls a new instance of the generator using the seed
            rollgap = rng.random() # This will generate a random float using the generator and the above seed
            checkgap = 0
            for cc in range(len(cum_gap)):
                if checkgap == 1:
                    break
                elif cc == 0 and rollgap <= cum_gap[cc]:
                    gaps.append('-')
                    gap_prob+=np.log(np.float(gap_vector[k][cc]))
                    checkgap+=1
                    break
                elif rollgap > cum_gap[cc-1] and rollgap <= cum_gap[cc] and cc == 1:
                    for i in range(len(pp_vector[k])):
                        position = i
                        hi_prob = max(pp_vector[k])
                        res_count+=1
                        ave_prob+=np.float(hi_prob)
                        if pp_vector[k][position] == hi_prob:
                            for key, val in AA_dict.items():
                                if val == position:
                                    gaps.append(key)
                                    log_prob+=np.log(np.float(pp_vector[k][position]))
                                    total_prob+=np.log(np.float(gap_vector[k][cc]))
                                    checkgap+=1
                                    break
                            break
                else:
                    pass
        newline = 69
        prob = ave_prob/np.float(res_count)
        prob_sum = log_prob + gap_prob
        total_prob+= prob_sum
        for n in range(len(gaps)):
            if n == 0:
                print(gaps[n])
            elif n == newline:
                print(gaps[n]+'\r')
                newline+=70
            else:
                print(gaps[n],end="")
        print("")
        print('ave prob:', "%11.8f " % prob)
        print('log prob:', "%11.8f " % log_prob)
        print('gap prob:', "%11.8f " % gap_prob)
        print('total prob:', "%11.8f " % total_prob)
        print("")
        print("")
    if args.debug is True:
        stdout_fileno = sys.stdout
        sys.stdout = open('indel_seeds.txt', 'w')
        print('Seeds Used for In/Del Sampling')
        for seed_gap in seed_listgap:
            print(seed_gap)
        sys.stdout.close()
        sys.stdout = stdout_fileno
else:
    print("")
    print('In/Dels were not sampled.')
    print("")

