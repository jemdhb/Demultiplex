#!/bin/bash Python
import bioinfo
import numpy as np
import argparse
def get_args():
    parser = argparse.ArgumentParser(description="Demultiplex interpreter")
    parser.add_argument("-f", "--filename",
                     help="input file",
                     required=False, type=str)
    return parser.parse_args()
args = get_args()

file_name=args.filename
run=file_name[file_name.rfind("_")+1:file_name.rfind(".")]

if run=="R2" or run=="R3":
    array_len=8 #indexes
else:
    array_len=101 #biological
q=np.zeros(array_len,)
for line_num,line in enumerate(open(file_name,"r")):
    line=line.strip()
    for pos,char in enumerate(line):
        #running sum
        q[pos]+=bioinfo.convert_phred(char)
q=[sum/line_num for sum in q]
print("per base read: ",q)
print(np.mean(q))
