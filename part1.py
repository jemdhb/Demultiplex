#!/bin/bash Python
import bioinfo
import numpy as np
import argparse
import gzip   

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
with gzip.open(file_name,"rb") as zipf:
    line_num=-1
    while True:
        line_num+=1
        line=zipf.readline().decode().strip()
        if line=="":
            break            
        if line_num%4==3:
            for pos,char in enumerate(list(line)):
                #running sum
                q[pos]+=bioinfo.convert_phred(str(char))
q=[sum/(line_num/4) for sum in q]
#print("read len of ",file_name,len(q))
print(f"index\tmean")
for i, item in enumerate(q):
    print(f"{i}\t{item}")
print(f"avg\t{np.mean(q)}")

#/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
#/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
#/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
#/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz