#!/bin/bash Python
import bioinfo
import numpy as np
import argparse
import gzip   

#get filename
def get_args():
    parser = argparse.ArgumentParser(description="Demultiplex interpreter")
    parser.add_argument("-f", "--filename",
                     help="input file",
                     required=False, type=str)
    return parser.parse_args()
args = get_args()
file_name=args.filename

#determine read-type from filename
run=file_name[file_name.rfind("_")-2:file_name.rfind("_")]

if run=="R2" or run=="R3":
    array_len=8 #indexes
else:
    array_len=101 #biological

#will populate with running sums as we move through file
mean_qualities=np.zeros(array_len,)
#figured rb was slightly faster
with gzip.open(file_name,"rb") as zipf:
    line_num=-1
    while True:
        line_num+=1
        line=zipf.readline().decode().strip()
        #if eof
        if line=="":
            break     
        #if quality line       
        if line_num%4==3:
            for pos,char in enumerate(list(line)):
                #running sum
                mean_qualities[pos]+=bioinfo.convert_phred(str(char))

#dividing running sums by record num (line_num/4) to obtain means
mean_qualities=[sum/(line_num/4) for sum in mean_qualities]

#whats written to the text file
print(f"index\tmean")
for i, item in enumerate(mean_qualities):
    print(f"{i}\t{item}")
print(f"avg\t{np.mean(mean_qualities)}")

#FILE PATHS
#/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
#/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
#/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
#/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz