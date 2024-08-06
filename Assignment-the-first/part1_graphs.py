#!/usr/bin/env python
import argparse
import matplotlib.pyplot as plt

#get file name
def get_args():
    parser = argparse.ArgumentParser(description="Demultiplex interpreter")
    parser.add_argument("-f", "--filename",
                     help="input file",
                     required=False, type=str)
    return parser.parse_args()

args = get_args()
file_name=args.filename

first_line=True
index=[] #pos
means=[] #mean by pos

for line in open(file_name):
    #skip header
    if first_line:
        first_line=False
        continue
    #isolate relevant info
    qual_by_pos_info=line.split()
    index.append(int(qual_by_pos_info[0]))
    means.append(float(qual_by_pos_info[1]))

#subplot so I can remove spines
fig,ax=plt.subplots()
    
plt.plot(index,means, color="mediumaquamarine",alpha=0.8,linewidth=2,linestyle="--")

#detailed labeling and naming
ax.set(title="Average Phred33 Quality by Position",
       ylabel="Quality",xlabel="Position")
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False) 

plt.savefig(file_name[:file_name.index(".")]+".png")
