import bioinfo
import numpy as np
import argparse
import gzip   
import matplotlib.pyplot as plt
index=[]
means=[]
first_line=True
for line in open("R1_OLD_means.txt"):
    if first_line:
        first_line=False
        continue
    l=line.split()
    index.append(int(l[0]))
    means.append(float(l[1]))

fig,ax=plt.subplots()
    
plt.bar(x=index,height=means, color="mediumaquamarine",alpha=0.8,width=8)
#detailed labeling and naming
ax.set(title="Average Phred33 Quality by Position",
       ylabel="Quality",xlabel="Position")
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False) 
plt.savefig("TEST.png")
