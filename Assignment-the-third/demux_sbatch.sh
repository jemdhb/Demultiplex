#!/bin/bash
#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --mem=100GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --job-name=demux                  #optional: job name
conda activate bgmp_demultiplex

/usr/bin/time -v python part3.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz\
 -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
 -r3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \
 -r4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz