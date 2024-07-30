# Demultiplexing P1

The following was done in an interactive node with the following memory specifications
```srun -A bgmp -p bgmp -N 1 -c 4 --mem=100G -t 8:00:00 --pty bash```

Determined which files are indexes and which are biological reads with:

```bash
zcat 1294_S1_L008_R1_001.fastq.gz | head -n 4 #biological
zcat 1294_S1_L008_R2_001.fastq.gz | head -n 4 #index
zcat 1294_S1_L008_R3_001.fastq.gz | head -n 4 #index
zcat 1294_S1_L008_R4_001.fastq.gz | head -n 4 #biological
```

## Make into a table

### Determined the read lengths in the file

```bash
zcat 1294_S1_L008_R4_001.fastq.gz | head -n 4 | grep "^@K00337" -A 1 | \
grep -vE "^@K00337|--" | awk '{print length}' | sort -n | uniq -c

zcat 1294_S1_L008_R1_001.fastq.gz | head -n 4 | grep "^@K00337" -A 1 | \
grep -vE "^@K00337|--" | awk '{print length}' | sort -n | uniq -c
```

With this I found all my biological reads have length of 101

```bash
zcat 1294_S1_L008_R2_001.fastq.gz | head -n 4 | grep "^@K00337" -A 1 | \
grep -vE "^@K00337|--" | awk '{print length}' | sort -n | uniq -c

zcat 1294_S1_L008_R3_001.fastq.gz | head -n 4 | grep "^@K00337" -A 1 | \
grep -vE "^@K00337|--" | awk '{print length}' | sort -n | uniq -c
```

And my indexed reads have a length of 8

### Determined phred encoding

I isolated my qualities with

```bash
zcat 1294_S1_L008_R1_001.fastq.gz | head -n 40 | grep "^@K00337" -B 1 | grep -vE "^@K00337|--" > qualities_R1.txt
zcat 1294_S1_L008_R2_001.fastq.gz | head -n 40 | grep "^@K00337" -B 1 | grep -vE "^@K00337|--" > qualities_R2.txt
zcat 1294_S1_L008_R3_001.fastq.gz | head -n 40 | grep "^@K00337" -B 1 | grep -vE "^@K00337|--" > qualities_R3.txt
zcat 1294_S1_L008_R4_001.fastq.gz | head -n 40 | grep "^@K00337" -B 1 | grep -vE "^@K00337|--" > qualities_R4.txt
```

R1 average: `45.02145214521452`
R2 average: `35.203125`
R3 average: `33.44642857142857`
R4 average: `51.603960396039604`

These averages correspond with Phred-33 encoding.

I now felt confident moving forward and calculating my per nucleotide distributions

My strategy was to build a one dimensional array the size of the record, iterate through the files creating a running sum of the qualities found at that position. Once you've exited the loop, you can divide each running sum by the number of records to recieve the mean quality at that position shown below:

```bash
q=[sum/line_num for sum in q]
print("read len of ",file_name,len(q))
print(np.mean(q))
```

I wrote these results to a file before creaitng my figures to ensure I could do figure creation multiple times 

I ran my python code in the sbatch script `p1_sbatch.sh` with the following memory specifications and arguments:

```bash
#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=4                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --job-name=stats                  #optional: job name
python part1.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz > R1_means.txt
python part1.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz > R2_means.txt
python part1.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz > R3_means.txt
python part1.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz > R4_means.txt
```

I failed to time these but R1 took over 2:20 minutes.

# Part 2:

I brainstormed my approach for this in `Assignment-the-first/Answers.md`:

After recieving some feedback and reading other's code I changed the following:

1. When I read Claire's pseudocode, I realized I handled quality checking way too early which complicated things.
    1. I changed my logic to creating my records BEFORE verifying any barocde quality
2. For writing to my files I initally did the following

    ```bash
    try:
        open(file,"a")
    except:
        open(file, "w")
    ```

    1. After listening to Leslies lecture, I realized though not as slow as constantly writing, this would take an estimated 42 hours!
    2. I created the functions 
3. Created my reverse compliment function `reverse_compliment`
4. Deleted my  `determine_output_file` as I realized my if statements served the same purpose.
5. Updated my `format_fastq` to add information I was previously missing in my headers.
6. Had planned sliding qualities but while in the testing phase I kept my quality threshold to a global 30.
7. Changed my naming so all hopped and unknown barcodes are grouped in two files.

To create my second test files 100 sequences:
```bash
(base) (base) [jujo@n0349 Demultiplex]$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | head -n 400000 | tail -n 400 > r1_mid_test.fq
(base) (base) [jujo@n0349 Demultiplex]$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | head -n 40
0000 | tail -n 400 > r2_mid_test.fq
(base) (base) [jujo@n0349 Demultiplex]$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | head -n 40
0000 | tail -n 400 > r3_mid_test.fq
(base) (base) [jujo@n0349 Demultiplex]$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz | head -n 40
0000 | tail -n 400 > r4_mid_test.fq
```

To create my first test file 2 sequences:
```bash
(base) (base) [jujo@n0349 Demultiplex]$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | head -n 8 > r1_mid_test.fq
(base) (base) [jujo@n0349 Demultiplex]$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | head -n 8 > r2_mid_test.fq
(base) (base) [jujo@n0349 Demultiplex]$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | head -n 8  > r3_mid_test.fq
(base) (base) [jujo@n0349 Demultiplex]$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz | head -n 8 > r4_mid_test.fq
```