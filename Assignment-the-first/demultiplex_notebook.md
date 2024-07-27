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
zcat 1294_S1_L008_R4_001.fastq.gz | grep "^@K00337" -A 1 | grep -vE "^@K00337|--" | awk '{print length}' | sort -n | uniq -c
zcat 1294_S1_L008_R1_001.fastq.gz | grep "^@K00337" -A 1 | grep -vE "^@K00337|--" | awk '{print length}' | sort -n | uniq -c
```

CAN ASSUME ALL SAME LENGTH, SIMPLIFY

With this I found all my biological reads have length of 101

```bash
zcat 1294_S1_L008_R2_001.fastq.gz | grep "^@K00337" -A 1 | grep -vE "^@K00337|--" | awk '{print length}' | sort -n | uniq -c
zcat 1294_S1_L008_R3_001.fastq.gz | grep "^@K00337" -A 1 | grep -vE "^@K00337|--" | awk '{print length}' | sort -n | uniq -c
```

And my indexed reads have a length of 8

### Determined phred encoding

I isolated my qualities with

```bash
zcat 1294_S1_L008_R2_001.fastq.gz | head -n 40 | grep "^@K00337" -B 1 | grep -vE "^@K00337|--" > qualities_R2.txt
```

R1 average: `45.02145214521452`
R2 average: `35.203125`
R3 average: `33.44642857142857`
R4 average: `51.603960396039604`

R1-R2
R4-R3

if barcode unknown on one pair:
keep the other
if read fails quality check on one pair:
keep the other 