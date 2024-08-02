# Assignment the First

## Part 1

1. Be sure to upload your Python script. Provide a link to it here: 
[part1.py](https://github.com/jemdhb/Demultiplex/blob/master/Assignment-the-first/part1.py)
[part1_graphs.py](https://github.com/jemdhb/Demultiplex/blob/master/Assignment-the-first/part1_graphs.py)


| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | R1  | 101 | 33 |
| 1294_S1_L008_R2_001.fastq.gz | R2 | 8 | 33 |
| 1294_S1_L008_R3_001.fastq.gz | R3 | 8 | 33 |
| 1294_S1_L008_R4_001.fastq.gz | R4 | 101 | 33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
R1

![text](results/R1_means.png)

R2

![text](results/R2_means.png)

R3

![text](results/R3_means.png)

R4

![text](results/R4_means.png)

## Part 2

1. Define the problem
    1. We need to demultiplex are four read files, separating properly paired barcode reads with unknown or hopped indices. Unknown is defined as unknown barcodes, or low quality barcodes. Hopped is when two different, known barcodes are assigned to the same record. 
2. Describe output
    1. My output will be n>=6 FASTQ files each with unique relevant labels
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode

### See my #5 function definitions if the pseudocode makes no sense :heart: :duck: :heart:

```bash
def demultiplex(str F1, str I1, str F2, str I2):
    header_info=[]
    sequence_info=[]

    #open files
    with open F1, I1, F2, I2 as f1, i1, f2, i2:
        #going line by line (fn=biological record, in=index record)
        for line_f1, line_f2, line_i1, line_i2 in f1, i1, f2, i2:

            #new record, write header info
            if line_f1, line_f2, line_i1, line_i2 are headers:
                header_info=[]
                header_info.append(line_f1, line_f2, line_i1, line_i2)

            #new record, write sequence info
            if line_f1, line_f2, line_i1, line_i2 are a sequence line:
                sequence_info=[]

                #lines match and indices match
                if line_i1 == line_i2.rc and verify_index_match(line_i1, line_i2): 

                    #since indices match only need to check one
                    if search_for_index_match(line_i1)==True:
                        #if indices are matching and valid, can append all data
                        sequence_info.append(line_f1, line_f2, line_i1, line_i2)

                #lines dont match, may be able to save data with valid indices
                elif line_i1 == line_i2.rc:

                    #if lines dont match and indices are invalid
                    if search_for_index_match(line_i1)==False and\
                    search_for_index_match(line_i1)==False:
                        continue #no record to write

                    #need to check separately since indices dont match
                    #if index exists, record line-index pair
                    if search_for_index_match(line_i1) is True:
                        sequence_info.append(line_f1, line_i1)
                    if search_for_index_match(line_i2) is True:
                        sequence_info.append(line_f2, line_i2)

            if line_f1, line_f2, line_i1, line_i2 are a quality line:

                #data did not survive previous checks at sequence level
                if sequence_info or header_info == []:
                    continue #skip this record

                #if some surviving data from this record
                f1_quality=calculate_per_base_quality(line_f1)
                i1_quality=calculate_per_base_quality(line_i1)
                f2_quality=calculate_per_base_quality(line_f2)
                i2_quality=calculate_per_base_quality(line_i2)

                pull out information from header_info and sequence_info into named
                variables formatted f1_index, f1_quality, etc

                #if both parts of the pair are high quality
                if check_quality_thresholds(f1_quality)==True and check_quality_thresholds(i1_quality)==True:
                    #write to fastq
                    fastq_file=determine_output_file(f1_header, f1_index, f1_sequence, f1_quality)
                    fastq_record=format_fastq(f1_header, f1_index, f1_sequence, f1_quality)
                    fastq_file.write(fastq_record)

                #if both parts of the pair are high quality
                if check_quality_thresholds(f2_quality)==True and check_quality_thresholds(i2_quality)==True:
                    #write to fastq
                    fastq_file=determine_output_file(f2_header, f2_index, f2_sequence, f2_quality)
                    fastq_record=format_fastq(f2_header, f2_index, f2_sequence, f2_quality)
                    fastq_file.write(fastq_record)

```

5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement

### Functions reflect current state of my Part 3 code NOT the above pseudocode

```bash
def create_index_set():
    """
    Populates our global variable holding all possible barcodes in a set from indexes.txt
    """
    return None

create_index_set()
VALID_INDICES_SET==("AACAGCGA",...,"TGTTCCGT")
```

```bash
def verify_index_match(str index1, str index2):
    """
    function to determine if index1 and index1 are matching
    """
    return index1.strip().lower()==index2.strip().lower()

verify_index_match("AA","AA")==True
verify_index_match("AA","AT")==False
```

```bash
def search_for_index_match(str index):
    """check to see if index is in our known list of indices (VALID_INDICES_SET)
    """
    return index in VALID_INDICES_SET

search_for_index_match("NNN")==False
search_for_index_match("AACAGCGA")==True
```

```bash
def format_fastq(header="", index1="", index2="", sequence="", quality=""):
    """Format our known, high quality read into a fastq record, appending on the 
    added index information
    """
    return formatted_fastq #str version of our record

format_fastq("@ex","aa","aa","tttgggccc","IIIIII")==\
"""
@ex aa-aa
tttggg
+
IIIIII
"""
```

```bash
def calculate_per_base_quality(list quality):
    """For each base in the quality string, calculate the phred33 quality score
    """
    return [bioinfo.convert_phred(qual) for qual in quality] #list of quality floats

calculate_per_base_quality(["II"])==[ 40.0, 40.0 ]
```

```bash
def check_quality_thresholds(list qualities,index=True):
    """iterate through qualities verifying that each bp meets our quality threshold
    (which will change whether were examining index data or biological read data)
    if any bp is below the threshold, return False. If you fully iterate through the
    read, return True
    """
    return boolean #^see above for specifications

check_quality_thresholds([2,30,2,2,2,2,2,2])==False
check_quality_thresholds([30,30,40,40,42,42,42,42])==True
```

```bash
def reverse_complement(str seq):
    """reverse complement a sequence
    """
    return rc_seq #our string reverse complimented

reverse_complement("AA")=="TT"
```

```bash
def open_files_of_interest(set barcodes_set):
    """open all output file handles for easy access during input file traversal
    """
    return all_fh #dictionary where key is formatted barcode_FN and value is relevant fh

ALL_FILE_HANDLES= open_files_of_interest(VALID_INDICES_SET)== \
                                            {"AACAGCGA_F1": open("AACAGCGA_F1.fq","w"),
                                            "unknown_F1": open("unknown_F1.fq","w"),
                                            "hopped_F2":open("hopped_F2.fq","w"),...}
```

```bash
def close_files_of_interest():
    """close all fh in file handle dictionary created by open_files_of_interest
    """
    return None

close_files_of_interest(ALL_FILE_HANDLES)#all files in ALL_FILE_HANDLES will be closed
```

```bash
def demultiplex(str file1, str index1, str file2, str index2):
    """driver function to read through our four files and perform ~demultiplexing~
    """
    return None #(but files will be written and populated, see #3 for more details)
```
