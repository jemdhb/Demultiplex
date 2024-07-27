# Assignment the First

## Part 1

1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz |  |  | 33 |
| 1294_S1_L008_R2_001.fastq.gz |  |  | 33 |
| 1294_S1_L008_R3_001.fastq.gz |  |  | 33 |
| 1294_S1_L008_R4_001.fastq.gz |  |  | 33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. **YOUR ANSWER HERE**
    3. **YOUR ANSWER HERE**

## Part 2

1. Define the problem
    1. We need to demultiplex are four read files, throwing out any biological-index pairs.
    with unknown or hopped indices
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
                    format_fastq(f2_header, f2_index, f2_sequence, f2_quality)
                    fastq_file.write(fastq_record)

```

5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement

```bash
def demultiplex(file1, index1, file2, index2):
"""driver function to read through our four provided files and perform demultiplexing"""
    return True #if demultiplexing was successful
    return False #otherwise

def verify_index_match(index1, index2):
"""function to determine if index1 and index1 are matching
"""
return True #for match
return False #for mismatch

def search_for_index_match(index):
"""check to see if index is in our known list of indices
"""
return True #if match is found
return False #if not match found

def format_fastq(str header, str index, str sequence, str quality):
"""Format our known, high quality read into a fastq record, appending on the 
added index information
"""
return formatted_fastq # our formatted fastq record as a string

def calculate_per_base_quality(str quality):
"""For each base in the quality string, calculate the phred33 quality score
"""
return quality_lits #list of floats of quality scores at each index

def check_quality_thresholds(list qualities,index=True):
"""iterate through qualities verifying that each bp meets our threshold
(which will change whether were examining index data or biological read data)
if any bp is below the threshold, return False. If you fully iterate through the
read, return True
"""
return boolean #^see above for specifications

def determine_output_file(str header, str sequence, str index, str quality):
"""Based on the record information provided, determine unique output file name
"""
return output_file_name #string of output file name
```
