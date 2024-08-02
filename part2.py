#!/bin/bash Python
import bioinfo
import zipfile as gzip
import itertools
VALID_INDICES="/projects/bgmp/shared/2017_sequencing/indexes.txt"
VALID_INDICES_SET=set()

def create_index_set():
    first_line=True
    for line in open(VALID_INDICES,"r"):
        if first_line:
            first_line=False
            continue
        index=line.split()[-1].strip()
        VALID_INDICES_SET.add(index)
    return
def verify_index_match(index1, index2):
    """
    function to determine if index1 and index1 are matching
    """
    return index1.strip().lower()==index2.strip().lower()

def search_for_index_match(index):
    """check to see if index is in our known list of indices
    """
    return (index in VALID_INDICES_SET)

def format_fastq(header="", index1="", index2="", sequence="", quality=""):
    """Format our known, high quality read into a fastq record, appending on the 
    added index information
    """
    if index1=="": #if an unknown barcode
            formatted_fastq=header+" unknown\n"+sequence+"\n+\n"+quality

    formatted_fastq=header+" "+index1+" "+index2+"\n"+sequence+"\n+\n"+quality
    return formatted_fastq # our formatted fastq record as a string

def calculate_per_base_quality(quality):
    """For each base in the quality string, calculate the phred33 quality score
    """
    return [bioinfo.convert_phred(qual) for qual in quality]  

def check_quality_thresholds(qualities,index=True):
    """iterate through qualities verifying that each bp meets our threshold
    (which will change whether were examining index data or biological read data)
    if any bp is below the threshold, return False. If you fully iterate through the
    read, return True
    """
    thresholds=[20]*len(qualities) #trash for now
    for index,qual in enumerate(qualities):
        if qual<thresholds[index]:
            return False #found trash
        return True #found good sequence line

def reverse_complement(seq):
    """reverse complement a sequence
    """
    nucs={"A":"T","C":"G","T":"A","G":"C","N":"N"}
    rc_seq=""
    for char in seq.strip().upper():#since our dictionary is uppercase
        rc_seq=nucs[char]+rc_seq #build backwards 
    return rc_seq

def open_files_of_interest(barcodes_set):
    all_fh={}
    for barcode in barcodes_set:        
        file1_name=f"outfiles/{barcode}_{barcode}_R1.fastq"
        file2_name=f"outfiles/{barcode}_{barcode}_R2.fastq"
        all_fh[barcode+"_R1"]=open(file1_name,"w")
        all_fh[barcode+"_R2"]=open(file2_name,"w")
    #we already know the files names of these
    all_fh["unknown_R1"]=open("outfiles/unknown_R1.fastq","w")
    all_fh["unknown_R2"]=open("outfiles/unknown_R2.fastq","w")
    all_fh["hopped_R1"]=open("outfiles/hopped_R1.fastq","w")
    all_fh["hopped_R2"]=open("outfiles/hopped_R2.fastq","w")

    return all_fh


def close_files_of_interest():
    for fh in ALL_FILE_HANDLES.values():
        fh.close()
def extended_stats():
    #all hops
    l=list(itertools.permutations(list(VALID_INDICES_SET),2))
    #add matches
    for item in VALID_INDICES_SET:
        l.append((item,item))
    #add idk
    l.append("unk")
    return {key: 0 for key in l}
### See my #5 function definitions if the pseudocode makes no sense :heart: :duck: :heart:
def demultiplex(R1, I1, R2,  I2):
    header_info={"r1":"","i1":"","r2":"","i2":""}
    sequence_info={"r1":"","i1":"","r2":"","i2":""}
    stats_info=extended_stats()
    #open files
    with open(R1,"r") as r1, open(I1,"r") as i1, open(R2,"r") as r2, open(I2,"r") as i2:
        #going line by line (fn=biological record, in=index record)
        index=0
        for line_r1, line_r2, line_i1, line_i2 in zip(r1,r2, i1, i2):
            #new record, write header info
            if index%4==0:
                
                header_info={"r1":line_r1.strip(), "r2":line_r2.strip(), 
                             "i1":line_i1.strip(), "i2":line_i2.strip()}

            #new record, write sequence info
            if index%4==1:
                line_i2=reverse_complement(line_i2) #only this one?
                sequence_info={"r1":line_r1.strip(), "r2":line_r2.strip(),
                               "i1":line_i1.strip(), "i2":line_i2.strip()}

            if index%4==3:
                i1_quality=calculate_per_base_quality(line_i1.strip())
                i2_quality=calculate_per_base_quality(line_i2.strip())
                i2_quality.reverse() #because we rc the sequence line

                fastq_record1=format_fastq(header=header_info["r1"],
                                              index1=sequence_info["i1"],
                                              index2=sequence_info["i2"],
                                              sequence=sequence_info["r1"],
                                              quality=line_r1)
                
                fastq_record2=format_fastq(header=header_info["r2"],
                                              index1=sequence_info["i1"],
                                              index2=sequence_info["i2"],
                                              sequence=sequence_info["r2"],
                                              quality=line_r2)
                #NOW we have enough info to check everything
                #indices match
                if verify_index_match(sequence_info["i1"], sequence_info["i2"]): 
                    #if indices match, exist and are HQ
                    if search_for_index_match(sequence_info["i1"])==True and\
                        check_quality_thresholds(i2_quality) and\
                        check_quality_thresholds(i1_quality):
                        stats_info[(sequence_info["i1"],sequence_info["i1"])]+=1

                        ALL_FILE_HANDLES[f"{sequence_info["i1"]}_R1"].write(fastq_record1)
                        ALL_FILE_HANDLES[f"{sequence_info["i2"]}_R2"].write(fastq_record2)

                    else: #indices both match an unknown barcode, or are low quality
                        #so cant be confident in the match
                        stats_info["unk"]+=1

                        ALL_FILE_HANDLES["unknown_R1"].write(fastq_record1)
                        ALL_FILE_HANDLES["unknown_R2"].write(fastq_record2)

                #indices dont match, may be able to save data with valid indices
                elif sequence_info["i1"] != sequence_info["i2"]:                      
                    #for a barcode to hop both need to exist and be hq
                    if search_for_index_match(sequence_info["i1"]) and\
                        search_for_index_match(sequence_info["i2"]) and\
                        check_quality_thresholds(i1_quality) and\
                        check_quality_thresholds(i2_quality):
                        stats_info[(sequence_info["i1"],sequence_info["i2"])]
                        ALL_FILE_HANDLES["hopped_R1"].write(fastq_record1)
                        ALL_FILE_HANDLES["hopped_R2"].write(fastq_record2) 
                    
                    else:#otherwise UNKNOWN
                        stats_info["unk"]+=1

                        ALL_FILE_HANDLES["unknown_R1"].write(fastq_record1)
                        ALL_FILE_HANDLES["unknown_R2"].write(fastq_record2) 
                #add a catch all else
                else:
                    print("TRASH FOUND CHECK ERROR_LOG")
                    error_log.write(fastq_record1)
                    error_log.write(fastq_record2)

            index+=1
    return stats_info
create_index_set()


ALL_FILE_HANDLES=open_files_of_interest(VALID_INDICES_SET)
error_log=open("error_log.txt","w")
stats_info=demultiplex(R1="r1_test.fq",I2="r3_test.fq",I1="r2_test.fq",R2="r4_test.fq")
close_files_of_interest()

for barcodes,counts in stats_info.items():
    if barcodes == "unk":
        error_log.write(f"{barcodes}\t{counts}\n")

    else:
        error_log.write(f"{tuple(barcodes)[0]} {tuple(barcodes)[1]}\t{counts}\n")
