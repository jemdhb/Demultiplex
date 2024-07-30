#!/bin/bash Python
import bioinfo
import zipfile as gzip
VALID_INDICES="/projects/bgmp/shared/2017_sequencing/indexes.txt"

def verify_index_match(index1, index2):
    """
    function to determine if index1 and index1 are matching
    """
    return index1.strip().lower()==index2.strip().lower()

def search_for_index_match(index):
    """check to see if index is in our known list of indices
    """
    for line in open(VALID_INDICES,"r"):
        if line.split()[-1].strip()==index:
            return True
    return False #if never found

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
    thresholds=[26]*len(qualities) #trash for now
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

def open_files_of_interest():
    with open(VALID_INDICES,"r") as fh:
        all_barcodes=[]
        all_fh={}
        first_line=True
        for line in fh:
            if first_line:
                first_line=False
                continue
            barcode=line.split()[-1].strip()
            all_barcodes.append(barcode)
            file1_name=f"outfiles/{barcode}_{barcode}_F1.fastq"
            file2_name=f"outfiles/{barcode}_{barcode}_F2.fastq"
            all_fh[file1_name]=open(file1_name,"w")
            all_fh[file2_name]=open(file2_name,"w")
    #we already know the files names of these
    all_fh["outfiles/unknown_F1.fastq"]=open("outfiles/unknown_F1.fastq","w")
    all_fh["outfiles/unknown_F2.fastq"]=open("outfiles/unknown_F2.fastq","w")
    all_fh["outfiles/hopped_F1.fastq"]=open("outfiles/hopped_F1.fastq","w")
    all_fh["outfiles/hopped_F2.fastq"]=open("outfiles/hopped_F2.fastq","w")

    return all_barcodes, all_fh

ALL_BARCODE_FILES, ALL_FILE_HANDLES=open_files_of_interest()

def close_files_of_interest():
    for fh in ALL_FILE_HANDLES.values():
        fh.close()

### See my #5 function definitions if the pseudocode makes no sense :heart: :duck: :heart:
def demultiplex(F1, I1, F2,  I2):
    header_info={"f1":"","i1":"","f2":"","i2":""}
    sequence_info={"f1":"","i1":"","f2":"","i2":""}

    #open files
    with open(F1,"r") as f1, open(I1,"r") as i1, open(F2,"r") as f2, open(I2,"r") as i2:
        #going line by line (fn=biological record, in=index record)
        index=0
        for line_f1, line_f2, line_i1, line_i2 in zip(f1,f2, i1, i2):
            #new record, write header info
            if index%4==0:
                
                header_info={"f1":line_f1.strip(), "f2":line_f2.strip(), 
                             "i1":line_i1.strip(), "i2":line_i2.strip()}

            #new record, write sequence info
            if index%4==1:
                line_i2=reverse_complement(line_i2) #only this one?
                sequence_info={"f1":line_f1.strip(), "f2":line_f2.strip(),
                               "i1":line_i1.strip(), "i2":line_i2.strip()}

            if index%4==3:
                f1_quality=calculate_per_base_quality(line_f1)
                i1_quality=calculate_per_base_quality(line_i1)
                i2_quality=calculate_per_base_quality(line_i2)
                i2_quality.reverse() #because we rc the sequence line

                fastq_record1=format_fastq(header=header_info["f1"],
                                              index1=sequence_info["i1"],
                                              index2=sequence_info["i2"],
                                              sequence=sequence_info["f1"],
                                              quality=line_f1)
                fastq_record2=format_fastq(header=header_info["f2"],
                                              index1=sequence_info["i1"],
                                              index2=sequence_info["i2"],
                                              sequence=sequence_info["f2"],
                                              quality=line_f2)

                #NOW we have enough info to check everything
                #indices match
                if verify_index_match(sequence_info["i1"], sequence_info["i2"]): 

                    #since indices match only need to check one
                    #if indices match, exist and are HQ
                    if search_for_index_match(sequence_info["i1"])==True and\
                        check_quality_thresholds(f1_quality) and\
                        check_quality_thresholds(i1_quality):
                        outfile1=f"outfiles/{sequence_info["i1"]}_{sequence_info["i1"]}_F1.fastq"
                        outfile2=f"outfiles/{sequence_info["i2"]}_{sequence_info["i2"]}_F2.fastq"
                        ALL_FILE_HANDLES[outfile1].write(fastq_record1)
                        ALL_FILE_HANDLES[outfile2].write(fastq_record2)
                    else: #indices both match an unknown barcode, or are low quality
                        #so cant be confident in the match
                        outfile1=f"outfiles/unknown_F1.fastq"
                        outfile2=f"outfiles/unknown_F2.fastq"
                        ALL_FILE_HANDLES[outfile1].write(fastq_record1)
                        ALL_FILE_HANDLES[outfile2].write(fastq_record2)

                #indices dont match, may be able to save data with valid indices
                elif sequence_info["i1"] != sequence_info["i2"]:                      
                    #for a barcode to hop both need to exist and be hq
                    if search_for_index_match(sequence_info["i1"]) and\
                        search_for_index_match(sequence_info["i2"]) and\
                        check_quality_thresholds(i1_quality) and\
                        check_quality_thresholds(i2_quality):

                        outfile1=f"outfiles/hopped_F1.fastq"
                        outfile2=f"outfiles/hopped_F2.fastq"
                        ALL_FILE_HANDLES[outfile1].write(fastq_record1)
                        ALL_FILE_HANDLES[outfile2].write(fastq_record2) 
                    
                    else:#otherwise UNKNOWN
                        outfile1=f"outfiles/unknown_F1.fastq"
                        outfile2=f"outfiles/unknown_F2.fastq"
                        ALL_FILE_HANDLES[outfile1].write(fastq_record1)
                        ALL_FILE_HANDLES[outfile2].write(fastq_record2) 
            index+=1
demultiplex(F1="r1_mid_test.fq",I2="r3_mid_test.fq",I1="r2_mid_test.fq",F2="r4_mid_test.fq")
close_files_of_interest()
