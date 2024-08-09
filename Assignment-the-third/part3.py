#!/usr/bin/env python
import bioinfo
import gzip
import itertools
import argparse
import matplotlib.pyplot as plt
VALID_INDICES="/projects/bgmp/shared/2017_sequencing/indexes.txt"
VALID_INDICES_SET=set()

#get filename
def get_args():
    parser = argparse.ArgumentParser(description="Demultiplex interpreter")
    parser.add_argument("-r1", help="first record file",
                     required=True, type=str)
    parser.add_argument("-r2", help="first index file",
                     required=True, type=str)
    parser.add_argument("-r3", help="second index file",
                     required=True, type=str)
    parser.add_argument("-r4", help="second record file",
                     required=True, type=str)
    return parser.parse_args()
args = get_args()


def create_index_set():
    """populate the global set VALID_INDICES_SET with all valid barcodes from 
    VALID_INDICES (index.txt)
    """
    first_line=True
    for line in open(VALID_INDICES,"r"):
        #skip header
        if first_line:
            first_line=False
            continue
        #isolate barcode and add to set
        index=line.split()[-1].strip()
        VALID_INDICES_SET.add(index)
    return

def verify_index_match(index1:str, index2:str):
    """function to determine if index1 and index1 match

    Returns:
        boolean 
    """
    return index1.strip().lower()==index2.strip().lower()

def search_for_index_match(index: str):
    """check to see if index is in our known list of indices

    Returns:
        boolean
    """
    return (index in VALID_INDICES_SET)

def format_fastq(header: str ="", index1="", index2="", sequence="", quality="") -> str:
    """Format our known, high quality read into a fastq record, appending on the 
    added index information

    Returns:
        string containing the formatted fastq
    """
    #otherwise append barcodes
    formatted_fastq=header+" "+index1+" "+index2+"\n"+sequence+"\n+\n"+quality
    return formatted_fastq # our formatted fastq record as a string

def calculate_per_base_quality(quality: str):
    """For each base in the quality string, calculate the phred33 quality score

    Returns:
        list of float-based quality scores
    """
    return [bioinfo.convert_phred(qual) for qual in quality]  

def check_quality_thresholds(qualities:list,index=True):
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

def reverse_complement(seq:str):
    """reverse complement a sequence
    Returns:
        rc of input string
    """
    nucs={"A":"T","C":"G","T":"A","G":"C","N":"N"}
    rc_seq=""
    for char in seq.strip().upper():#since our dictionary is uppercase
        rc_seq=nucs[char]+rc_seq #build backwards 
    return rc_seq

def open_files_of_interest(barcodes_set:set):
    """open all possible output files of interest

    Args:
        barcodes_set (set): input set which determines dictionary content

    Returns:
        dictionary containing all relevant OPEN file handles
    """
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
    """close all open file handles in dictionary

    Returns:
        None
    """
    for fh in ALL_FILE_HANDLES.values():
        fh.close()

def create_extended_stats_dict():
    """perform statistics of %mapped to every barcode

    Returns:
        dictionary: containing all of our barcodes, which we will populate with 
        our counts as we run through the files.
    """
    #all hopped barcodes
    l=list(itertools.permutations(list(VALID_INDICES_SET),2))
    #all matched barcodes
    for item in VALID_INDICES_SET:
        l.append((item,item))
    #all idk barcodes
    l.append("unk")
    #set default value to zero 
    return {key: 0 for key in l}

def perform_extended_stats(stats_dict:dict, stats_list:list):
    """create our extended statistics to populate extended_stats.txt

    Args:
        stats_dict (dict): dictionary with detailed barcode stats
        stats_list (list): list with categorical stats
    """
    #sort by value (so top matches on top of file)
    stats_dict={k: v for k, v in sorted(stats_dict.items(), 
                                        key=lambda item: item[1], reverse=True)}
    #total number of records found
    total_num_records=sum(list(stats_dict.values()))

    #create percentages formatted 1.23%
    hopped_perc="%.4f" % ((stats_list[2]/total_num_records)*100)
    unknown_perc="%.4f" % ((stats_list[1]/total_num_records)*100)
    matched_perc="%.4f" % ((stats_list[0]/total_num_records)*100)

    #write major record categories
    stats_log=open("outfiles/extended_stats.txt","w")
    stats_log.write(f"Total Number of Records: {sum(stats_list)}\n\nRecord breakdown by category:\n")
    stats_log.write(f"Unknown Records: {unknown_perc}% ({stats_list[1]})\nHopped Records: {hopped_perc}% ({stats_list[2]})\nMatched Records: {matched_perc}% ({stats_list[0]})\n\n")
    stats_log.write(f"Detailed record breakdown:\n") 

    #show deatiled breakdown of categories 
    for barcodes,counts in stats_dict.items():
        percentage="%.2f" % ((counts/total_num_records)*100)
        #only separating for formatting reasons
        if barcodes == "unk":
            stats_log.write(f"Unknown records: {percentage}% ({counts})\n")

        else:
            stats_log.write(f"{tuple(barcodes)[0]} {tuple(barcodes)[1]} records: {percentage}% ({counts})\n")
    stats_log.close()

def graph_barcodes(file_name="outfiles/extended_stats.txt"):
    """perform optional graphing for barcodes stats

    Args:
        file_name (str, optional): file with our graphing information.
        Defaults="outfiles/extended_stats.txt".

    Returns:
        None: 
    """
    fh=open(file_name,"r")
    #skip general counts and header
    relevant=False
    #for % mapped
    mapped_counts=[]
    hopped_counts=[]
    #for barcode
    mapped_cat=[]
    hopped_cat=[]
    for line in fh:
        if "Detailed" in line:
            relevant=True
            continue
        if relevant:
            #only include barcode barcode
            theo_cat=line[:line.index("rec")].strip()
            #split the list to ensure its the matched barcode 
            l=line.split()
            #if these are two matched barcodes, then grab the floats
            if l[0]==l[1]:
                mapped_cat.insert(0,theo_cat)
                #grab % mapped from the list 
                #change to -2 for final file
                mapped_counts.insert(0,float(line.split()[-2][:-1]))
            elif l[0]!=l[1] and float(line.split()[-1][1:-1]) > 0 and l[0]!="Unknown":
                #want top values last so use insert
                hopped_cat.append(theo_cat)
                hopped_counts.append(float(line.split()[-1][1:-1]))
    #MAPPED
    #so barcodes are not cut off
    plt.rcParams.update({'font.size': 20, "figure.figsize" : (30,20)})
    #so I can remove my axes
    fig,ax=plt.subplots()
    ax.set(title="% mapped, by barcode",ylabel="Barcode",xlabel="% mapped")
    plt.barh(mapped_cat,mapped_counts,alpha=0.8,color="mediumturquoise")
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False) 
    #place with the other outfiles
    plt.savefig("outfiles/percent_mapped_by_barcode.png")

    #HOPPED
    plt.rcParams.update({'font.size': 20, "figure.figsize" : (30,20)})
    #so I can remove my axes
    fig,ax=plt.subplots()
    ax.set(title="hopped counts, by barcode",ylabel="Barcode",xlabel="number mapped")
    #too many hopped, only get the top 20
    plt.barh(hopped_cat[0:20][::-1],hopped_counts[0:20][::-1],alpha=0.8,color="mediumturquoise")
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False) 
    #place with the other outfiles
    plt.savefig("outfiles/hopped_counts_by_barcode.png")

def demultiplex(R1, I1, R2,  I2):
    """perform demultiplexing on our four files

    Args:
        R1 (str): File, R1 first biological files
        I1 (str): file R2, first index files
        R2 (str): File R4, second biological file
        I2 (str): File R2, second index file.

    Returns:
        stats_info: dictionary of statistics for per-barcode counts.
    """
    error_log=open("outfiles/error_log.txt","w")

    #tracks relevant header info, one record at a time
    header_info={"r1":"","i1":"","r2":"","i2":""}
    #tracks relevant seq info, one record at a time
    sequence_info={"r1":"","i1":"","r2":"","i2":""}
    #tracks global barcode counts (all possible permuations)
    stats_info=create_extended_stats_dict()
    #track barcode counts by category
    unknown_count=0
    barcode_count=0
    hopped_count=0

    #open files
    with gzip.open(R1,"rt") as r1, gzip.open(I1,"rt") as i1,\
         gzip.open(R2,"rt") as r2, gzip.open(I2,"rt") as i2:
        #going line by line (rn=biological record, in=index record)
        index=0
        for line_r1, line_r2, line_i1, line_i2 in zip(r1,r2, i1, i2):
            #new record, write header info
            if index%4==0:
                
                header_info={"r1":line_r1.strip(), "r2":line_r2.strip()}

            #new record, write sequence info
            if index%4==1:
                line_i2=reverse_complement(line_i2) 
                sequence_info={"r1":line_r1.strip(), "r2":line_r2.strip(),
                               "i1":line_i1.strip(), "i2":line_i2.strip()}

            #if quality line
            if index%4==3:
                i1_quality=calculate_per_base_quality(line_i1.strip())
                i2_quality=calculate_per_base_quality(line_i2.strip())
                i2_quality.reverse() #because we rc the sequence line

                #create fastq records, then decide where we will write them
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

                #case 1 indices match
                if verify_index_match(sequence_info["i1"], sequence_info["i2"]): 
                    #if indices match, exist and are HQ
                    if search_for_index_match(sequence_info["i1"])==True and\
                        check_quality_thresholds(i2_quality) and\
                        check_quality_thresholds(i1_quality):

                        stats_info[(sequence_info["i1"],sequence_info["i1"])]+=1
                        barcode_count+=1
                        ALL_FILE_HANDLES[f"{sequence_info["i1"]}_R1"].write(fastq_record1)
                        ALL_FILE_HANDLES[f"{sequence_info["i2"]}_R2"].write(fastq_record2)

                    #indices both match an unknown barcode, or are low quality
                    else: 
                        stats_info["unk"]+=1
                        ALL_FILE_HANDLES["unknown_R1"].write(fastq_record1)
                        ALL_FILE_HANDLES["unknown_R2"].write(fastq_record2)
                        unknown_count+=1

                #indices dont match
                elif sequence_info["i1"] != sequence_info["i2"]:                      
                    #for a barcode to hop both need to exist and be hq
                    if search_for_index_match(sequence_info["i1"]) and\
                        search_for_index_match(sequence_info["i2"]) and\
                        check_quality_thresholds(i1_quality) and\
                        check_quality_thresholds(i2_quality):

                        stats_info[(sequence_info["i1"],sequence_info["i2"])]+=1
                        ALL_FILE_HANDLES["hopped_R1"].write(fastq_record1)
                        ALL_FILE_HANDLES["hopped_R2"].write(fastq_record2) 
                        hopped_count+=1

                    #otherwise UNKNOWN
                    else:
                        stats_info["unk"]+=1
                        unknown_count+=1
                        ALL_FILE_HANDLES["unknown_R1"].write(fastq_record1)
                        ALL_FILE_HANDLES["unknown_R2"].write(fastq_record2) 

                #if our logic somehow missed a case
                else:
                    print("TRASH FOUND CHECK ERROR_LOG")
                    error_log.write(fastq_record1)
                    error_log.write(fastq_record2)

            index+=1
    error_log.close()
    return perform_extended_stats(stats_info,[barcode_count,unknown_count,hopped_count])

#read in ...indexes.txt
create_index_set()
#open all possible file combinations
ALL_FILE_HANDLES=open_files_of_interest(VALID_INDICES_SET)
#actually perform demultiplexing 
stats_info=demultiplex(R1=args.r1,I1=args.r2,I2=args.r3,R2=args.r4)
#close completed files of interest
close_files_of_interest()
#optional graphing
graph_barcodes()

# small test files
# stats_info=demultiplex(R1="my_test_files/r1_mid_test.fq",
#                      I2="my_test_files/r3_mid_test.fq",
#                       I1="my_test_files/r2_mid_test.fq",
#                       R2="my_test_files/r4_mid_test.fq")
