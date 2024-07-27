#!/usr/bin/env python

# Author: <Julia Jones> <jujo@uoregon.edu>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.4"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = None
RNA_bases = None

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter)-33

def qual_score(phred_score: str) -> float:
    """Convert a phred33 encoded string to a list of quality values 
    then compute the average of those values"""
    scores=[]
    for item in phred_score:
        scores.append(convert_phred(item))
    return sum(scores)/len(scores)

def validate_base_seq():
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    pass

def gc_content():
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    pass

def calc_median(lst):
    '''Given a sorted list, returns the median value of the list'''
    #if empty list, nothing to calculate
    if len(lst)==0:
        return None
    #find middle point, _f vs _i determines odd vs even list
    mid_point_f=(len(lst)-1)/2
    mid_point_i=int(mid_point_f)
    
    #if list is an even number, calc avg between the two mid points for median
    if mid_point_i<mid_point_f:
        return (lst[mid_point_i]+lst[mid_point_i+1])/2
    #otherwise just return item at mid_point
    return lst[mid_point_i]
        

def oneline_fasta(input_file:str):
    """remove new lines from sequence data in fasta record

    Args:
        input_file (str): name of file
        dir (str, optional): path to file (needed bc where my talapas terminal
        is opened). Default="PS/ps6-jemdhb/".
    """
    opened_file=open(input_file,"r")
    #create unique name
    short_name=input_file[:input_file.rfind(".")]
    output_file_name=short_name+"_collapsed.fa"
    output_file=open(output_file_name,"w")
    prev=""
    for line in opened_file:
        curr=line
        #if not near a header, strip newlines
        if ">" not in curr and ">" not in prev:
            prev=prev.strip()
        output_file.write(prev)
        prev=curr
    return output_file_name

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")
