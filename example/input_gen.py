
#!/usr/bin/env python3

import argparse
import os
import sys
import random as rnd

# Default params ---------------------------------------------------------------

NSEQ = 1

# Input params -----------------------------------------------------------------

def params():

    obj_parse_args = argparse.ArgumentParser(description="""Creates
        couple of random sequences and prints out in a file.""")
    
    obj_parse_args.add_argument("-l", "--len",
                                type = int,
                                required = True,
                                help = "Length of sequences")
    obj_parse_args.add_argument("-n", "--nseq",
                                type = int,
                                required = False,
                                default = NSEQ,
                                help = "Number of couple of sequences")

    args = obj_parse_args.parse_args()

    return args

# Parameters check -------------------------------------------------------------

def checkParams(parameters):

    if parameters.len < 0:
        print("Invalid parameter for len, neative len no accepted!")
        sys.exit(1)
        
    if parameters.nseq is None:
        print("""No parameter for --nseq,
                default will  be used instead: {}""").format(NSEQ)
        parameters.nseq = NSEQ

    else:
        if parameters.nseq < 0:
            print("""Invalid parameter for len,
                default will be used instead: {}""").format(NSEQ)
            parameters.nseq = NSEQ
        
    return parameters

# Main -------------------------------------------------------------------------

def SeqsCreator(seq_len):
    bases = ["A","C","T","G"]
    tmp_seq = ""
    out_seqs = []
    
    for i in range(2):
        for i in range(seq_len):
            tmp_base = rnd.choice(bases)
            tmp_seq += tmp_base
        out_seqs.append(tmp_seq)
        tmp_seq = ""

    return out_seqs
            
def PrintOut(out_seqs, out):
    fh = open(out, "a+")
    fh.write("\t".join(out_seqs))
    fh.close()

# Program caller ---------------------------------------------------------------

if __name__ == "__main__":
    args = params()
    args = checkParams(args)
    
    for i in range(args.nseq):
        out_seqs = SeqsCreator(args.len)
        if i > 0:
            out_seqs[0] = "\n" + out_seqs[0]
            
        PrintOut(out_seqs, "./random_seqs.txt")
