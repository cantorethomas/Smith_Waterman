#!/usr/bin/env python

# Sequence aligner based on Smith Waterman algorithm. 
import argparse
import os
import sys

# Default Parameters -------------------------------------------------------------------------------

MATCH = +3
MISMATCH = -3
GAP = -2

# Input parameters ---------------------------------------------------------------------------------

def params():
    obj_parse_args = argparse.ArgumentParser(description="""Given two sequences
        or a file containing sequences return the alignments.
        The sequences in the input file have to be tab separated""")
    obj_parse_args.add_argument("--seq1",
                                type = str,
                                help = "The first sequence to align")
    obj_parse_args.add_argument("--seq2",
                                type = str,
                                help = "The second sequence to align")
    obj_parse_args.add_argument("-i", "--input",
                                type = str,
                                help = "The file given in input")
    obj_parse_args.add_argument("-m", "--match",
                                default = MATCH,
                                type = float,
                                required = False,
                                help = "Score of the match (default = {})".format(MATCH))
    obj_parse_args.add_argument("-s", "--mismatch",
                                default = MISMATCH,
                                type = float,
                                required = False,
                                help = "Score of the mismatch (default = {})".format(MISMATCH))
    obj_parse_args.add_argument("-g", "--gap",
                                default = GAP,
                                type = float,
                                required = False,
                                help = "Score of the gap (default = {})".format(GAP))
    obj_parse_args.add_argument("--minscore",
                                type = float,
                                help = "The minimum score")
    obj_parse_args.add_argument("--minlength",
                                type = float,
                                help = "The minimum length")
    obj_parse_args.add_argument("--numresult",
                                type = int,
                                help = "The number of alignments to return")

    args = obj_parse_args.parse_args()
    return args
                                
# Check for Parameters -----------------------------------------------------------------------------                              
                           
def check_params(parameters):
    parameters.seqs = []
    input_f = False
    if parameters.input is not None:
        if os.path.exists(parameters.input):
            if os.path.isfile(parameters.input):
                size = len(parameters.seqs)
                parameters.seqs += [tuple(line.strip().split("\t")) for line
                    in open(parameters.input) if len(line.strip().split("\t"))==2]
                input_f = True
                if len(parameters.seqs) <= size:
                    print("""the file could be empty or wrongly formatted: \nseq1\tseq2
                        \nseq3\tseq4\n..\n(tab separated)""")
            else:
                print("input is not a file")
        else:
            print("input file doesn't exist!")

    if ((parameters.seq1 is not None) and (parameters.seq2 is not None)):
        parameters.seqs.append((parameters.seq1, parameters.seq2))

    if ((parameters.seq1 is None) or (parameters.seq2 is None)) and input_f is False:
        print("You have to give in input at least two sequences or a file of sequences!")
        if not parameters.seqs:
            print("You did not give me anything...")
            sys.exit(1)

    if parameters.match <= 0:
        print("provided value for matches is not correctm default will be used instead")
        parameters.match = MATCH

    if parameters.mismatch > 0:
        print("provided value for mismatches is not correct, default will be used instead!")
        parameters.mismatch = MISMATCH

    if parameters.gap > 0:
        print("provided value for gap is not correct, default will be used instead")
        parameters.gap = GAP

    if parameters.minscore is not None:
        if parameters.minscore < 0:
            print("provided value for minimum score is negative, no minimum score will be used")
            parameters.minscore = None

    if parameters.minlength is not None:
        if parameters.minlength < 0:
            print("provided value for minimum length is negative, no minimum length will be used")
            parameters.minlength = None

    return parameters
            
# Accessory Functions ----------------------------------------------------------

# Find max #
def findMax(matrix):
    """ Finds max in score matrix, if more than one,
        stores them in a list"""
    max_so_far = 0
    pos = []
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if matrix[i][j] > max_so_far:
                max_so_far = matrix[i][j]
                pos = [i,j]
                
    return pos, max_so_far

# Score dictionary #
def CreateSD(matrix):
    """ Creates a dictionary with all the scores and
        relative positions inside of the matrix """
    sd = {}
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if matrix[i][j] != 0:
                if matrix[i][j] not in sd:
                    sd[matrix[i][j]] = [[i,j]]
                else:
                    sd[matrix[i][j]].append([i,j])

    scores = sorted(sd, reverse = True)
    return sd, scores
                
# Check Events #
def checkEvents(pos, matrix):
    """ checks the max val for gap or mm from current pos.
        Stores data in a composed by three vals:
        1) Score inside of matrix
        2) Position considered
        3) Type of event """
    
    g1 = [matrix[pos[0]-1][pos[1]], [pos[0]-1, pos[1]], "G1"]
    g2 = [matrix[pos[0]][pos[1]-1], [pos[0], pos[1]-1], "G2"]
    mm = [matrix[pos[0]-1][pos[1]-1], [pos[0]-1, pos[1]-1], "MI"]

    dict_moves = {g1[0]:g1[1:], g2[0]:g2[1:], mm[0]:mm[1:]}
    max_move = max(dict_moves) # if want to implement bifurcation, start here

    next_pos = dict_moves[max_move][0]
    event = dict_moves[max_move][1]
    
    return next_pos, event

# Print out #
def printOut(stru, strd, events,i):
    print("> Alignment {}\n".format(i+1))
    print("{}\n{}\n{}\n".format(stru, "".join(events), strd))

# Matrix Construction ------------------------------------------------------------------------------

def ConstructMatrix(str1, str2, M, MM, G):
    """Constructing score matrix with M, MM, GG as penalties
    given in input. str1, str2 two sequence to align."""

    ls1 = len(str1)
    ls2 = len(str2)
    matrix = [[0]*(ls2+1) for i in range(ls1+1)]

    # filling the matrix
    
    for i in range(1,ls1+1): 
        for j in range(1,ls2+1):
            
            if str1[i-1] == str2[j-1]:
                matrix[i][j] = (matrix[i-1][j-1] + M)
                
            else:
                g1 = matrix[i-1][j] + G
                g2 = matrix[i][j-1] + G
                mm = matrix[i-1][j-1] + MM
                matrix[i][j] = max(g1,g2,mm,0)

    return matrix
    
# Trace back of Matrix -----------------------------------------------------------------------------

def TraceBack(matrix, str1, str2, pos):
    """ Trace back function from max val to 0: prints all possible
        aligments, considering every different starting point """
        
    events = ""
    stru = ""
    strd = ""

    pos_val = matrix[pos[0]][pos[1]]
    
    while pos_val != 0:
        base1 = str1[pos[0]-1]
        base2 = str2[pos[1]-1]

        if base1 == base2: #If same base, calls for match 
            events = "|" + events
            stru = base1 + stru
            strd = base2 + strd

            pos = [pos[0]-1, pos[1]-1]

        else:
            pos, tmp_event = checkEvents(pos, matrix)

            if tmp_event == "G1":
                stru = base1 + stru
                strd = "-" + strd
                events = " " + events
                
            elif tmp_event == "G2":
                stru = "-" + stru
                strd = base2 + strd
                events = " " + events
                
            else:
                stru = base1 + stru
                strd = base2 + strd
                events = ":" + events
        
        pos_val = matrix[pos[0]][pos[1]]

    return(stru, strd, events)

# Program Call -------------------------------------------------------------------------------------

if __name__ == "__main__":
    args = params()
    args = check_params(args)
    for seqs in args.seqs:
        str1 = seqs[0]
        str2 = seqs[1]
        
        matrix = ConstructMatrix(str1, str2,
                                 args.match,
                                 args.mismatch,
                                 args.gap)


        if args.numresult == None:
            pos, score = findMax(matrix)
            
            if ((args.minscore == None) or (args.minscore <= score)):
                stru, strd, events = TraceBack(matrix, str1, str2, pos)
                
                if ((args.minlength == None) or (args.minlength <= len(stru))):
                    print("SCORE: {}\n{}\n{}\n{}\n\n".format(score, stru, "".join(events), strd))

                else:
                    print("Resulting alignment length < minlength value: {} < {}".format(len(stru),
                            args.minlength))

            else:
                print("Resulting alignment score < minscore value: {} < {}".format(score,
                            args.minscore))
        
        else:

            score_dict, scores = CreateSD(matrix)
            count = args.numresult
            positions = [] #list of numresult positions

            i = 0
            j = 0
            while count != 0:
                positions.append(score_dict[scores[i]][j])
                count -= 1
                if len(score_dict[scores[i]]) > (j+1):
                    j += 1
                else:
                    if len(score_dict) > (i+1):
                        i += 1
                        j = 0

            tot_print = 0 # Number of alignment that passed the constrictions
            
            for pos in positions:
                score = matrix[pos[0]][pos[1]]

                if ((args.minscore == None) or (args.minscore <= score)):
                    stru, strd, events = TraceBack(matrix, str1, str2, pos)
                
                    if ((args.minlength == None) or (args.minlength <= len(stru))):
                        print("\nSCORE: {}\n{}\n{}\n{}".format(score, stru, "".join(events), strd))
                        tot_print += 1

# SOME NOTES FOR FUTURE USES #
# The algorithm was planned to work divided mainly into Score Matrix creation and Trace Back.
# If in input just one sequence is required, the accessory function findMax() is used. Else,
# A score_dictionary object is created through specific function and only the higher requsted
# positions are used to call iteratively traceback. The last part works on calling the program
# depending on which restrictions are used and print out (restrictions: --minlength, --minscore)
                
                
                
                
                    
            
    
                

        
        
