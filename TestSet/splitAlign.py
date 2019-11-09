#!/usr/bin/env python

"""
A script to take an alignment and split it into two halves. Goal is to generate positive control data for protein interaction detection.

:Arguments:
Input_MSA.fasta (the alignment to be split)

:Keyword Arguments:
--output, -o     output prefix

:Example:
>>> ./splitAlign.py  eggNOG_aligns/nrdA_trimmed_ENOG4105BZH.fasta -o nrdA_split

:By: Kim Reynolds
:On: 10/3/2019

"""
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment", help='Input Sequence Alignment')
    parser.add_argument("-o","--output", dest = "outPrefix", help="Output File Prefix")
    options = parser.parse_args()

    if options.outPrefix is None:
        outPrefix = options.alignment[0:4]
    else:
        outPrefix = options.outPrefix

    #read in alignment and determine length
    filelines = open(options.alignment, 'rb').readlines()
    headers = list(); sequences = list(); notfirst = 0
    for line in filelines:
        if line[0] == '>':
            if notfirst > 0: sequences.append(seq.replace('\n','').upper())
            headers.append(line[1:].replace('\n',''))
            seq = ''; notfirst = 1
        elif line != '\n': seq += line
    sequences.append(seq.replace('\n','').upper())


    #split alignment and write out
    Npos = len(sequences[0])
    splitSite = int(Npos/2)

    f1 = open(outPrefix+'_half1.fasta', 'wb')
    f2 = open(outPrefix+'_half2.fasta', 'wb')

    for i,hd in enumerate(headers):
        f1.write('>'+hd+'\n')
        f1.write(sequences[i][0:splitSite]+'\n')
        f2.write('>'+hd+'\n')
        f2.write(sequences[i][splitSite:]+'\n')

    f1.close()
    f2.close()
