# -*- coding: utf-8 -*-
"""
Prepare sequence-structure dataset to do blast
Input: One file with multiple sequences and structure data, 'cas2.3line.txt'
Output: Files named "sequence#N.fasta", each file contains 2 lines, the first line is seuqnece id, the second line is the sequences.
Process: read raw_input ->save the first two lines ->new file

"""
def blastParserRaw(filename):
    import os
    with open(filename) as f:
        data = f.read().splitlines()
    path = "/BioinfoProject/KB8025/dataset/SequencesToBlast"
    os.makedirs(path,exist_ok=True);
    os.chdir(path)
    for i in range(0,len(data),3):
        a = str(int(i/3))
        result = open("Sequence"+a+".fasta","w")
        result.write(data[i])
        result.write("\n")
        result.write(data[i+1])
        result.write("\n")
        result.close()
blastParserRaw('/BioinfoProject/KB8025/dataset/cas2.3line.txt')
