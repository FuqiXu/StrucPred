import os
import numpy as np
import pandas as pd
from pandas.core.frame import DataFrame

path = os.getcwd()

def test_fasta_read(filename):
# Read fasta file
    seqID,seq = [],[]
    with open(filename) as f:
        data = f.read().splitlines()
        for i in range(0,len(data),2):
            if(data[i][0]=='>'):
                seqID.append(data[i])
                seq.append(data[i+1])
    seqData = DataFrame({
            "seqID":seqID,
            "seq":seq
            })
    
    for i in range(len(seqData.seq)):
        a = list(seqData.seq[i])
        seqData.seq[i]=a
    return seqData

def seq_converter(seq):
    # Vectorize data
    aa_dic = {  'A':[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],
                'R':[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],
                'N':[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],
                'D':[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],
                'C':[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],
                'Q':[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],
                'E':[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,],
                'G':[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,],
                'H':[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,],
                'I':[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,],
                'L':[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,],
                'K':[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,],
                'M':[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,],
                'F':[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,],
                'P':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,],
                'S':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,],
                'T':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,],
                'Y':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,],
                'W':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,],
                'V':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,],}
    
    for i in range(len(seq)):
            for j in range(len(seq[i])):
                if seq[i][j] in aa_dic:
                    seq[i][j] = aa_dic.get(seq[i][j])
                else:
                    seq[i][j] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,]
    return seq

def test_vectors(filename):
    data = test_fasta_read(filename)
    # Converting residues from letters into numbers
    seq_converter(data.seq)
    return data

def test_window(windowsize,seqData):
    # Add slide window to evaluate the environment's impact on topology
    for i in range(len(seqData)):
        seqFirst=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,]     
        seqLast=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,]
        halfwin = int((windowsize-1)/2)
        for j in range(halfwin):
            seqData.seq[i].append(seqLast)
            seqData.seq[i].insert(0,seqFirst)
           
    # Creating a slide window.The basic element in one window is #windowsize*AA   
    for m in range(len(seqData)):
        seq_single = []
        len_m = len(seqData.seq[m])-windowsize+1
        for p in range(len_m):
            temp = []
            for n in range(windowsize):
                temp.extend(seqData.seq[m][p+n])
            seq_single.append(temp)    
        seqData.seq[m]=seq_single
    return seqData

def test_parser(filename,windowsize):
    testSeq = []
    for i in range(len(test_fasta_read(filename).seq)):
        testSeqSingle=''.join(test_fasta_read(filename).seq[i])
        testSeq.append(testSeqSingle) 
    
    testBinary = test_vectors(filename)
    testWind = test_window(windowsize,testBinary)
    testData = []
    for i in range(len(testWind)):
        seqData = testWind.iloc[i]
        testData.append(seqData)    
    return testSeq,testData




