import os
import numpy as np
import pickle

def Seq_converter_PSSM(windowsize,filename):
    with open(filename, 'rb') as fp:
        #load PSSM data one sample at a time
        data = pickle.load(fp) 
        #Add head and tail to build a window
        halfwindow = int((windowsize-1)/2)    
        for i in range(halfwindow):
            data.append(data[-1])
            data.insert(windowsize,data[0])
        
        #Use slide window
        windowscale = int((windowsize-1)/2)
        seqLength = len(data)
        Seq = []
        
        #for i in range(len(data)):
         #   seqLength = len(data)
        for j in range(windowscale,seqLength-windowscale):
            windowSeq = data[j-windowscale:j+windowscale+1]
            windowSeq = np.hstack(windowSeq)
            Seq.append(windowSeq)
        Seq = np.asarray(Seq)
        return Seq
        
def UsePSSMDirectory(windowsize,directory):
    PSSMcombinedWindow = []
    for filename in os.listdir(directory):
        if filename.endswith(".pssmparsed"):
            a = Seq_converter_PSSM(windowsize,filename)
            PSSMcombinedWindow.append(a)
    flat_PSSM = [item for sublist in PSSMcombinedWindow for item in sublist]   
    return flat_PSSM

def UsePSSM(windowsize):
    path = os.getcwd()
    path2 = path+'\\psiblast_parsed'
    os.chdir(path2)
    return UsePSSMDirectory(windowsize,path2)
