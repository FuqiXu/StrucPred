import os
import numpy as np
import pandas as pd
from pandas.core.frame import DataFrame
import pickle
from sklearn import svm

from predictor import binary_rawdata
from predictor import data_window
from predictor import data_svm

path = os.getcwd()

###### Parsing single train sequence ######

def pssm_parser(seqfile,pssmfile,windowsize):
    # Read a single sequence file
    seqID1, seq1= [], []
    with open(seqfile) as f:
        data = f.read().splitlines()
        for i in range(len(data)):
            if i%3 == 1:
               seq1.append(data[i])
            if i%3 == 0:
               seqID1.append(data[i])
             
    seqData1 = {
            "seqID":seqID1,
            "seq":np.array(seq1),
            }
    seqData = DataFrame(seqData1)
    
    for i in range(len(seqData.seq)):
        a = list(seqData.seq[i])
        seqData.seq[i]=a
    
    # For each single sequence file, vercterize amino acid.    
    pssmcsv = pd.read_csv(pssmfile, header=None)
    pssm = pssmcsv.transpose()
    
    data=[]
    col=pssm.shape[1]
    for i in range(col):
         data.append(pd.Series.tolist(pssm[i]))
    
    #adding head and tails to build a window
    seqFirst=data[0]      
    seqLast=data[-1]
    halfwin = int((windowsize-1)/2)
    for i in range(halfwin):
            data.append(seqLast)
            data.insert(0,seqFirst)
        
    # Creating a slide window
    seq= []
    for i in range(len(data)-2*halfwin):
        temp = []
        for n in range(windowsize):
            temp.extend(data[i+n])
        #print(temp)
        #print('\n')
        seq.append(temp)    

    return seq

def mult_seq_pssm_parser(seqDirectory,pssmDirectory,windowsize):
    dataSeq = []
    for seqfilename in os.listdir(seqDirectory):
        seqfilepath = os.path.join(seqDirectory, seqfilename)
        pssmfilepath = os.path.join(pssmDirectory,seqfilename+'.pssm.csv')
        dataSeq.extend(pssm_parser(seqfilepath,pssmfilepath,windowsize))
    return dataSeq

if __name__ == "__main__":
    windowsize = 15  
    
    print("Parsing data...")
    dataBinary = binary_rawdata("data/trainset.dat")
    
    print("Adding window...")
    dataWind = data_window(windowsize,dataBinary)
    
    print("SVM prediction preparing...")
    dataSVM = data_svm(dataWind)
    dataSeq = mult_seq_pssm_parser('pssm/Sequences',
                                   'pssm/pssmMatrix',windowsize)
    dataStruc = pd.Series.tolist(dataSVM.seqTopo)

    print("Model building...")
    clf = svm.LinearSVC(max_iter = 500, dual = False)
    clf.fit(dataSeq,dataStruc)
	
    print("Saving models...")
    filepath = os.path.join('models', 'linsvm_pssm.pkl')
    if not os.path.exists('models'):
        os.makedirs('models')
    with open(filepath,'wb')as f:
	    pickle.dump(clf,f)
	
    print("Model Built!")
    pass
