# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 20:32:51 2018

@author: Fuqi Xu
"""
from pandas.core.frame import DataFrame
import numpy as np
import os
import pandas
from sklearn.model_selection import cross_validate
from sklearn import svm
from sklearn.externals import joblib

path = os.getcwd()

###### Parsing train dataset######

# Read data from 3 line fasta file and store them in a data frame
def rawtoframe(filename):
    seqID1, seq1, seqTopo1= [], [], []
    with open(filename) as f:
        data = f.read().splitlines()
        for i in range(len(data)):
            if i%3 == 1:
               seq1.append(data[i])
            if i%3 == 2:
               seqTopo1.append(data[i])
            if i%3 == 0:
               seqID1.append(data[i])
              
    seqData1 = {
            "seqID":seqID1,
            "seq":seq1,
            "seqTopo":seqTopo1
            }
    seqData = DataFrame(seqData1)

    # Convert every sequence and sequence topology from list to arrays
    for i in range(len(seqData.seq)):
        a = list(seqData.seq[i])
        seqData.seq[i]=a
    for i in range(len(seqData.seqTopo)):
        a = list(seqData.seqTopo[i])
        seqData.seqTopo[i]=a
        
    return seqData

# Vectorising data.
def seq_converter(seq):
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
                'V':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,]}

    for i in range(len(seq)):
        for j in range(len(seq[i])):
            if seq[i][j] in aa_dic:
                seq[i][j] = aa_dic.get(seq[i][j])
                
    return (seq)

def topo_converter(seqTopo):
    for i in range(len(seqTopo)):
        for j in range(len(seqTopo[i])):
            if seqTopo[i][j]=='H':
                seqTopo[i][j] = 0
            if seqTopo[i][j]=='E':
                seqTopo[i][j] = 1
            if seqTopo[i][j]=='C':
                seqTopo[i][j] = 2
    
    return seqTopo

# Read raw data into a dataframe and convert them into vectors.
def binary_rawdata(filename):
    data = rawtoframe(filename)
    
    # Converting residues from letters into numbers
    seq_converter(data.seq)
    topo_converter(data.seqTopo)
    return data


########### adding winodws #############

# Add slide window to evaluate the environment's impact on topology
def data_window(windowsize,data):
    # Adding head and tails in protein sequence data.
    for i in range(len(data)):
        seqFirst=data.seq[i][0]      
        seqLast=data.seq[i][-1]
        halfwin = int((windowsize-1)/2)
        for j in range(halfwin):
            data.seq[i].append(seqLast)
            data.seq[i].insert(0,seqFirst)
        
    # Creating a slide window.The basic element in one window is #windowsize*AA   
    for m in range(len(data)):
        seq_single = []
        for p in range(len(data.seqTopo[m])):
            temp = []
            for n in range(windowsize):
                temp.extend(data.seq[m][p+n])
            seq_single.append(temp)    
        data.seq[m]=seq_single
        
    return data

# Transfering data into a binary array to be used in svm
def data_svm(data):
    sequence = []
    structure = []
    data.seq = np.array(data.seq)
    data.seqTopo = np.array(data.seqTopo)
    for i in range(len(data)):
        for j in range(len(data.seq[i])):
            sequence.append(data.seq[i][j])
        for k in range(len(data.seqTopo[i])):
            structure.append(data.seqTopo[i][k])
    dataSVM = DataFrame({
            'seq':sequence,
            'seqTopo':structure
                           })
    return dataSVM


### Test file parser（3 line fasta format) ###
# Store each sequence infomation in a dataframe, and keep them in a list
def test_fasta(filename,windowsize):
    # Storing original format of sequences.
    testSeq = []
    for i in range(len(rawtoframe(filename).seq)):
        testSeqSingle=''.join(rawtoframe(filename).seq[i])
        testSeq.append(testSeqSingle) 
    
    realStruc= topo_converter(rawtoframe(filename).seqTopo)
    
    testBinary = binary_rawdata(filename)
    testWind = data_window(windowsize,testBinary)
    testData = []
    for i in range(len(testWind)):
        seqData = testWind.iloc[i]
        testData.append(seqData)
        
    return testSeq,testData,realStruc


### Save prediction result ###
def sav_pred(prediction,testdata,testSeq,testfilename):
    # Create a prediction result folder a .dat file"
    filepath = os.path.join('result', 'pred.dat')
    if not os.path.exists('result'):
        os.makedirs('result')
    f = open(filepath, "a")
    
    # Change numberic prediction to letters
    for i in range(len(prediction)):
        # Save prediction result to file
        f.write(testdata[i].seqID)
        f.write("\n")
        f.write(testSeq[i])
        f.write("\n")
        
        predStruc = []
        pred = []
        for j in range(len(prediction[i])):
            if prediction[i][j]==0:
                predStruc.append ('H')
            if prediction[i][j]==1:
                predStruc.append ('E')
            if prediction[i][j]==2:
                predStruc.append ('C')
        pred = str.join("",predStruc)
        f.write(pred)
        f.write("\n")
    f.close()


### evaluation ###
def performance(pred,real):
    import math
    f = open("result/evaluation.dat",'a')
       
    # Q3
    correct = 0
    total = 0
    for i in range(len(pred)):
        for j in range((len(pred[i]))):
            total+=1
            if pred[i][j]==real[i][j]:
                correct+=1
    q3 = correct/total
    f.write("Q3:"+str(format(q3, '.0%')))
    f.write('\n')
    
    # Q(x); 0-H;1-E;2-C
    for x in [0,1,2]:
        correctx = 0
        totalx = 0
        for i in range(len(pred)):
            for j in range((len(pred[i]))):
                if real[i][j]==x:
                    totalx+=1
                    if pred[i][j]==real[i][j]:
                        correctx+=1
        qx = correctx/totalx
        
        convert = {0:'H',1:'E',2:'C'}
        if x in convert:
            rep = convert[x] 
        
        f.write("Q("+rep+'):'+str(format(qx, '.0%')))  
        f.write('\n')
        
    # Corrlation coefficient(C(H),C(E),C(C))
    structures = [0,1,2]
    for x in structures:
        realx = 0
        predx = 0
        prednotx = 0
        Nopredx=0
        NopredNotx=0
        totalx = 0
        for i in range(len(pred)):
            for j in range((len(pred[i]))):
                if real[i][j]==x:
                    realx+=1
                    if pred[i][j]==x:
                        predx +=1
                    else:
                        prednotx+=1       
                if real[i][j]!=x:
                    if pred[i][j]==x:
                        Nopredx+=1
                    if pred[i][j]!=x:
                        NopredNotx+=1
        Px = predx
        Rx = NopredNotx
        Ux = prednotx
        Ox = Nopredx
        
        Cx = format((Px*Rx-Ux*Ox)/
                    (math.sqrt((Px+Ux)*(Px+Ox)*(Rx+Ux)*(Rx+Ox))), '.0%')
        
        convert = {0:'H',1:'E',2:'C'}
        if x in convert:
            rep = convert[x] 
        f.write("C("+rep+'):'+str(Cx)) 
        f.write('\n')

    f.close()
### pssm ###

### novel sequence parser ###  

def use(windowsize,trainfile,testfile):
    print("Parsing data...")
    dataBinary = binary_rawdata(trainfile)
    
    print("Adding window...")
    dataWind = data_window(3,dataBinary)
    
    print("SVM prediction preparing...")
    dataSVM = data_svm(dataWind)
    dataSeq = pandas.Series.tolist(dataSVM.seq)
    dataStruc = pandas.Series.tolist(dataSVM.seqTopo)
      
    print("Model building...")
    clf = svm.SVC(kernel='linear', C=1, random_state=0)
    clf.fit(dataSeq,dataStruc)
    
    print("Preparing test data...")
    testSeq,testData,realStruc = test_fasta(testfile,3)  
    
    print("Predicting...")
    preds = []
    for i in range(len(testData)):
        pred = clf.predict(testData[i].seq)
        preds.append(pred)  
    
    print("Saving prediction...")
    sav_pred(preds,testData,testSeq,testfile)
    
    print("Cross validating...")
    scoring = ['precision_macro', 'recall_macro']
    scores = cross_validate(clf, dataSeq, dataStruc, scoring=scoring,cv=5, 
                            return_train_score=False)
    sorted(scores.keys())
    scores['test_recall_macro']
    DataFrame.from_dict(data=scores, orient='index').to_csv("result/cross_validation_score.csv")
    
    print("Evaluating performance...")
    performance(preds,realStruc)
      
    print("Saving models...")
    filepath = os.path.join('models', 'clf_win.pkl')
    if not os.path.exists('models'):
        os.makedirs('models')
    joblib.dump(clf, filepath)
    
    print("Done!")
    pass

if __name__ == "__main__":   
    for i in range(3,4,2):
        use(i,"data/trainset.dat","data/testset.dat")  
        f = open("result/evaluation.dat",'a')
        f.write("windowsize:"+str(i)+'\n')
        f.close()
