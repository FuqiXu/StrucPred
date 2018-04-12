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
    for i in range(len(rawtoframe("data/test1.txt").seq)):
        testSeqSingle=''.join(rawtoframe("data/test1.txt").seq[i])
        testSeq.append(testSeqSingle) 
     
    testBinary = binary_rawdata(filename)
    testWind = data_window(windowsize,testBinary)
    testData = []
    for i in range(len(testWind)):
        seqData = testWind.iloc[i]
        testData.append(seqData)
        
    return testSeq,testData


### Save prediction result ###
def sav_pred(prediction,testdata,testSeq,testfilename):
    # Create a prediction result folder a .dat file"
    filepath = os.path.join('result', 'pred.dat')
    if not os.path.exists('result'):
        os.makedirs('result')
    f = open(filepath, "w")
    
    # Change numberic prediction to letters
    predStruc = []
    preds = []
    for i in range(len(prediction)):
        for j in range(len(prediction[i])):
            if prediction[i][j]==0:
                predStruc.append ('H')
            if prediction[i][j]==1:
                predStruc.append ('E')
            if prediction[i][j]==2:
                predStruc.append ('C')
        pred = str.join("",predStruc)
        preds.append(pred)
        
    # Save prediction result to file
    for m in range(len(testdata)):
        f.write(testdata[m].seqID)
        f.write("\n")
        f.write(testSeq[m])
        f.write("\n")
        f.write(preds[m])
        f.write("\n")
    f.close()

'''

### svm model ###

    clf = svm.SVC()
    linSVC(trainSeq,trainTopo,testSeq)


    #scores = cross_val_score(clf, trainSeq, trainTopo, cv=5, verbose=40, n_jobs=-1)
    #print('scores')
'''
'''
#using random forest classifier)
def randomforest(trainSeq,trainTopo,testseq):
    X, y = make_classification(trainSeq, trainTopo)
    rf_clf = RandomForestClassifier(max_depth=2, random_state=0)
    rf_clf.fit(X, y)
    print("prediction_random Forest:")
    print(rf_clf.predict(testseq))
    prediction = open('prediction.txt','a')
    PredTopo = Topo_Converter_Rev(rf_clf.predict(testseq))
    prediction.write(PredTopo)
    prediction.close()
    scores = cross_val_score(rf_clf, trainSeq, trainTopo, cv=5, verbose=40, n_jobs=-1) 
    print('scores')
    print(scores)
    return rf_clf    
'''

### evaluation ###
### pssm ###
### output formating###

### novel sequence parser ###  


if __name__ == "__main__":
    print("Parsing data...")
    dataBinary = binary_rawdata("data/test2.dat")
    
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
    testSeq,testData = test_fasta("data/test1.txt",3)
    
    print("Predicting...")
    preds = []
    for i in range(len(testData)):
        pred = clf.predict(testData[i].seq)
        preds.append(pred)
    
    print("Saving prediction...")
    predResult = sav_pred(preds,testData,testSeq,'test1')
    
    print("Cross validating...")
    scoring = ['precision_macro', 'recall_macro']
    scores = cross_validate(clf, dataSeq, dataStruc, scoring=scoring,cv=5, 
                            return_train_score=False)
    sorted(scores.keys())
    scores['test_recall_macro']

    DataFrame.from_dict(data=scores, orient='index').to_csv("result/cross_validation_score.csv")
    #for key, value in sorted(scores.items()):
     #   score.write(str(key)+'\t'+str(value)+'\n')
    #score.close()
  
    print("Done!")

   
