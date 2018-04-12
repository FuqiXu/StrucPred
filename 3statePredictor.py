# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 20:32:51 2018

@author: Fuqi Xu
"""
from pandas.core.frame import DataFrame
import numpy as np
import os

path = os.getcwd()

###### Parsing train dataset######
print("parsing data...")

# Read data from 3 line fasta file and store them in a data frame
def rawtoframe(filename):
    seqID1, seq1, seqTopo1= [], [], []
    
    # Read the three line data into three lists
    with open(filename) as f:
        data = f.read().splitlines()
        for i in range(len(data)):
            if i%3 == 1:
               seq1.append(data[i])
            if i%3 == 2:
               seqTopo1.append(data[i])
            if i%3 == 0:
               seqID1.append(data[i])
              
    # Store all data into a dataframe
    seqData1 = {
            "seqID":seqID1,
            "seq":seq1,
            "seqTopo":seqTopo1
            }
    seqData = DataFrame(seqData1)

    # Convert every sequence and sequence topology from list to arrays, so that the sequence can be processed residue by residue.
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
    
    # Converting residues into numbers
    seq_converter(data.seq)
    topo_converter(data.seqTopo)
    
    return data


########### adding winodws #############
print("adding window...")

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
        
    # Creating a slide window.The basic element in one window is #windowsize*Amino acids. 
    # Scanning every protein   
    for m in range(len(data)):
        seq_single = []
        #Scanning every residue
        for p in range(len(data.seqTopo[m])):
            temp = []
            for n in range(windowsize):
                temp.append(data.seq[m][p+n])
            seq_single.append(temp)    
        data.seq[m]=seq_single
        
    return data

# transfering data into a binary array to be used in svm
def data_svm(data):
    seq = []
    struc = []
    for i in range(len(data)):
        for j in range(len(data.seq[i])):
            seq.append(data.seq[i][j])
        for k in range(len(data.seqTopo[i])):
            struc.append(data.seqTopo[i][k])
    dataSVM = DataFrame({
            'seq':seq,
            'struc':struc
                           })
    return dataSVM
### Cross validation ### 

'''
### cross validation ###
def LooCV(x,y):
    loo = LeaveOneOut()
    for train, test in loo.split(x, y):
        return train, test

### svm model ###
# SVC prediction，using one-versus-one classifier, n_class*(n-1)/2 classifiers are built
def svmSVC(seq,topo,testSeq):
    clf = svm.SVC()
    clf.fit(seq,topo)
	 
    a = clf.predict(testSeq)
    prediction = open('prediction.txt','a')
    PredTopo = Topo_Converter_Rev(a)
    prediction.write(PredTopo)
    prediction.write("\n")
    prediction.close()
    #scores = cross_val_score(clf, trainSeq, trainTopo, cv=5, verbose=40, n_jobs=-1)
    #print('scores')
    #print(scores)
    return clf

#Linear SVC, using “one-versus-rest" classifiers
def linSVC(trainSeq,trainTopo,testSeq):
    lin_clf = svm.LinearSVC(C=1.0).fit(trainSeq,trainTopo)
    prediction = open('prediction.txt','a')
    PredTopo = Topo_Converter_Rev(lin_clf.predict(testSeq))
    prediction.write(PredTopo)
    prediction.write("\n")
    prediction.close()
    #scores = cross_val_score(lin_clf, trainSeq, trainTopo, cv=5, verbose=40, n_jobs=-1)
    #print('scores')
    #print(scores)
    return lin_clf

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


### evaluation ###
### pssm ###
### output formating###
### Test file parser（3 line fasta format) ###
def testseq(windowsize,filename):
    
### novel sequence parser ###    
'''
### main ###
dataBinary = binary_rawdata("data/test1.txt")
dataWind = data_window(5,dataBinary)
dataSVM = data_svm(dataWind)
'''
from sklearn.preprocessing import MultiLabelBinarizer

from sklearn.model_selection import cross_validate
from sklearn.metrics import recall_score
from sklearn import svm
scoring = ['precision_macro', 'recall_macro']
clf = svm.SVC(kernel='linear', C=1, random_state=0)
scores = cross_validate(clf, data.seq, data.seqTopo, scoring=scoring,cv=5, 
                        return_train_score=False)
sorted(scores.keys())
scores['test_recall_macro']
'''