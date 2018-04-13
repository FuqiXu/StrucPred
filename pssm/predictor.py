'''
Function : predict with cross validation scores
ToUseCV-Input : 3 line data file, windowsize,svmMethod, cvMethod.(The default windowsize is 13, svmMethod = svmSVC, cvMethod = .....manually)
ToUseCV-Output: print cross validation scores
'''
# import data parser
from parser_RawData_Numberize import TopoArray
from window_size_SeqStructure import WindowSeq
from parser_RawData_Numberize import Parser
from parser_OriginalDataset_withWindow import TestSeqFile
from parser_OriginalDataset_withWindow import TestTopoFile
#import cross validation functions
from crossValidation_Kfold_LOO import KfoldCV
from crossValidation_Kfold_LOO import LooCV
#import other packages
from sklearn import svm
import numpy as np
import os 
from parser_PSSMtoSVM_multipleFiles import UsePSSM
from parser_PredictionResult import Topo_Converter_Rev 
# import svm functions
from predictionFunction_SVM_RF_SDT import svmSVC
from predictionFunction_SVM_RF_SDT import linSVC
from predictionFunction_SVM_RF_SDT import randomforest
from predictionFunction_SVM_RF_SDT import SimpleDeci
from predictionEval import score

path = os.path.abspath('.')
os.chdir(path)

def PredictorSVM(trainfilename,testfilename,windowsize, svmfunction = linSVC,cvfunc=KfoldCV):
    os.chdir(path)    
    #parse data
    trainData = Parser(trainfilename)       
    topoData = TopoArray(trainData) 
    seqData = WindowSeq(windowsize,trainData)
    testSeq = TestSeqFile(windowsize,testfilename)
    testTopoReal = TestTopoFile(testfilename)
    
    #cross validation
    #please check the cross validation function. I didn't succeed to use cross validation in my predictor.
    #lasso = linear_model.Lasso()
    #topoData_pred = cross_val_predict(lasso, seqData, topoData)
    #print(topoData_pred)
    
    seqData2 = UsePSSM(windowsize)
    os.chdir(path)
    seqPSSM = []
    # the orignal pssm sequence data set is 3D, (#samples,#AA,20*windowsize)
    # to use pssm data in svm, we need to convert it to (#AAtotal,20*windowsize)
    for i in seqData2:
        seqPSSM.append(i)
    print('builiding model')
    model = svmfunction(seqData,topoData,testSeq)
    # I don't know why my PSSM prediction predicts worse than seqData, so I used the sequence without PSSM
    print("predicting")
    a = model.predict(testSeq)
    Topo_Converter_Rev(a)
    #score prediction
    score(testTopoReal,a)
    pass
        #if testTopoReal[i]==a[i] and testTopoReal[i]=

#split the last 30 sequences as test set
PredictorSVM('..\\dataset\\trainset.fasta','..\\dataset\\testset.fasta',windowsize=8)
 

