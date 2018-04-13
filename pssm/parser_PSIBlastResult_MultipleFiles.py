# read file as an array
#The drawback of these function is that :
# 1. In converting from amino acid from letter format "A,K,M" to binary format (0,0,1)...) It requires the amino acid to present in the right order instead of using dictionary.key. I don't have enough time to create a dictionary and convert it.
# 2. Also, in the pssm trainset, all sequences files must be assigned in the right order, otherwsie it won't correspond to the structure. This can also be improved (although I didn't), we can use the file/sequence name as key and use a dictionary to do this.
import os
import pickle

def PSSMparser(filename):
    import numpy as np
    with open (filename) as f:
        data = f.read().splitlines()
        for i in range(2,len(data)):
            data[i] = data[i].split()
        # only keep the pssm information
        for i in range(3,len(data)):
            del data[i][0:22]
            data[i] = data[i][0:20]
        # only keep line amino acid and amino acid value
        pssm = []
        for i in range(2,len(data)-6):
            pssm.append(data[i])
        pssm[0]=pssm[0][20:40]
        for j in range(len(pssm)):
            pssm[j] = np.asarray(pssm[j])
        pssmAA = []
        for i in range(1,len(pssm)):
            pssmAA.append(pssm[i])
        pssm = []
        for i in pssmAA:
            intPSSM = list(map(int, i))
            pssm.append(intPSSM)
        
        with open(filename+"parsed", 'wb') as fp:
            pickle.dump(pssm, fp)
        pass
        # remember to convert to array

def UniPSSMparser(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".pssm"):
            PSSMparser(filename)
    pass

UniPSSMparser('/psiblast_rawresult')
    
    

