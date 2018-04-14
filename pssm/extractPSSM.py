# run this file in the directory to be proceeded.
import os
import pickle
import numpy as np
import  csv

if __name__ == "__main__":
    path = os.getcwd()
    
    # read all psiblast pssm result files
    for filename in os.listdir(path):
        if filename.endswith(".pssm"):
            with open (filename) as f:
                data = f.read().splitlines()
                for i in range(2,len(data)):
                    data[i] = data[i].split()
        
                # only keep pssm
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
                
        with open(filename+".csv", 'w', newline='') as f:
            wr = csv.writer(f)
            wr.writerows(pssm)
            
print("All PSSMs extracted")


    

