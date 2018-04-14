# The secondary structure data is downloaded at https://cdn.rcsb.org/etl/kabschSander/ss.txt.gz
# To convert the sturcture into 3 state, we keeped the H-helix definition and the E -sheet definition. And converted all the rest structure into C.
def sequenceConverter(filename):
    swiss = open("/BioinfoProject/KB8025/dataset/swissprot.txt",'w')
    with open(filename) as f:
        data = f.read().splitlines()
        for i in range(0,len(data),4):
            swiss.write(data[i]+ '\n')
            swiss.write(data[i+1]+ '\n')
            a = list((data[i+3]))
            for j in range(len(a)):
                if a[j] =='H' or a[j]=='E':
                    pass
                else: 
                    a[j] = 'C'
            b= ''.join(a)
            swiss.write(b)
            swiss.write('\n')
sequenceConverter('/BioinfoProject/KB8025/dataset/ss.txt')