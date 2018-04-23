with open('data/testset.dat') as f:
    with open('data/testset_fasta.dat','w') as g:
        data = f.read().splitlines()
        for i in range(len(data)):
            if i%3 == 1:
               g.write(data[i])
               g.write('\n')
            if i%3 == 2:
               pass
            if i%3 == 0:
               g.write(data[i])
               g.write('\n')  
