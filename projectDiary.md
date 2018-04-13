# 20 Feb

1. I started to __write a parser__ to analyse the data feature
2. Project __file management__: creating folders for notes, which is the most important thing for me
3. Learning __git__: I finally figured out how to use git pull. Also I learned git rm, git clone,ect.(See in _gitNotes_)

# 21 Feb

1.I deleted the__ .gitignore__ file and lost everything.

```git
git rm -r --cached .
git add .
git commit -m "Removing all files in .gitignore"
```
A SAD BRAND NEW START...

2.I couldn't __assign two values from a list to one key__ in a dictionary.
I tried to use dict(zip()) and also dict.append() The error message was "lists are unhashable"

And I am now wondering __why I need to use a dictionary__. There is not an obvious benefit of using dictionary compared with using lists.
Also when should we __turn different amino acids into numberic vectors__? 

# 22 Feb 

1.Preparing for the __convolution neural network__ journal club

I understood the mathematical meaning of convolution,but still couldn't understand the CNN
I also found a 3d visualization website of CNN in imagine recognization.It was so cool!

[3D CNN visualization](http://scs.ryerson.ca/~aharley/vis/conv/)

2.Writing the project plan

3.I found 5 papers of SVM in protein structure prediction and read two of them.

# 25 Feb
1.Learning svm tutorial

# 26 Feb
Q1: what x, y represent (3rd block in ski-learn.py)? 
Q2: why use float instead of integer in clf.predict (3rd block in ski-learn.py)?
Q3: why my clf.support returns a reversed list compared to the example code on the website?
Q4: what is muliclass classification
Q5: differences between 'one vs one' and 'one vs the rest'

A5: The difference is the number of classifiers you have to learn, which strongly correlates with the decision boundary they create.

Assume you have NN different classes. One vs all will train one classifier per class in total NN classifiers. For class ii it will assume ii-labels as positive and the rest as negative. This often leads to imbalanced datasets meaning generic SVM might not work, but still there are some workarounds.

In one vs one you have to train a separate classifier for each different pair of labels. This leads to N(N−1)2N(N−1)2 classifiers. This is much less sensitive to the problems of imbalanced datasets but is much more computationally expensive.

# 27 Feb
1.Prepare for the presentation (see /KB8024/notes/journal/)
2.An improved parser
3.Preprocessing
4.OneHotEncoder

# 1 March
Finish the most primary version of predictor

# 2 March
1. Journal Club
2. Presentation

# 3 March
1. write a peer presentation evaluation

# 6 March
1. Change data input to dataframe format
2. succeedd in reading multiple sequence information
3. succeedd in onehotencode and labelencode (in a brutal way)
4. added to the head and tail of my data according to the window size

# 7 March
1. read paper: Chen, Ke, Lukasz Kurgan, and Jishou Ruan. "Optimization of the sliding window size for protein structure prediction." Computational Intelligence and Bioinformatics and Computational Biology, 2006. CIBCB'06. 2006 IEEE Symposium on. IEEE, 2006.
2. slide along the window
 
# 8 March
1. attending the presentation
2. Cross validation
3. tried array, matrix.flatten, and many other methods to adjust the shape of the input data.
 
# 9 March
1. solved the input data shape problems
2. built crossvalidation, predictor, testdataparser, trainsetpreprocess,windowsize files
3. JC club and group meeting

# 11 March
1. solved the windowsize problems
2. Tested windowsize optimization, windowsize = (3,21,step = 2)
3. finished crossvalidation, including K-fold, and LOO method. Read papers and found the best split could be 7-fold.

# 12 March
1. tried different svm functions， SVC，NuSVC and Linear SVC

# 13 March
1. finished psi-blast
2. read about PSSM