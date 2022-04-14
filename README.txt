These are the two programs needed to apply the method described by Montgomery Slatkin
in the 2022 paper, "Joint estimation of selection intensity and mutation rate under
balancing selection with application to HLA." Genetics XX:XX-XX.

The program MWspectrum.c generates the normalized frequency spectrum, denoted by phi_i
in the paper, using the method developed by Muirhead and Wakeley in their 2009
Genetics paper. MWspectrum takes no input parameters. It generates a large output file
containing phi_i (called a in the program). The file is named MWspectrum.out and is
approximately 73.5 mb. This file has to be created before the data analysis program
(called dataanalysis.c) can be run. 

dataanalysis.c carries out the analysis described in Genetics paper. The input file 
containing the data to be analyzed is on the command line. The output is in a file
with the same name with .output appended. 

The input file contains allele counts for each locus on each line. The counts can be 
separated by spaces or tabs but not carriage returns. There may be one or more zero
counts. The data for each locus is separated by a carriage return. The program will
read each line until the end of file is reached. The program will also write the results
for each sample to the screen.

For each line the program produces 9 numbers:
1: k, the number of alleles
2: the number of singletons (the number of alleles found in only one copy)
3: the number of doubletons (the number of alleles found in two copies)
4: the number of tripletons (the number of alleles found in three copies)
5: F, the sum of squares of allele frequencies
6: P, the probability from the Ewens-Watterson test. P is the fraction of 1000
   values of F generated under neutrality, that are less than the observed value.
7: Estimated value of theta (theta hat in the paper).
8: Estimated value of S (S hat in the paper).
9: Product of k and F (i. e. the product of the values in column 1 and column 5).

The file 1000Gdata contains the allele counts obtained from the supplement of
the Gourraud et al. (2014) paper. The values generated are those in Table 1 of
Slatkin (2022). 

Comments

All the calculations are deterministic except for the results of the Ewens-Watterson
test. The random number seed is initialized from the time() function in C. Therefore
the P values for the same data set will differ slightly when the program is run more
than once. 

The estimates of S and theta are obtained from a grid search as described in the paper.
The file MWspectrum.out calculates the frequency spectrum  for 0≤S≤1000 in steps of 20
and 1≤theta≤40 in steps of 1. I found this range of values suitable for analyzing the 1000
Genomes data. The programs could be adjusted to allow for different ranges and step
sizes but you would probably need my help in doing that. Please contact me (slatkin@berkeley.edu)
if you want to make such changes.
