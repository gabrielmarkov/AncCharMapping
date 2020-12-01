
# Script corresponding to Fig. 4 in Girard et al., paper in prep for Frontiers in Plant Science

# loading the ape and phytools libraries (R version 3.2.2)
library(ape)
library(phytools)

# setting working directory. Warning: the path has to be adjusted on your own computer!
setwd( "/YourPath/YourWorkingDirectory")

# reading a tree without branch lengths
tree<-read.nexus("cox3_tree_30sp.txt"); 

# reading the table with late SSR pathway presence/absence data
trait <- read.csv("LateSSR_30sp.csv", header =T, row.names = 1);
# coding: Y = Late SSR present; N = Late SSR absent; U = unknown

# transforming the data table into a vector format.
vec_trait<-trait[,1]
names(vec_trait)<-rownames(trait)
head(vec_trait)

# stochastic mapping under the symmetrical model 
simmaptrees <- make.simmap(tree,vec_trait, model="SYM", nsim= 1000, Q="empirical", pi="estimated")

# Below is the expected output:

Using pi estimated from the stationary distribution of Q assuming a flat prior.
pi =
       N        U        Y 
0.333333 0.333333 0.333333 

make.simmap is sampling character histories conditioned on the transition matrix
Q =
            N          U         Y
N -13.0346723  0.4258012  12.60887
U   0.4258012 -0.4258012   0.00000
Y  12.6088711  0.0000000 -12.60887
(estimated using likelihood);
and (mean) root node prior probabilities
pi =
        N         U         Y 
0.3333333 0.3333333 0.3333333 
Done.

# end of expected output


# computing the state frequencies from the stochastic maps for each internal nodes. Posterior probabilities are illustrated by pies. 
simmap_summary<-describe.simmap(simmaptrees, plot=TRUE);
simmap_summary
1000 trees with a mapped discrete character with states:
 N, U, Y 

trees have 97.015 changes between states on average

changes are of the following types:
       N,U    N,Y   U,N U,Y    Y,N Y,U
x->y 1.259 47.253 0.516   0 47.987   0

mean total time spent in each state is:
             N          U         Y   total
raw  3.7545392 0.40146671 3.7720441 7.92805
prop 0.4735766 0.05063877 0.4757846 1.00000


simmap_summary$ace
       N     U     Y
30 0.430 0.125 0.445
31 0.515 0.003 0.482
32 0.498 0.000 0.502
33 0.476 0.000 0.524
34 0.495 0.001 0.504
35 0.464 0.000 0.536
36 0.442 0.000 0.558
37 0.477 0.000 0.523
38 0.374 0.000 0.626
39 0.033 0.000 0.967
40 0.406 0.000 0.594
41 0.534 0.000 0.466
42 0.549 0.000 0.451
43 0.571 0.000 0.429
44 0.574 0.000 0.426
45 0.625 0.000 0.375
46 0.166 0.000 0.834
47 0.161 0.000 0.839
48 0.179 0.000 0.821
49 0.767 0.003 0.230
50 0.804 0.000 0.196
51 0.786 0.000 0.214
52 0.677 0.000 0.323
53 0.460 0.000 0.540
54 0.447 0.000 0.553
55 0.506 0.000 0.494
56 0.701 0.000 0.299
57 0.462 0.000 0.538


# plotting posterior density of stochastic mapping on the tree and saving it into a pdf file.
pdf("Simmap_SYM_LateSSR_30sp.pdf", height=7, width=6)
plot(simmap_summary, cex=0.6, fsize=1);
dev.off()
