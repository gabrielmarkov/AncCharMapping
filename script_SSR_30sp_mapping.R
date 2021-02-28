
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
N -81.65273  0.0000000  81.6527308
U   0.00000 -0.4347289   0.4347289
Y  81.65273  0.4347289 -82.0874597
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

trees have 5.297 changes between states on average

changes are of the following types:
       N,U    N,Y   U,N U,Y    Y,N Y,U
x->y   0 1.286   0 0.286 2.61 1.115

mean total time spent in each state is:
             N          U         Y   total
raw  0.9565393 0.30185465 6.6696560 7.92805
prop 0.1206525 0.03807426 0.8412732 1.00000


simmap_summary$ace
       N     U     Y
30 0.332 0.085 0.583
31 0.001 0.001 0.998
32 0.001 0.000 0.999
33 0.000 0.000 1.000
34 0.001 0.000 0.999
35 0.001 0.000 0.999
36 0.000 0.000 1.000
37 0.000 0.000 1.000
38 0.000 0.000 1.000
39 0.000 0.000 1.000
40 0.000 0.000 1.000
41 0.000 0.000 1.000
42 0.000 0.000 1.000
43 0.000 0.000 1.000
44 0.001 0.000 0.999
45 0.000 0.000 1.000
46 0.000 0.000 1.000
47 0.000 0.000 1.000
48 0.000 0.000 1.000
49 0.000 0.000 1.000
50 0.000 0.000 1.000
51 0.000 0.000 1.000
52 0.001 0.000 0.999
53 0.000 0.000 1.000
54 0.000 0.000 1.000
55 0.000 0.000 1.000
56 0.006 0.000 0.994
57 0.000 0.000 1.000

# numbers in first column refer to node labels, to plot them type:
plot(tree)
nodelabels()
# node 37 refers to the last common ancestor of Laminariales and Ectocarpales, discussed in the manuscript

# plotting posterior density of stochastic mapping on the tree and saving it into a pdf file.
pdf("Simmap_SYM_LateSSR_30sp.pdf", height=7, width=6)
plot(simmap_summary, cex=0.6, fsize=1);
dev.off()
