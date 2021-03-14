
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
# character encoding: Y = Late SSR present; N = Late SSR absent; U = unknown

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

trees have 614.77 changes between states on average

changes are of the following types:
     N,U    N,Y U,N   U,Y     Y,N   Y,U
x->y   0 305.83   0 0.511 307.164 1.265

mean total time spent in each state is:
             N          U         Y   total
raw  3.7818801 0.40812895 3.7380410 7.92805
prop 0.4770253 0.05147911 0.4714956 1.00000


simmap_summary$ace
       N     U     Y
30 0.457 0.119 0.424
31 0.531 0.003 0.466
32 0.520 0.000 0.480
33 0.472 0.000 0.528
34 0.484 0.000 0.516
35 0.518 0.000 0.482
36 0.521 0.000 0.479
37 0.514 0.000 0.486
38 0.494 0.000 0.506
39 0.616 0.000 0.384
40 0.493 0.000 0.507
41 0.506 0.000 0.494
42 0.489 0.000 0.511
43 0.515 0.000 0.485
44 0.506 0.000 0.494
45 0.499 0.000 0.501
46 0.509 0.000 0.491
47 0.515 0.000 0.485
48 0.513 0.000 0.487
49 0.513 0.004 0.483
50 0.476 0.000 0.524
51 0.485 0.000 0.515
52 0.514 0.000 0.486
53 0.535 0.001 0.464
54 0.519 0.000 0.481
55 0.500 0.000 0.500
56 0.481 0.000 0.519
57 0.514 0.000 0.486

# numbers in first column refer to node labels, to plot them type:
plot(tree)
nodelabels()
# node 37 refers to the last common ancestor of Laminariales and Ectocarpales, and node 49 is the last 
# common ancestor of Chordariaceae. Probability states at both nodes are discussed in the manuscript.

# plotting posterior density of stochastic mapping on the tree and saving it into a pdf file.
pdf("Simmap_SYM_LateSSR_30sp.pdf", height=7, width=6)
plot(simmap_summary, cex=0.6, fsize=1);
dev.off()
