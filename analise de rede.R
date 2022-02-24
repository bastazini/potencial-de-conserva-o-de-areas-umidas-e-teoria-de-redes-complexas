###packages
require(bipartite)
require(igraph)
require(influential)
require(vegan)

rm(list=ls())
###Import data
#Load data
rede=read.table(pipe("pbpaste"), sep="\t", header=T,row.names=1);rede
dim(rede)

##Calculate habitat diversity using equivalent numbers 
habitat=read.table(pipe("pbpaste"), sep="\t", header=T);habitat
eq.numbers=exp(vegan::diversity(habitat, "shannon")); eq.numbers 
##copy and paste to excel, to organize data
clipr::write_clip(eq.numbers)

####calculating nod redundancy
library(networktools)

teste <- goldbricker(t(rede), threshold=0.30)
plot(teste)
teste1 <- net_reduce(data=t(rede), badpairs=teste,method="best_goldbricker")




#### Calcularing IVI
My_graph <- graph_from_incidence_matrix(rede,weighted = TRUE)
Graph_IVI <- ivi(graph = My_graph, mode = "all")
ivi=Graph_IVI[1:48]
clipr::write_clip(ivi.aves)
cent_network.vis(graph = My_graph, cent.metric = Graph_IVI,
                 legend.title = "Valor Integrado de Influencia",
                 plot.title = "", layout= "kk", dist.power=2, legend.position="right", boxed.legend=FALSE, show.labels=FALSE, 
                 node.shape=c("circle", "circle", "circle", "circle", "circle", "circle", 
                              "circle", "circle", "circle", "circle", "circle", "circle", "circle", 
                              "circle", "circle", "circle", "circle", "circle", "circle", "circle", 
                              "circle", "circle", "circle", "circle", "circle", "circle", "circle", 
                              "circle", "circle", "circle", "circle", "circle", "circle", "circle", 
                              "circle", "circle", "circle", "circle", "circle", "circle", "circle", 
                              "circle", "circle", "circle", "circle", "circle", "circle", "circle", 
                              "square", "square", "square", "square", "square", "square", "square", 
                              "square", "square", "square", "square", "square", "square", "square", 
                              "square", "square", "square", "square", "square", "square", "square", 
                              "square", "square", "square", "square", "square", "square", "square", 
                              "square", "square", "square", "square", "square", "square", "square", 
                              "square", "square", "square"))
               
                 ###Order IVI
#negative sign serves as a trick to rank IVI on descending order
ivi.order=rank(-ivi);ivi.order
edit(ivi.order)

###Robustness
##specify "ext.row" based on the IVI rank vector
rob =second.extinct(rede, participant="lower", method="external",ext.row= c(14, 39.5, 13, 5, 33, 27, 41, 10, 9, 1, 22, 37, 25, 
                                                                               7, 26, 38, 21, 39.5, 35.5, 28, 24, 42.5, 42.5, 18, 35.5, 4, 20, 
                                                                               47.5, 45, 29, 44, 47.5, 15, 46, 32, 23, 30, 2, 8, 19, 6, 17, 
                                                                           11, 31, 3, 34, 12, 16))
rob.random=  second.extinct(rede, participant="lower", method="random", nrep=1000,details=FALSE)  
###Robustness
robustness(rob)
robustness(rob.random)
### this fucntion plots the ATC and calculates its exponet
###This is a modification from slope.bipartite taken from Bastazini et al. 2019. Environmental Conservation 46.1 (2019): 52-58.
source("fit.hyperbolic.R")
par(mfrow=c(1,2))
fit.hyperbolica(rob)
fit.hyperbolica(rob.random)
