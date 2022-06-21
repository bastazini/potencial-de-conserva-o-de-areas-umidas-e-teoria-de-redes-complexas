###packages
require(bipartite)
require(igraph)
require(influential)
require(vegan)
require(networktools)

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

####Node redundancy
teste <- goldbricker(t(rede), threshold = 0.25, corMin = 0.5)
plot(teste)
summary(teste)
teste1 <- net_reduce(data=t(rede), badpairs=teste,method="best_goldbricker")

####Assortativity_degree
assortativity_degree
assortativity_degree(My_graph,directed = FALSE)


#### Calcularing IVI
My_graph <- graph_from_incidence_matrix(rede,weighted = TRUE)
Graph_IVI <- ivi(graph = My_graph, mode = "all")
ivi=Graph_IVI[1:46]
clipr::write_clip(ivi)

#### Plotting (change shape parameter for each network)
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
rob =second.extinct(rede, participant="lower", method="external",ext.row= c(17, 18, 10, 30, 12, 2, 8, 11, 7, 19, 3, 4, 16, 22, 
                                                                            6, 1, 28.5, 28.5, 23.5, 15, 26.5, 26.5, 5, 13, 14, 9, 25, 23.5, 
                                                                            20, 21))
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
