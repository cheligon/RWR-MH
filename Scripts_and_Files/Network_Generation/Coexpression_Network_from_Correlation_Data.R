################################################################################################################
## 0.- Prepare environment
################################################################################################################

## Delete objects and variables & free memory
rm(list=ls())
gc()

## Load libraries
library(igraph)

## Set work dir depending on computer (= node) used
node.name <- Sys.info()[[4]]
if (node.name == "Cedre-24"){
  setwd(dir = '/home/chris/CeDRE/RWR-MH/Scripts_and_Files/Network_Generation/')
} else if (node.name == "cedre-14a"){
  setwd(dir = '/home/heligon/RWR-MH/Scripts_and_Files/Network_Generation/')
}

## We get the date of execution. Then we will paste it to the name of the network. 
date <- Sys.Date()

## We create a directory called networks_files to save the files generated.
data.dir <- paste("networks_files_human", date, sep = "_")

if (!file.exists(data.dir)) {
  dir.create(data.dir)
}


################################################################################################################
## 1.- Load Correlation Data
################################################################################################################
matrix_file <- './networks_files_human_2021-11-22/Correlation_matrix_2021-11-22.RData'

load(file = matrix_file) # CR

################################################################################################################
## 2.- We transform the matrix into a network, remove self-loops and isolated components.
################################################################################################################

## Build network from matrix
net <- igraph::graph.adjacency(CR, mode = 'undirected', weighted = TRUE)

## Simplify network
net <- igraph::simplify(net, remove.multiple = TRUE, remove.loops = TRUE)

## Remove isolated vertices
isolates <- which(degree(net, mode = c("all")) == 0) - 1
net <- delete.vertices(net, names(isolates))
  
print(paste("Number of Edges: ", ecount(net),sep="",collapse = NULL))
print(paste("Number of Nodes: ", vcount(net),sep="", collapse = NULL))
graph_name <- paste("networks_files/Co-Expression_",date,".gr",sep="",collapse = "")

write.graph(net, graph_name, format=c("ncol"),  weights=NULL)
