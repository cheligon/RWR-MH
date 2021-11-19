## Clears R environment and memory
rm(list=ls())
gc()

## Set working directory depending on computer used
node.name <- Sys.info()[[4]]
if (node.name == "Cedre-24"){
  setwd(dir = '/home/chris/CeDRE/RWR-MH/Scripts_and_Files/Network_Generation/')
} else if (node.name == "cedre-14a"){
  setwd(dir = '/Projects_Space/People/Christophe//RWR-MH/Scripts_and_Files/Network_Generation/')
}

################################################################################################################
## 1.- We load the common network generation functions, and we install/load the packages needed. 
################################################################################################################
source("Network_Generation_Utils_CH.R")
source("Network_Generation_Utils.R")
#CRAN.packages <- c("igraph", "reshape")
#bioconductor.packages <- c("biomaRt")
#install.packages.if.necessary(CRAN.packages,bioconductor.packages)
library('igraph')
library('reshape')
library('biomaRt')

# We get the date of execution. Then we will paste it to the name of the network. 
date <- Sys.Date()

# We create a directoy called networks_files to save the files generated.
data.dir <- paste("networks_files_human", date, sep = "_")

if (!file.exists(data.dir)) {
  dir.create(data.dir)
}

################################################################################################################
## 2.- We Fetch RNA-expression data from Protein Atlas. We save Ensemble ID, Gene Name, Sample and FPKM value 
################################################################################################################
## Tissue.  62 tissues on the 16nov2021
tissue.url <- "http://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip"
tissue.expression.table <- get.data.from.zip.at.url(tissue.url, data.dir)
head(tissue.expression.table)
summary(tissue.expression.table)

## Cell type (single cell).  51 cell types on the 16nov2021
celltype.url <- "http://www.proteinatlas.org/download/rna_single_cell_type.tsv.zip"
celltype.expression.table <- get.data.from.zip.at.url(celltype.url, data.dir)
head(celltype.expression.table)
summary(celltype.expression.table)

## Cell Lines.  69 cell lines on the 16nov2021
celline.url <- "http://www.proteinatlas.org/download/rna_celline.tsv.zip"
celline.expression.table <- get.data.from.zip.at.url(celline.url, data.dir)
head(celline.expression.table)
summary(celline.expression.table)
celline.expression.table.2 <- celline.expression.table[-c(4,5)] # removes unnecessary columns
summary(celline.expression.table.2)

names.columns <- c('Gene ID', 'Gene.name', 'Sample', 'Value')
names(tissue.expression.table) <- names.columns
names(celltype.expression.table) <- names.columns
names(celline.expression.table.2) <- names.columns

## We merge records from tissue and cell lines.
data_Coexpr <- rbind(celline.expression.table.2, celltype.expression.table, tissue.expression.table)


################################################################################################################
## 3.- We Check that the HGNC symbols of our Co-expression data are in Biomart. (There are some non-coding RNA) 
## We use a Biomart query.
################################################################################################################
ensemble_HGNC <- check.HGNC.symbols(unique(data_Coexpr$Gene.name))

data3_Coexpr <- data_Coexpr[which(data_Coexpr$Gene.name %in% ensemble_HGNC),]

not.in.Biomart.symbols <- nrow(data_Coexpr) - nrow(data3_Coexpr)

################################################################################################################
## 4.- We transform our data into a matrix in order to calculate the Correlation.
################################################################################################################
## Output a matrix with Gene.name as first column and, in every other column compute the mean Value for each given Sample  
## Then transform "Gene.name"" column into row names
matrix_Coexpr<-cast(data3_Coexpr, Gene.name ~ Sample,mean,value="Value") 
rownames(matrix_Coexpr)<-matrix_Coexpr[,1]
matrix_Coexpr[,1]<-NULL
matrix_Coexpr<-as.matrix(matrix_Coexpr)

### We Compute Spearman correlation (correlation between ranks of values)  
### We need to take the transposed matrix (rows <-> columns)  
### Values are rounded to the 2nd decimal
CR<-round(cor(t(matrix_Coexpr), method="spearman"),2)

################################################################################################################
## 5.- We Select correlation threshold and build a network.
################################################################################################################
### The correlation threshold is set to 0.7. All pairs of genes that have a correlation above this threshold 
### (positive or negative) are considered as linked in the Co-Expression network.  

#Gets geneIDs of genes with score > 0
inds2 <- which(abs(CR[1:nrow(CR),]) >0, arr.ind=TRUE)
rnames2 <- rownames(CR)[inds2[,1]]
cnames2 <- colnames(CR)[inds2[,2]]

#Get values of correlations > 0
values2 <- CR[which(abs(CR) > 0)]

#Free some memory
rm(inds2)
rm(CR)
gc()

#Generates table of co-expression [Gene1, Gene2]
#as data frame to allow different data types to be stored
CT.df <- as.data.frame(cbind(rnames2, cnames2))

#Free some memory
rm(CT)
rm(cnames2)
rm(rnames2)
gc()

#Add weight(=absolute spearman correlation) to the data frame
CT.df2 <- cbind(CT.df, as.numeric(abs(values2)))

#free some memory
rm(CT.df)
gc()

#Rename columns
colnames(CT.df2) <- c('Gene1', 'Gene2', 'weight')

#Free some memory
rm(values2)
gc()

CT.dummy <- read.csv('CT.dummy')
CT.dummy <- CT.dummy[-1]
head(CT.dummy)
colnames(CT.dummy) <- c('Gene1', 'Gene2', 'weight')
head(CT.dummy)
g <- graph_from_data_frame(CT.dummy)
is_weighted(g)
g <- igraph::simplify(g,remove.multiple = TRUE, remove.loops = TRUE)
isolates <- which(degree(g, mode = c("all")) == 0) - 1
g <- delete.vertices(g, names(isolates))

print(paste("Number of Edges: ", ecount(g),sep="",collapse = NULL))
print(paste("Number of Nodes: ", vcount(g),sep="", collapse = NULL))
graph_name <- paste("networks_files_CH/Co-Expression_",date,".gr",sep="",collapse = "")

write.graph(g, graph_name, format=c("ncol"))
