#!/usr/bin/Rscript
############################
library('bio3d')
library('NACEN')
library('igraph')
############################

############################
dsspfile <- "/usr/bin/mkdssp"
MT_betweeness <- c()
MT_closeness <- c()
MT_eigenvector <- c()
Polarity <- c()
Hydrophobicity <- c()

# MT_degree <- c()
# MT_clustering <- c()
# MT_pagerank <- c()
############################

############################
# Wild Type
data <- "/home/wangjingran/APMA/data/pten.pdb"
Net <- suppressMessages(NACENConstructor(PDBFile=data,WeightType = "Polarity",exefile = dsspfile,plotflag=F))
Net_hydro <- suppressMessages(NACENConstructor(PDBFile=data,WeightType = "Hydrophobicity",exefile = dsspfile,plotflag=F))
NetP <- suppressMessages(NACENAnalyzer(Net$AM,Net$NodeWeights))
WT_polarity <- Net$NodeWeights
WT_hydrophobicity <- Net_hydro$NodeWeights
WT_polarity <- as.vector(WT_polarity)
WT_hydrophobicity <- as.vector(WT_hydrophobicity)

net <- NetP$Edgelist
network <- c(net[,1],net[,2])
network <- graph(network)
result <- NetP$NetP

WT_betweeness <- result$B
WT_closeness <- result$C
WT_eigenvector <- suppressWarnings(evcent(network,scale=F)$vector)
# WT_degree <- result$K

# WT_clustering <- suppressWarnings(transitivity(network,type="localundirected"))
# WT_clustering[is.na(clustering)] <- 0
# WT_pagerank <- suppressWarnings(page_rank(network, damping = 0.999)$vector)
############################

############################
for(i in 1:5){
	data <- paste(c("/home/wangjingran/APMA/FoldX/Cancer_",i,".pdb"),collapse="")

    # Use NACEN to calculate nodewights based on Polarity and Hydrophobicity
	Net <- suppressMessages(NACENConstructor(PDBFile=data,WeightType = "Polarity",exefile = dsspfile,plotflag=F))
    Net_hydro <- suppressMessages(NACENConstructor(PDBFile=data,WeightType = "Hydrophobicity",exefile = dsspfile,plotflag=F))

    # Fetch NodeWeights
    polarity <- Net$NodeWeights
    hydrophobicity <- Net_hydro$NodeWeights

    # Change into vectors
    polarity <- as.vector(polarity)
    hydrophobicity <- as.vector(hydrophobicity)

    # Calculate the difference value of MT and WT
    polarity <- polarity - WT_polarity
    hydrophobicity <- hydrophobicity - WT_hydrophobicity

    # Use global polarity and hydrophobicity change
    polarity <- sum(polarity)
    hydrophobicity <- sum(hydrophobicity)

    # Record
    Polarity <- c(Polarity, polarity)
    Hydrophobicity <- c(Hydrophobicity, hydrophobicity)

    # Calculate graph centrality
	NetP <- suppressMessages(NACENAnalyzer(Net$AM,Net$NodeWeights))
	net <- NetP$Edgelist
	result <- NetP$NetP
    
    # Degree
	# degree <- result$K
    # MT_degree <- cbind(MT_degree,degree)

    # Calculate betweenness
	betweeness <- result$B
    MT_betweeness <- cbind(MT_betweeness,betweeness)

    # Calculate closeness
	closeness <- result$C
	MT_closeness <- cbind(MT_closeness,closeness)
	
    # Calculate eigenvetor
	network <- c(net[,1],net[,2])
    network <- graph(network)
    ev <- suppressWarnings(evcent(network,scale=F)$vector)
    MT_eigenvector <- cbind(MT_eigenvector,ev)

    # tr <- suppressWarnings(transitivity(network,type="localundirected"))
    # tr[is.nan(tr)] <- 0
    # pg <- suppressWarnings(page_rank(g, damping = 0.999)$vector)
    # MT_clustering <- cbind(MT_clustering,tr)
    # MT_pagerank <- cbind(MT_pagerank, pg)
}
############################

############################
# Mean graph centrality
MT_betweeness <- rowMeans(MT_betweeness)
MT_closeness <- rowMeans(MT_closeness)
MT_eigenvector <- rowMeans(MT_eigenvector)
# MT_degree <- rowMeans(MT_degree)
# MT_clustering <- rowMeans(MT_clustering)
# MT_pagerank <- rowMeans(MT_pagerank)



Betweeness <- MT_betweeness - WT_betweeness
Closeness <- MT_closeness - WT_closeness
Eigenvector <- MT_eigenvector - WT_eigenvector
# Degree <- MT_degree - WT_degree
# Clustering_coefficient <- MT_clustering - WT_clustering
# PageRank <- MT_pagerank - WT_pagerank

AA_web_1 <- cbind(Betweeness, Closeness)
AA_web_1 <- cbind(AA_web_1, Eigenvector)
AA_web_2 <- cbind(Polarity, Hydrophobicity)

# AA_web <- cbind(AA_web, Degree)
# AA_web <- cbind(AA_web, Clustering_coefficient)
# AA_web <- cbind(AA_web, PageRank)

write.table(AA_web_1,"/home/wangjingran/APMA/data/AAWeb/Cancer.txt",sep="\t",row.names = FALSE)
dNW_file_path <- "/home/wangjingran/APMA/data/AAWeb/dNodeWeight.txt"
if (!file.exists(dNW_file_path)) {
    write.table(AA_web_2, file=dNW_file_path, sep="\t", row.names=FALSE, col.names=TRUE, append=FALSE)
} else {
    write.table(AA_web_2, file=dNW_file_path, sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE)
}
############################
