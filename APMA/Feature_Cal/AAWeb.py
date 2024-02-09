import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
bio3d = importr("bio3d")
NACEN = importr("NACEN")
igragh = importr("igraph")
from .. import category
from .. import Mut_PDB
from .. import WT_PDB
set_category = list(set(category))
for cate in set_category:
    element_count = category.count(cate)
    r_code = f'''
    dsspfile <- "../data/dssp-3.0.0.exe"
    MT_degree <- c()
    MT_betweeness <- c()
    MT_closeness <- c()
    MT_eigenvector <- c()
    MT_clustering <- c()
    for(i in 1:{element_count}){{
	    data <- paste(c("{Mut_PDB}/{cate}_",i,".pdb"),collapse="")
	    Net <- NACENConstructor(PDBFile=data,WeightType = "Polarity",exefile = dsspfile,plotflag=T)
	    NetP <- NACENAnalyzer(Net$AM,Net$NodeWeights)
	    result <- NetP$NetP
	    degree <- result$K
	    betweeness <- result$B
	    closeness <- result$C
	    MT_degree <- cbind(MT_degree,degree)
	    MT_betweeness <- cbind(MT_betweeness,betweeness)
	    MT_closeness <- cbind(MT_closeness,closeness)
        net <- NetP$Edgelist
	    network <- c(net[,1],net[,2])
	    network <- graph(network)
	    eigenvector <- evcent(network,scale=F)$vector
	    clustering <- transitivity(network,type="localundirected")
	    MT_eigenvector <- cbind(MT_eigenvector,eigenvector)
	    MT_clustering <- cbind(MT_clustering,clustering)
    }}
    MT_degree <- rowMeans(MT_degree)
    MT_betweeness <- rowMeans(MT_betweeness)
    MT_closeness <- rowMeans(MT_closeness)
    MT_eigenvector <- rowMeans(MT_eigenvetor) 
    MT_clustering <- rowMeans(MT_clustering)
    # 计算野生型
    data <- paste(c("{Mut_PDB}",".pdb"),collapse="")
	Net <- NACENConstructor(PDBFile=data,WeightType = "Polarity",exefile = dsspfile,plotflag=T)
	NetP <- NACENAnalyzer(Net$AM,Net$NodeWeights)
	result <- NetP$NetP
	degree <- result$K
	betweeness <- result$B
	closeness <- result$C
	net <- NetP$Edgelist
	network <- c(net[,1],net[,2])
	network <- graph(network)
	eigenvector <- evcent(network,scale=F)$vector
	clustering <- transitivity(network,type="localundirected")
	# 整合指标
	Betweeness <- MT_betweeness - betweeness
	Closeness <- MT_closeness - closeness
	Degree <- MT_degree - degree
	Eigenvector <- MT_eigenvector - eigenvector
    Clustering <- MT_clustering - clustering
    
    AA_web <- cbind(Betweeness, Closeness)
    AA_web <- cbind(AA_web, Degree)
    AA_web <- cbind(AA_web, Eigenvector)
    AA_web <- cbind(AA_web, Clustering)
    write.table（AA_web,"../data/AAWeb/{cate}.txt",sep="\t")
    '''


robjects.r(r_code)


