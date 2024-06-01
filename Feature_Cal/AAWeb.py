import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
bio3d = importr("bio3d")
NACEN = importr("NACEN")
igraph = importr("igraph", robject_translations={'.env':'__env'})
def AAWEB(route,t, category, Mut_PDB,WT_PDB,data_rote):
    '''
    Function to calculate the Amino Acids Web Features for several given mutation type pdb files
    unluckily, the NACEN package is base on R, so you should prepare it before calling this function
    :param route:
    :param t:
    :param category:
    :param Mut_PDB:
    :param WT_PDB:
    :return: several files with Amino Acids Web features
    '''
    element_count = category.count(t)
    r_code = f'''
dsspfile <- "{route}"
MT_degree <- c()
MT_betweeness <- c()
MT_closeness <- c()
MT_eigenvector <- c()
MT_clustering <- c()
for(i in 1:{element_count}){{
	data <- paste(c("{Mut_PDB}/{t}_",i,".pdb"),collapse="")
	Net <- suppressMessages(NACENConstructor(PDBFile=data,WeightType = "Polarity",exefile = dsspfile,plotflag=F))
	NetP <- suppressMessages(NACENAnalyzer(Net$AM,Net$NodeWeights))
	net <- NetP$Edgelist
	result <- NetP$NetP
	degree <- result$K
	betweeness <- result$B
	closeness <- result$C
	MT_degree <- cbind(MT_degree,degree)
	MT_betweeness <- cbind(MT_betweeness,betweeness)
	MT_closeness <- cbind(MT_closeness,closeness)
	
	network <- c(net[,1],net[,2])
    network <- graph(network)
    ev <- suppressWarnings(evcent(network,scale=F)$vector)
    tr <- suppressWarnings(transitivity(network,type="localundirected"))
    tr[is.nan(tr)] <- 0
    MT_eigenvector <- cbind(MT_eigenvector,ev)
    MT_clustering <- cbind(MT_clustering,tr)
}}
MT_degree <- rowMeans(MT_degree)
MT_betweeness <- rowMeans(MT_betweeness)
MT_closeness <- rowMeans(MT_closeness)
MT_eigenvector <- rowMeans(MT_eigenvector)
MT_clustering <- rowMeans(MT_clustering)

data <- "{WT_PDB}"
Net <- suppressMessages(NACENConstructor(PDBFile=data,WeightType = "Polarity",exefile = dsspfile,plotflag=F))
NetP <- suppressMessages(NACENAnalyzer(Net$AM,Net$NodeWeights))
net <- NetP$Edgelist
network <- c(net[,1],net[,2])
network <- graph(network)
result <- NetP$NetP
degree <- result$K
betweeness <- result$B
closeness <- result$C
eigenvector <- suppressWarnings(evcent(network,scale=F)$vector)
clustering <- suppressWarnings(transitivity(network,type="localundirected"))
clustering[is.na(clustering)] <- 0

Betweeness <- MT_betweeness - betweeness
Closeness <- MT_closeness - closeness
Degree <- MT_degree - degree
Eigenvector <- MT_eigenvector - eigenvector
Clustering_coefficient <- MT_clustering - clustering

AA_web <- cbind(Betweeness, Closeness)
AA_web <- cbind(AA_web, Degree)
AA_web <- cbind(AA_web, Eigenvector)
AA_web <- cbind(AA_web, Clustering_coefficient)

write.table(AA_web,"{data_rote}/{t}.txt",sep="\\t",row.names = FALSE)
'''
    robjects.r(str(r_code))

def data_AAW_gener(position, category):
    AAW_data = []
    for i in range(len(position)):
        current_cate = category[i]
        current_data = []
        with open(f"/home/wangjingran/APMA/data/AAWeb/{current_cate}.txt", 'r') as file:
            for line in file:
                columns = line.strip().split('\t')
                current_data.append(columns)
        current_data_features = current_data[position[i]]
        AAW_data.append(current_data_features)
    return AAW_data

