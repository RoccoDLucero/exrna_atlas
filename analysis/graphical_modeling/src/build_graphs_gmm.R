# Build graph using a Gaussian Markov Model
#Oscar, you may need to alter this script so it works with the
# data frames you generated.

library(huge)
library(igraph)

#Supply a data frame of samples x RNA counts with metadata perhaps
#atlas.rpm.lst <- readRDS("./atlas.all.counts.and.biogps.meta")
#rm(atlas.rpm.lst)
#rna.types <- names(atlas.rpm.lst)

#Subset your data and remove metadata fields before graph estimation
z <- atlas.rpm.lst$miRNA
z.cond <- z[z$condition == "Parkinson's Disease",]

#This is to reformat Rocco's data frames. 
rownames(z.cond) <- z.cond$`Sample Name`
z.cond <- z.cond[,-1]
z.cond.meta <- z.cond[,(ncol(z.cond)-16):ncol(z.cond)]
z.cond <- z.cond[,1:(ncol(z.cond)-16)]

#Remove all columns that have no variance
z.cond <- z.cond[,which(round(apply(z.cond, 2, sd), 4) != 0)] 

#Our data must be transformed to standard normal for the modeling
#We use the non-paranormal transformation offered in the HUGE package
z.npn <- huge.npn(z.cond) 

#Perform Graphical LASSO
z.npn.inv.cov <- huge(x = z.npn, method = "glasso", cov.output = F, nlambda = 15)

#Change this according to selected data
#saveRDS(z.npn.inv.cov,"./miRNA_Parkinson_only_precision_mtx")
#z.npn.inv.cov <- readRDS("./miRNA_Parkinson_only_precision_mtx")  
#z.npn.inv.cov <- readRDS("./miRNA_Alzheimer_only_precision_mtx")  

#Select the graphical model with optimal regularization parameter using empirical bayes:
weighted.graph <- huge.select(z.npn.inv.cov,criterion = "ebic")
grph.opt <- weighted.graph$opt.icov
grph.opt@Dimnames <- list(dimnames(z.cond)[[2]],dimnames(z.cond)[[2]])

#Now select subsets of variables to visualize
#Here I find rna present at an average rpm threshold in 
#This will be a function that gets all RNA species with mean read counts at or above a threshold
#Appropriately rename variables
var.means <- by(data = z.cond, INDICES = z.cond.meta$biofluid_name  , FUN = colMeans)
mean.reads <- 15
var.means.thresh <-lapply(var.means,function(x){x[x > mean.reads]})
vmt.names <- lapply(var.means.thresh,names)
vmt.names <- vmt.names[!sapply(vmt.names,is.null)]
vmt.intersect <- Reduce(intersect, vmt.names) #Select miRNA appearing in all biofluids considered
#Subset the graph here:
grph.intsct <-(grph.opt[vmt.intersect,vmt.intersect])
grph.intsct@Dimnames[[1]][grph.intsct@i[which(grph.intsct@x > .1 & grph.intsct@x < 1 )]]

#Turn the adjacency matrix graph structure into a true graph structure and export
#for visualization with CYTOSCAPE.   
grph.intsct.gml <- graph.adjacency(adjmatrix = grph.intsct,mode = "upper",weighted = T,diag =  F)
#write.graph(graph = grph.intsct.gml,file = "./output/Parkinson.graphml", format = "graphml")
write.graph(graph = grph.intsct.gml,file = "./output/Alzheimer.graphml", format = "graphml")


a <- sample(1:ncol(z.cond),1)
hist(z.cond[,a])
hist(z.npn[,a])

