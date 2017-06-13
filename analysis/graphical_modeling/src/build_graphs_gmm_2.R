# Build graph using a Gaussian Markov Model

source("./exrna_atlas_functions.R")
library(huge)
library(igraph)
library(tsne)

atlas.rpm.lst <- readRDS("./atlas.all.counts.and.biogps.meta")
#rm(atlas.rpm.lst)
rna.types <- names(atlas.rpm.lst)

#For when we want to look at covariance of differnt RNA species
z <- atlas.rpm.lst$miRNA
z.healthy <- z[z$condition == "Healthy Control",]
rownames(z.healthy) <- z.healthy$`Sample Name`
z.healthy <- z.healthy[,-1]
z.healthy.meta <- z.healthy[,(ncol(z.healthy)-16):ncol(z.healthy)]
z.healthy <- z.healthy[,1:(ncol(z.healthy)-16)]
z.healthy <- z.healthy[,which(round(apply(z.healthy, 2, sd), 4) != 0)] #Remove all columns that have no variance

 levels(z.healthy.meta$biofluid_name)
#Perform Graphical LAsso
healthy.by.biofluid.graphs <- list(CSF = NULL, PLASMA=NULL, SERUM = NULL, SALIVA= NULL)
i <- 0
for(bf in levels(z.healthy.meta$biofluid_name)[c(2,4,5,6)]){
    i <- i+1
    z.healthy.bf <- z.healthy[which(z.healthy.meta$biofluid_name == bf),]
    z.healthy.bf <- z.healthy.bf[,which(round(apply(z.healthy.bf, 2, sd), 4) != 0)]
    z.healthy.bf.npn <- huge.npn(z.healthy.bf)
    z.healthy.bf.npn.inv.cov <- huge(x = z.healthy.bf.npn,
                              method = "glasso",
                              cov.output = F,
                              nlambda = 15)
    
    healthy.by.biofluid.graphs[[i]] <- z.healthy.bf.npn.inv.cov

}
saveRDS(healthy.by.biofluid.graphs,"./miRNA_Healthy_by_biofluid_precision_mtx")
#saveRDS(z.healthy.npn.inv.cov,"./miRNA_Healthy_only_precision_mtx")
#Select the graphical model with optimal regularization parameter using empirical bayes:
all.healthy.weighted.graph <- huge.select(z.healthy.npn.inv.cov,criterion = "ebic")

###################################################################################################
#This will be a function that gets all RNA species with mean read counts at or above a threshold
#Appropriately rename variables
zz <- by(data = z.healthy,INDICES = z.healthy.meta$biofluid_name  , FUN = colMeans)
mean.reads <- 15
zzz <-lapply(zz,function(x){x[x > mean.reads]})
zzz.names <- lapply(zzz,names)
zzz.names <- zzz.names[!sapply(zzz.names,is.null)]
zzz.intersect <- Reduce(intersect, zzz.names) #Select miRNA appearing in all biofluids considered
zzz.intersect
###################################################################################################
#Now select graph subset based on the mutual intersect of miRNA at threshold in all 4 biofluids
grph.opt <- all.healthy.weighted.graph$opt.icov
grph.opt@Dimnames <- list(dimnames(z.healthy)[[2]],dimnames(z.healthy)[[2]])
grph.intsct <-(grph.opt[zzz.intersect,zzz.intersect])
hist((grph.intsct@x),100)


grph.intsct@Dimnames[[1]][grph.intsct@i[which(grph.intsct@x > .1 & grph.intsct@x < 1 )]]
grph.intsct.gml <- graph.adjacency(adjmatrix = grph.intsct,mode = "upper",weighted = T,diag =  F)


#grph <- grph.intsct.gml
grph <- list(est.igCSF,est.igPLS,est.igSER,est.igSAL)[[1]]

uniq.edges.csf <- simplify(Reduce(f = difference,x = list((est.igCSF),(est.igPLS),(est.igSER),(est.igSAL))))
sort(abs(E(uniq.edges.csf)$weight),decreasing = T)
uniq.edges.pls <- simplify(Reduce(f = difference,x = list((est.igPLS),(est.igSER),(est.igSAL),(est.igCSF))))

common.nodes.all.bf <- Reduce(f = intersection,x = list(V(est.igCSF),V(est.igPLS),V(est.igSER),V(est.igSAL)))

num.edges <- length(E(grph)$weight)


grph <- uniq.edges.csf
grph <- simplify(grph)
thr <- quantile( abs(E(grph)$weight),probs = seq(0,1,.001))
T.M50 <-  matrix(nrow = length(thr),ncol = 4,
                 dimnames = list(NULL,c("threshold_weight","M50a","M50b","Number_of_clusters")) ) 
for(cur.thr in thr){
    
    #Definition of  M50: 50% of nodes are in modules of this size or greater  
    grph.cpy <- delete.edges(grph,which(abs(E(grph)$weight) < cur.thr))
    num.c <- no.clusters(grph.cpy)
    c <- sort(components(grph.cpy)$csize,decreasing = T)
    
    l <- c()
    for(mod in c){l <- c(l,rep(mod,mod))}
    m50a <- median(l)
    m <- sum(c)/2
    m50b<- which(cumsum(c) > m)[1]
    #print(c)[1:15]
    idx <- which(thr == cur.thr)
    T.M50[idx,] <- c(cur.thr,m50a,m50b,num.c)
 
}

par(mar = c(5,5,5,7))
plot(x = T.M50[,1],y = T.M50[,2],xlab = "Minimum abslolute edge weight for inclusion",
     ylab =  NA, ylim = c(0,1000),
     main = "M50 vs. edge-weight threshold\nfor clusters based on 101 co-expressed miRNA\n(CSF,Saliva,Serum,Plasma)",
     sub = "(biofluid mean RPM >= 15 for all miRNA)",
     col = c("purple2"), cex = .8, pch = 20)
points(x = T.M50[,1], y = T.M50[,4], col = "black", cex = .5, pch = 16 )
axis(side = 4)
mtext(side = 2, line = 3, 'M50 value')
mtext(side = 4, line = 3, 'Number of modules')
legend(x = .6,y = 90,legend=c("M50", "N modules"),
       pch=c(20, 16), col=c("purple3", "black"))

#Identify Thresholds that give us various clustering
thresholds <- T.M50[which(T.M50[,4]>10 &  T.M50[,4]<= 100),]
if(is.null(dim(thresholds))){
    thresholds <- thresholds[1]
}else{
        thresholds <- thresholds[match(unique(thresholds[,4]), thresholds[,4]),][,1]
}    

#For the chosen set of thresholds get all of the clusters present at each threshold
my.cluster.sets <- list()
for(t in 1:length(thresholds)){
   
    grph.cpy <- delete.edges(grph,which(abs(E(grph)$weight) < thresholds[t]))
    my.cluster.sets[[t]] <- decompose.graph(grph.cpy,min.vertices = 2)
    names(my.cluster.sets)[t] <- thresholds[t]
    
}





saveRDS(my.cluster.sets,"./output/AtlasSecondOrderModuleAnalysis/SAL_modules")
#use these clusters to test various hypotheses
###################################################################################################
###################################################################################################
#Now select graph subset based on the complement of the mutual intersect
#grph.opt <- all.healthy.weighted.graph$opt.icov
grph.opt@Dimnames <- list(dimnames(z.healthy)[[2]],dimnames(z.healthy)[[2]])
set.seed(seed = 20170324)
subthresh.miRNA <- (setdiff(colnames(z.healthy),zzz.intersect))
subthresh.miRNA <- sample(subthresh.miRNA,101,replace = F)
grph.intsct.cmplmnt <-(grph.opt[subthresh.miRNA,subthresh.miRNA])
hist((grph.intsct.cmplmnt@x),100)

grph.intsct.cmplmnt@Dimnames[[1]][grph.intsct.cmplmnt@i[which(grph.intsct.cmplmnt@x > .1 & grph.intsct.cmplmnt@x < 10 )]]
grph.intsct.cmplmnt.ig <- graph.adjacency(adjmatrix = grph.intsct.cmplmnt,mode = "upper",weighted = T,diag =  F)
num.edges <- length(E(grph.intsct.cmplmnt.ig)$weight)

thr <- quantile( abs(E(grph.intsct.cmplmnt.ig)$weight),probs = seq(0,1,.001))
T.M50.cmp <-  matrix(nrow = length(thr),ncol = 4,
                 dimnames = list(NULL,c("threshold_weight","M50a","M50b","Number_of_clusters")) ) 
for(cur.thr in thr){
    
    #Definition of  M50: 50% of nodes are in modules of this size or greater  
    grph.cpy <- delete.edges(grph.intsct.cmplmnt.ig,which(abs(E(grph.intsct.cmplmnt.ig)$weight) < cur.thr))
    num.c <- no.clusters(grph.cpy)
    c <- sort(components(grph.cpy)$csize,decreasing = T)
    
    l <- c()
    for(mod in c){l <- c(l,rep(mod,mod))}
    m50a <- median(l)
    m <- sum(c)/2
    m50b<- which(cumsum(c) > m)[1]
    #print(c)[1:15]
    idx <- which(thr == cur.thr)
    T.M50.cmp[idx,] <- c(cur.thr,m50a,m50b,num.c)
    
}

par(mar = c(5,5,5,7))
plot(x = T.M50.cmp[,1],y = T.M50.cmp[,2],xlab = "Minimum abslolute edge weight for inclusion",
     ylab =  NA, ylim = c(0,100),
     main = "M50 vs. edge-weight threshold\nfor clusters based on complement of 101 co-expressed miRNA\n(CSF,Saliva,Serum,Plasma)",
     sub = "(biofluid mean RPM < 15 for all miRNA)",
     col = c("purple2"), cex = .8, pch = 20)
points(x = T.M50.cmp[,1], y = T.M50.cmp[,4], col = "black", cex = .5, pch = 16 )
axis(side = 4)
mtext(side = 2, line = 3, 'M50 value')
mtext(side = 4, line = 3, 'Number of modules')
legend(x = 1.0,y = 900,legend=c("M50", "N modules"),
       pch=c(20, 16), col=c("purple3", "black"))

#Identify Thresholds that give us various clustering
thresholds <- T.M50.cmp[which(T.M50.cmp[,4]>10 &  T.M50.cmp[,4]<= 40),]
thresholds <- thresholds[match(unique(thresholds[,4]), thresholds[,4]),][,1]
###################################################################################################
#For the chosen set of thresholds get all of the clusters present at each threshold
my.cluster.sets.cmplmt <- list()
for(t in 1:length(thresholds)){ 
    grph.cpy.cmplmt <- delete.edges(grph.intsct.cmplmnt.ig,which(abs(E(grph.intsct.cmplmnt.ig)$weight) < thresholds[t]))
    my.cluster.sets.cmplmt[[t]] <- decompose.graph(grph.cpy.cmplmt,min.vertices = 2)
    names(my.cluster.sets.cmplmt)[t] <- thresholds[t]
    
}

saveRDS(my.cluster.sets.cmplmt,"./output/test_clusters_cmplmt_30_thresholds_2")
#use these clusters to test various hypotheses
###################################################################################################
###################################################################################################
###################################################################################################
wkt.com <-(walktrap.community(grph.cpy))
edge.betweenness.community(grph.cpy)
cluster.distribution(grph.cpy)
components(grph.cpy)
no.clusters(grph.cpy)

write.graph(graph = grph.intsct.gml,file = "./output/common_101_miRNA_trsh_15.uppr.graphml", format = "graphml")
write.graph(graph = grph.cpy,file = "./output/jnk.graphml", format = "graphml")




###################################################################################################
my.list <- vector("list",length(levels(z.healthy$biofluid_name)))
names(my.list) <- levels(z.healthy$biofluid_name)
if(F){
for(bf in levels(z.healthy$biofluid_name)){
        datfr <- z.healthy.npn[z.healthy$biofluid_name == bf,]
        if(!any(dim(datfr)== 0 )){
            datfr <- datfr[,which(round(apply(datfr, 2, sd), 4) != 0)] #Ensure that no column is empty 
            print(bf)
            print(dim(datfr))
            my.list[[bf]]$sample_names <- rownames(datfr)
            my.list[[bf]]$rna_names <- colnames(datfr)
            datfr <- huge.npn(datfr) #Perform NPN normalization of the subset
            hist(datfr[,3])
            my.list[[bf]] <- huge(datfr,method = "glasso")
        }
        rm(datfr)
        
}
saveRDS(my.list,"./npnGlassoGraphsHealthyByBiofluid")
}

if(F){
    ###THIS USES BIC METHOD####
my.ebic.optimized.list <- vector("list", length(my.list))
names(my.ebic.optimized.list) <- names(my.list)
for(x in names(my.list)){
    if(!is.null(my.list[[x]])){
        
        opt <- huge.select(est = my.list[[x]])
        
        opt$opt.icov@Dimnames[[1]]  <- dimnames(opt$data)[[2]]
        opt$opt.icov@Dimnames[[2]]  <- dimnames(opt$data)[[2]]
        
        my.ebic.optimized.list[[x]][["huge.selected.est"]] <- opt
        ig <- graph.adjacency(adjmatrix = as.matrix(opt$opt.icov),
                              mode = "undirected", weighted = TRUE)
        my.ebic.optimized.list[[x]][["igraph.est"]] <- ig
    }
    
}
saveRDS(my.ebic.optimized.list,"./eBIC_optimalGraphs_ExRNA_Atlas_Healthy_miRNA")
}

if(F){
###THIS USEs STARS METHOD####
my.stars.optimized.list <- vector("list", length(my.list))
names(my.stars.optimized.list) <- names(my.list)
for(x in names(my.list)){
    if(!is.null(my.list[[x]])){
        x
        my.list[[x]]$lambda
        my.list[[x]]$sparsity
        my.list[[x]]$path
        opt <- huge.select(est = my.list[[x]],criterion = "stars",rep.num = 20, stars.thresh = .05)
        #min(opt$sparsity[which((opt$sparsity > .10 & opt$sparsity < .30))])
        
        #opt$opt.icov@Dimnames[[1]]  <- dimnames(opt$data)[[2]]
        #opt$opt.icov@Dimnames[[2]]  <- dimnames(opt$data)[[2]]
        
        my.stars.optimized.list[[x]][["huge.selected.est"]] <- opt
        ig <- graph.adjacency(adjmatrix = as.matrix(opt$opt.icov),
                              mode = "undirected", weighted = TRUE,
                              add.colnames = dimnames(opt$data)[[2]],
                              add.rownames = dimnames(opt$data)[[2]])
        my.stars.optimized.list[[x]][["igraph.est"]] <- ig
    }

}
saveRDS(my.stars.optimized.list,"./stars_optimalGraphs_ExRNA_Atlas_Healthy_miRNA")
}
####################################################################################
#Subset the graph based on expression threshold

my.ebic.optimized.list <- readRDS("./output/AtlasSecondOrderModuleAnalysis/eBIC_optimalGraphs_ExRNA_Atlas_Healthy_miRNA")

tail(colnames(z.healthy))
s <- by(data = z.healthy,INDICES = z.healthy.meta$biofluid_name,FUN = colMeans, na.rm = T)
healthy.mean.expression <- colMeans(z.healthy)
length(healthy.mean.expression[healthy.mean.expression > 100])

est.igCSF <- my.ebic.optimized.list$`Cerebrospinal fluid`$igraph.est
est.igSER <- my.ebic.optimized.list$Serum$igraph.est
est.igPLS <- my.ebic.optimized.list$Plasma$igraph.est
est.igSAL <- my.ebic.optimized.list$Saliva$igraph.est

#identify the number of edges shared among the 4 graphs
xx <- simplify(Reduce(intersection, list(est.igCSF, est.igSER, est.igPLS, est.igSAL) ))
#Idenify whch vertices in the shared edges come from the 101 @ 15RPM
core.nodes <- intersect(ends(graph = xx,es = E(xx)),zzz.intersect)
length(core.nodes)

est.igCSF.ew <- as.data.frame(cbind( get.edgelist(est.igCSF) , round( E(est.igCSF)$weight, 3 )))
est.igSER.ew <- as.data.frame(cbind( get.edgelist(est.igSER) , round( E(est.igSER)$weight, 3 )))
est.igPLS.ew <- as.data.frame(cbind( get.edgelist(est.igPLS) , round( E(est.igPLS)$weight, 3 )))
est.igSAL.ew <- as.data.frame(cbind( get.edgelist(est.igSAL) , round( E(est.igSAL)$weight, 3 )))

#Convert factor values to numeric where appropriate
est.igCSF.ew$V3 <- as.numeric(as.character(est.igCSF.ew$V3))
est.igSER.ew$V3 <- as.numeric(as.character(est.igSER.ew$V3))
est.igPLS.ew$V3 <- as.numeric(as.character(est.igPLS.ew$V3))
est.igSAL.ew$V3 <- as.numeric(as.character(est.igSAL.ew$V3))

est.igCSF.ew <- est.igCSF.ew[order(est.igCSF.ew[,3],decreasing = T),]
est.igSER.ew <- est.igSER.ew[order(est.igSER.ew[,3],decreasing = T),] 
est.igPLS.ew <- est.igPLS.ew[order(est.igPLS.ew[,3],decreasing = T),] 
est.igSAL.ew <- est.igSAL.ew[order(est.igSAL.ew[,3],decreasing = T),] 

#NEXT REMOVE ALL SELF EDGES
est.igCSF.ew <- est.igCSF.ew[est.igCSF.ew[,1] != est.igCSF.ew[,2],]
est.igSER.ew <- est.igSER.ew[est.igSER.ew[,1] != est.igSER.ew[,2],]
est.igPLS.ew <- est.igPLS.ew[est.igPLS.ew[,1] != est.igPLS.ew[,2],]
est.igSAL.ew <- est.igSAL.ew[est.igSAL.ew[,1] != est.igSAL.ew[,2],]

hist(est.igCSF.ew$V3)
hist(est.igSER.ew$V3)
hist(est.igPLS.ew$V3)
hist(est.igSAL.ew$V3)


#Find a way to uniquely identify edges irrespective of vertex order 
my.get.ordered.edge.names <- function(graph.ew){
    a <- graph.ew
    aa <- apply(a[,1:2],1,naturalsort)
    aaa <- apply(aa,2,function(x){paste(x[1],x[2],sep = "<-->")})
    a$EdgeNames <- aaa
    
    return(a)
}

a <- my.get.ordered.edge.names(est.igCSF.ew)
b <- my.get.ordered.edge.names(est.igSER.ew)
c <- my.get.ordered.edge.names(est.igPLS.ew)
d <- my.get.ordered.edge.names(est.igSAL.ew)

a <- a[order(abs(a$V3),decreasing = T),]
b <- b[order(abs(b$V3),decreasing = T),]
c <- c[order(abs(c$V3),decreasing = T),]
d <- d[order(abs(d$V3),decreasing = T),]
tail(a)

edges.list <- list(a$EdgeNames,b$EdgeNames,c$EdgeNames,d$EdgeNames)
names(edges.list) <- c("CSF","SERUM","PLASMA","SALIVA")
el <- calculate.overlap(edges.list)
venn.diagram(edges.list,"./output/AtlasSecondOrderModuleAnalysis/Venn_by_biofluid.png",
             imagetype = "png", main = "Shared miRNA Modules for 4 biofluids in ExRNA Atlas",
             fill=c("darkmagenta", "blue", "yellow2", "orange2"), cex = 1.5,
             cat.fontface = 4, cat.cex = 1.5, category.names = names(edges.list))

a.top <- a[order(abs(a$V3),decreasing = T),][1:500,]
b.top <- b[order(abs(b$V3),decreasing = T),][1:500,]
c.top <- c[order(abs(c$V3),decreasing = T),][1:500,]
d.top <- d[order(abs(d$V3),decreasing = T),][1:500,]

edges.list.top <- list(a.top$EdgeNames,b.top$EdgeNames,c.top$EdgeNames,d.top$EdgeNames)
names(edges.list.top) <- c("CSF","SERUM","PLASMA","SALIVA")
el <- calculate.overlap(edges.list.top)
venn.diagram(edges.list.top,"./output/AtlasSecondOrderModuleAnalysis/VennTop500Edges_by_biofluid.png",
             imagetype = "png", main = "Shared miRNA Modules for 4 biofluids in ExRNA Atlas",
             fill=c("darkmagenta", "blue", "yellow2", "orange2"), cex = 1.5,
             cat.fontface = 4, cat.cex = 1.5, category.names = names(edges.list))

a.101 <- a[order(abs(a$V3),decreasing = T),][which(a$V1 %in% zzz.intersect & a$V2 %in% zzz.intersect),]
b.101 <- b[order(abs(b$V3),decreasing = T),][which(b$V1 %in% zzz.intersect & b$V2 %in% zzz.intersect),]
c.101 <- c[order(abs(c$V3),decreasing = T),][which(c$V1 %in% zzz.intersect & c$V2 %in% zzz.intersect),]
d.101 <- d[order(abs(d$V3),decreasing = T),][which(d$V1 %in% zzz.intersect & d$V2 %in% zzz.intersect),]

edges.list.101 <- list(a.101$EdgeNames,b.101$EdgeNames,c.101$EdgeNames,d.101$EdgeNames)
names(edges.list.top) <- c("CSF","SERUM","PLASMA","SALIVA")
el <- calculate.overlap(edges.list.top)
venn.diagram(edges.list.101,"./output/AtlasSecondOrderModuleAnalysis/VennAllSharedEdges_by_biofluid_for101miRNA.png",
             imagetype = "png", main = "Shared miRNA Modules for 4 biofluids in ExRNA Atlas",
             fill=c("darkmagenta", "blue", "yellow2", "orange2"), cex = 1.5,
             cat.fontface = 4, cat.cex = 1.5, category.names = names(edges.list))



names(zzz.names) <- names(edges.list.top) 
RPM15_miRNA <- calculate.overlap(zzz.names)
lapply(RPM15_miRNA,length)

venn.diagram(zzz.names, "./output/AtlasSecondOrderModuleAnalysis/Venn15RPM_4biofluids.png",
             imagetype = "png", main = "Shared miRNA for 4 biofluids in ExRNA Atlas",
             fill=c("darkmagenta", "blue", "yellow2", "orange2"), cex = 1.5,
             cat.fontface = 4, cat.cex = 1.5, category.names = names(zzz.names))



est.igCSF.ew$EdgeNames <- est.i
b <- (est.igPLS.ew)
aa <- apply(a[,1:2],1,naturalsort)
aaa <- apply(aa,2,function(x){paste(x[1],x[2],sep = "<-->")})
bb <- apply(b[,1:2],1,naturalsort)
bbb <- apply(bb,2,function(x){paste(x[1],x[2],sep = "<-->")})
a$Pairs <- aaa
b$Pairs <- bbb

c <- my.get.ordered.edge.names(est.igSAL.ew)
d <- my.get.ordered.edge.names(est.igSER.ew)

write.graph(est.igCSF,"./output/igCSF.graphml",format = "graphml")
