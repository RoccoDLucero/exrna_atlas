#Hi Rocco,
#Please focus on the graph in slide 4 that is based on 124 edges shared by all biofluids.
#For each of the 15 modules in the graph that have 3 or more nodes, do the following:
#   1. For all miRNAs within the module construct a vector of expression levels across tissues.
#   2. Calculate correlations of the vectors for all pairs of miRNAs within the module.
#   3. Calculate pairwise correlations of miRNAs within the module with miRNAs from other 14 modules.
#   4. Plot distributions in (2) and (3) and perform a test to examine any difference between (2) vs (3)
#The hypothesis here is that some modules will have higher intra-correlation than
# inter-correlation. How many (if any) modules show this pattern?
#################################################################################
pr.cors <- function(x,dat = dat1){
    #print(x[1])
    #print(x[2])
    a <- dat[,x[1]]
    b <- dat[,x[2]]
    #print(a)
    #print(b)
    return( list(x[1],x[2],cor(a,b)) )
}



#Load the Cell Atlas Raw miRNA counts and rectify data types.
Cell.Atlas <- read.delim2("./input/Merged_results_joel/exceRpt_miRNA_ReadsPerMillion.txt",
                          sep = "\t",row.names = 1,stringsAsFactors = F)
CA.meta <- read.delim2("./input/Merged_results_joel/exceRpt_sampleGroupDefinitions.txt",
                       sep = "\t",row.names = 1,stringsAsFactors = T)


SRA.Atlas <- read.delim2("./input/Brainspan+SRA_merged_results/exceRpt_miRNA_ReadsPerMillion.txt",
                         sep = "\t",row.names = 1,stringsAsFactors = F)
SRA.meta <-read.delim2("./input/Brainspan+SRA_merged_results/run_annotations_full--fixed.txt",
                       sep = "\t",row.names = 1,stringsAsFactors = T)

Cell.Atlas <- to.numeric.df(tFrame(Cell.Atlas))
SRA.Atlas <- to.numeric.df(tFrame(SRA.Atlas))

Cell.Atlas <- Cell.Atlas[,which(round(apply(Cell.Atlas, 2, sd), 4) != 0)] #Remove all columns that have no variance
SRA.Atlas <- SRA.Atlas[,which(round(apply(SRA.Atlas, 2, sd), 4) != 0)] #Remove all columns that have no variance

QN <- TRUE
if(QN){
    CA.Qnorm <- normalize.quantiles( as.matrix(Cell.Atlas) )
    Cell.Atlas.QN <- as.data.frame(CA.Qnorm)
    dimnames(Cell.Atlas.QN) <- dimnames(Cell.Atlas)
    Cell.Atlas <- Cell.Atlas.QN
    
    SRA.Qnorm <- normalize.quantiles( as.matrix(SRA.Atlas) )
    SRA.Atlas.QN <- as.data.frame(SRA.Qnorm)
    dimnames(SRA.Atlas.QN) <- dimnames(SRA.Atlas)
    SRA.Atlas <- SRA.Atlas.QN
}

NPN <- TRUE

if(NPN){
    Cell.Atlas <- huge.npn(Cell.Atlas)
    
    SRA.Atlas <- huge.npn(SRA.Atlas)
    
}

################################################################################

gph.file <- "./output/AtlasSecondOrderModuleAnalysis/grph.4biofluids.intersect.simplified"
grph.4biofluids.intersect.simplified <- readRDS(gph.file)
cmpnts <- clusters(grph.4biofluids.intersect.simplified)
top15.mods <- decompose(grph.4biofluids.intersect.simplified,min.vertices = 3)

#Set the input arguments and parameters
mods <- top15.mods
dat <-SRA.Atlas
brks <- 20
txp <- .3
par(mfrow = c(2,2))
combined.corrs <- list(in.ex = c(),in.in = c(), true.edges = c())
for(md in mods){
    #   1. For all miRNAs within the module construct a vector of expression levels across tissues.
    nds <- intersect(names(V(md)), colnames(dat))
    md <-  delete.vertices(md, setdiff(names(V(md)),nds ))
    if(length(E(md))==0){next}
    un.nds <- setdiff(colnames(dat),nds)
    dat1 <- dat[,nds]
    dat2 <- dat[,un.nds]
    #   2. Calculate correlations of the vectors for all pairs of miRNAs within the module.
    md.pairs <- as.matrix(ends(md,E(md)))
    im.edge.cors  <- apply(md.pairs,1,pr.cors, dat = dat1)
    im.edge.cors.vals <- as.numeric(unlist(im.edge.cors)[seq(from = 3,to = length(unlist(im.edge.cors)),by = 3)])
    im.all.cors <- as.numeric(combn(nds,2,FUN = function(x){cor(dat1[,x[1]], dat1[,x[2]])}))
     
    
    #   3. Calculate pairwise correlations of miRNAs within the module with miRNAs from other 14 modules.
    ex.md.mirs <- names(unlist(sapply(mods,V)))
    ex.md.mirs <- intersect(ex.md.mirs, colnames(dat2))
    in.ex.cors <- c()
    for(nd in nds){
        for(ex.nd in ex.md.mirs){
            a <- dat1[,nd]
            b <- dat2[,ex.nd]
            in.ex.cors <- c(in.ex.cors, cor(a,b))
        }
    }
    
    combined.corrs$in.ex <- c(combined.corrs$in.ex, in.ex.cors)
    combined.corrs$in.in <- c(combined.corrs$in.in, im.all.cors)
    combined.corrs$true.edges <- c(combined.corrs$true.edges, im.edge.cors.vals)
    
    #   4. Plot distributions in (2) and (3) and perform a test to examine any difference between (2) vs (3)
    if(F){

    ttl <- paste("Module size = ", length(nds))
    hist(in.ex.cors,breaks = brks,xlim = c(-1,1), col = rgb(0,0,1,txp), main = ttl)
    hist(im.all.cors, add = T, col = rgb(1,0,0,txp), breaks = brks)
    hist(im.edge.cors.vals, add = T, col = rgb(0,1,0,txp), breaks = brks)
    legend(x = 'topleft',legend = c('in/out', 'in/in','in/in true edges'),
           fill = c(rgb(0,0,1,txp),rgb(1,0,0,txp),col = rgb(0,1,0,txp)) )
    
    permTS(x = in.ex.cors, y = im.edge.cors.vals)
    
    }
}

prmtst.im.all  <-permTS(x = combined.corrs$in.ex, y = combined.corrs$in.in, alternative = "less")
prmtst.true.edge <- permTS(x = combined.corrs$in.ex, y = combined.corrs$true.edges, alternative = "less")


txp <- .5
ttl <- paste("All modules combined:\nmiRNA expression correlation in SRA samples")
plot(density(combined.corrs$true.edges), lwd = 4,col = rgb(0,1,0,txp), main = ttl)
lines(density(combined.corrs$in.ex),lwd = 4,col = rgb(0,0,1,txp))
lines(density(combined.corrs$in.in),lwd = 4,col = rgb(1,0,0,txp))
legend(x = 'topleft',legend = c('in/out', 'in/in','in/in true edges'),
       fill = c(rgb(0,0,1,txp),rgb(1,0,0,txp),col = rgb(0,1,0,txp)) )
text.default("topright",labels = c(paste("in/out permutation tes p:",),)

par(mfrow = c(1,1))
ttl <- paste("All modules combined")
hist(combined.corrs$in.ex,breaks = brks,xlim = c(-1,1), col = rgb(0,0,1,txp), main = ttl, freq = F)
hist(combined.corrs$in.in, add = T, col = rgb(1,0,0,txp), breaks = brks,freq = F)
hist(combined.corrs$true.edges, add = T, col = rgb(0,1,0,txp), breaks = brks, freq = F)
legend(x = 'topleft',legend = c('in/out', 'in/in','in/in true edges'),
       fill = c(rgb(0,0,1,txp),rgb(1,0,0,txp),col = rgb(0,1,0,txp)) )

