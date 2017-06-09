source("./exrna_atlas_functions.R")

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

SRA.meta <- SRA.meta[order(SRA.meta$tissue),]
SRA.meta <- SRA.meta[intersect(rownames(SRA.meta),rownames(SRA.Atlas)),]
SRA.Atlas <- SRA.Atlas[intersect(rownames(SRA.meta),rownames(SRA.Atlas)),]
SRA.Atlas <- SRA.Atlas[rownames(SRA.meta),]

Cell.Atlas <- Cell.Atlas[,which(round(apply(Cell.Atlas, 2, sd), 4) != 0)] #Remove all columns that have no variance
SRA.Atlas <- SRA.Atlas[,which(round(apply(SRA.Atlas, 2, sd), 4) != 0)] #Remove all columns that have no variance

QN <- T
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

NPN <- F

if(NPN){
    Cell.Atlas <- huge.npn(Cell.Atlas)
    
    SRA.Atlas <- huge.npn(SRA.Atlas)
    
}

#OPTIONAL SUBSETTING
SRA.wt <- list(1,0,"?")
SRA.Atlas <- SRA.Atlas[which(SRA.meta$wt %in% SRA.wt),]
SRA.meta <- droplevels.data.frame(SRA.meta[which(SRA.meta$wt %in% SRA.wt),])

################################################################################
mods <- readRDS("./output/AtlasSecondOrderModuleAnalysis/intersectGraphs/grph_4biofluids_intersect_simplified")
mods <- decompose.graph(mods,min.vertices = 3)
mod.order.sz <- order(sapply(sapply(mods,V),length),decreasing = T)
mods <-mods[mod.order.sz]
input.module <- names(unlist(sapply(mods,V)))
input.colors <- as.list(lapply((lapply(mods,V)),names))
mod.cols.df <-NULL
for(md in 1:15){ 
     a <- (cbind(input.colors[[md]],md) )
     mod.cols.df <- rbind(mod.cols.df,a)
}
mod.cols.df <- as.data.frame(mod.cols.df)
mColors <- with(mod.cols.df,data.frame(md = levels(md),
                                  color = I(coul[1:15])))
CA.side.Colrs <- mColors$color[match(mod.cols.df$md, mColors$md)]


ds <- 2
dat3 <- list(SRA.Atlas, Cell.Atlas)[[ds]]
input.module <- intersect(input.module,colnames(dat3))
non.input.module <- setdiff(colnames(dat3),input.module)
bio.miRNA.inside.module = dat3[,input.module]
bio.miRNA.outside.module = dat3[, non.input.module]

CorMatrix.outside = cor(bio.miRNA.inside.module, bio.miRNA.outside.module, use="complete")
CorMatrix.inside = cor(bio.miRNA.inside.module, bio.miRNA.inside.module, use="complete")
CorMatrix.inside.lower = CorMatrix.inside[lower.tri(CorMatrix.inside)]
CorMatrix.outside.list = as.data.frame(as.table(CorMatrix.outside))
CorMatrix.inside.lower.list = as.data.frame(as.table(CorMatrix.inside.lower))

#png("~/Documents/BCM/Lab/exRNA/miRNA_Table_March272017/Clustering.CSF.Module.1.png",1000,1000)

my_palette <- c(rev(brewer.pal(4,"Reds")),"white",brewer.pal(9,"Blues"))
plot(rep(1,20),col=my_palette,pch=19,cex=3)

CA.clst.seps <- c(2,6,9,37,38,43,50,55,77,104)
test = heatmap.2(as.matrix(CorMatrix.inside), Colv = F,Rowv = F,
                 trace="none",
                 col=   my_palette,#"redblue",# 
                 #colsep = CA.clst.seps , rowsep = (104 - CA.clst.seps),  #c(3,9,27,36,41,50), #
                 sepcolor = "black",
                 breaks = seq(min(CorMatrix.inside)-.05,1,length.out = 15),
                 symbreaks = T, symm = F, symkey = F,
                 RowSideColors = CA.side.Colrs , #rainbow(nrow(mod.cols.df)),
                 labCol = "", labRow = input.module,
                 margins=c(10,10))
if(F){
    test2 = heatmap.2(as.matrix(CorMatrix.outside),
                     trace="none",
                     col= "redblue",
                     #breaks = col_breaks,
                     labCol = FALSE,
                     margins=c(10,10))
}
#dev.off()

ks.test(CorMatrix.inside.lower, CorMatrix.outside)


#colnames(bio.miRNA.inside.module)[test$colInd]
#c(2,6,9,37,38,43,50,55,77,104)
my.CA.clsts <- list(c(2:6), c(9:37), c(50,55), c(55:77), c(77:104))
#my.CA.clsts <- list(c(4:9),c(28:36),c(41:50)) #c(3,9,27,36,41,50)
par(mfrow = c(length(my.CA.clsts)+1,1), xpd = T)
i <-0
for(clstr in rev(my.CA.clsts)){
    i <- i+1
    G1 <- colnames(bio.miRNA.inside.module)[test$colInd][clstr]
    bio.miRNA.inside.module.G1 = bio.miRNA.inside.module[,colnames(bio.miRNA.inside.module) %in% G1]
    if(ds == 1){
        tiss <-levels(SRA.meta$tissue)
        grp <- which(SRA.meta$tissue %in% tiss)
        tmp.df <- SRA.meta[grp,]
        
        pColors <- with(tmp.df,data.frame(tissue = levels(tissue),
                                        color = I(coul[1:length(tiss)])))
        pColrs <- pColors$color[match(tmp.df$tissue, pColors$tissue)]
        
        pColrs <- data.frame(subset(tmp.df, select = c(tissue)),
                   matchRetVal = match(tmp.df$tissue, pColors$tissue))
        
        if(i < length(my.CA.clsts)){
            par(mar = c(0.5,2,0.5,1))
            barplot(t(as.matrix(bio.miRNA.inside.module.G1[grp,])),las = 3,
                    beside = F, col = rainbow( length(G1)),
                    names.arg = gsub("","",SRA.meta$tissue[grp]),
                    cex.names = 1, xaxt = 'n', ann = F
            )
        }else{
            par(mar = c(0.5,2,1,0))
            plt <- barplot(t(as.matrix(bio.miRNA.inside.module.G1[grp,])),las = 3,
                    beside = F, col = rainbow( length(G1)),
                    names.arg = gsub("","",SRA.meta$tissue[grp]),
                    cex.names = 1
            )        
            rect(xleft = plt[1:length(plt)-1],ybottom =-200,xright = plt[2:length(plt)],ytop = -1,
                density = 100,col = pColrs,angle = 0)
        }    
    }else{
        tiss <-levels(CA.meta$sampleGroup)
        grp <- which(CA.meta$sampleGroup %in% tiss)
        if(i < length(my.CA.clsts)){
            par(mar = c(0.5,2,0.5,1))
            barplot(t(as.matrix(bio.miRNA.inside.module.G1[grp,])),las = 3,
                    beside = F, col = "black", #rainbow( length(G1)),
                    names.arg = gsub("Tissue_","",CA.meta$sampleGroup[grp]),
                    cex.names = 1, xaxt = 'n', ann = F
            )
        }else{
            par(mar = c(0.5,2,1,0))
            barplot(t(as.matrix(bio.miRNA.inside.module.G1[grp,])),las = 3,
                    beside = F, col = "black", #rainbow( length(G1)),
                    names.arg = gsub("Tissue_","",CA.meta$sampleGroup[grp]),
                    cex.names = 2
            )        
        }    
    }    
}



library(RColorBrewer)

# Classic palette BuPu, with 4 colors
coul = brewer.pal(12, "Paired") 

# I can add more tones to this palette :
coul = colorRampPalette(coul)(15)

# Plot it
pie(rep(1, length(coul)), col = coul , main="") 


