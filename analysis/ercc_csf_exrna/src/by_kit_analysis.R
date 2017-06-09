#Input file is in the Genboree Workbench under the group: 'UH2_CSF_Supplement'
#An archive file for the exceRpt output can be found in 'Files' under the
#Data base for 'Small RNA-seq Data'
#file structure produced by exceRpt.
#Download the file: 'UH2_CSF_Supp_smallRNA_exceRpt_v4.3.5_Results.zip'


#Create violin plots with
#Required Packages
#install.packages('VennDiagram')
#library(VennDiagram)
library(lattice)

#library(plyr)
#library(reshape)
source('./Rscripts/UH2.CSF.functions.0.2.R',echo = F)

#Point to the directory containing the result files for each sample
df.pathRoot <- './input/exceRpt_output/UH2_CSF_Supp_smallRNA_exceRpt_v4.3.5_CoreResults/CORE/'

####################################################################################################
#Once through before adding loops to cover all cases...(Turn this into a funciton call):
####################################################################################################
#Create list structure:
#Here I will analyze variability within RNA species
#so I choose this as the top level 
RKC <- populate.lst(targ.lst = UH2.CSF.R.K.C.rdct,df.path = df.pathRoot)

by.R.type_kits.count <- list()
for( R.type in names(RKC) ){
    #For Each RNA type: For each kit get all of the readcount data 
    kits.count <- list()
    for( K.type in names(RKC[[R.type]]) ){
    all.cond.counts.lst <- RKC[[R.type]][[K.type]]
        #We now have the list of data for conditions for a given kit and RNA type
        #We want to compile the data for read counts in each condition so first
        #we must identify all of the detected species
        #Then create a data frame with counts of each RNA in each condition
        unq <- c()
        all.df <- data.frame()
        ubiq.df <- data.frame()
        for( cond in names(RKC[[R.type]][[K.type]]) ){
            cond.df <- all.cond.counts.lst[[cond]]
            if(is.null(cond.df)){next}
            rnames <- rownames(cond.df)
            cond.df <- cond.df$uniqueReadCount
            names(cond.df) <- rnames
            cond.df <- as.data.frame(cond.df)
            colnames(cond.df) <- cond
            unq <- unique(c(unq,rnames))
            all.df <- merge.data.frame(cond.df,all.df,by = 0,all = T)
            rownames(all.df) <-all.df$Row.names
            all.df <- all.df[,c(2:ncol(all.df)),drop =F]
            if( all(colnames(all.df) %in% names(RKC[[R.type]][[K.type]])) ){
                ubiq.df <- all.df[complete.cases(all.df),]}
            
        }
        kits.count[[K.type]] <- list(all = all.df, ubiq = ubiq.df)
    }
    
    by.R.type_kits.count[[R.type]] <- kits.count  
}
#saveRDS(object = by.R.type_kits.count, file = "./output/by.R.type_kits.count.Rdat" )

####################################################################################################
#Get summary of diveristy of RNA species detected by each kit
#RNATYPE: mirvana detected,mirvana ubiq,ExorDetected, exor ubiq, Both detected, Total Detected
####################################################################################################
R.type.summary <- list()
for( R.type in names(RKC) ){
    K.nums <- list()
    for( K.type in names(RKC[[R.type]]) ){
        a.names <- rownames(by.R.type_kits.count[[R.type]][[K.type]]$all)
        u.names <- rownames(by.R.type_kits.count[[R.type]][[K.type]]$ubiq)
        a.cnt <- length(a.names)
        u.cnt <- length(u.names)
        K.nums[[K.type]] <- list(all = a.names, ubiq = u.names)
        
        
    }
    tot.cnt <-  length(unique(c(K.nums$exoRneasy$all,K.nums$miRvana$all)))
    inx.a.cnt <- length(intersect(K.nums$exoRneasy$all,K.nums$miRvana$all))
    inx.u.cnt <- length(intersect(K.nums$exoRneasy$ubiq,K.nums$miRvana$ubiq))
    K.nums <- list( ExoR = lapply(K.nums$exoRneasy,length),
                    MiRv = lapply(K.nums$miRvana,length))
    
    R.type.summary[[R.type]] <-list(total.count = tot.cnt,
                                    all.intersect = inx.a.cnt,
                                    ubiq.intersect = inx.u.cnt,
                                    kit.counts = K.nums
                                    )
    
}

zz <-data.frame()
for(x in names(RKC)){
    z <- unlist(R.type.summary[[x]])
    zz <- rbind(zz,z)
    
}
rownames(zz) <- names(RKC)
colnames(zz) <- c(names(unlist(R.type.summary[[x]])))
write.table(x = zz,file = './output/diveristy.txt',quote = F,
            sep = '\t',col.names = T)
####################################################################################################
# Where Possible get the correlations between kits
####################################################################################################


my.R.types <- c(5:8,10) 
pdf(file = './output/junkplot.pdf',width = 20,height = 20,onefile = T)
par(oma = c(10, 10, 0, 0)) # make room (i.e. the 4's) for the overall x and y axis titles

par(mar= c(rep(0.5,4))) 
par(mfrow = c(length(RKC[my.R.types]), length(RKC[[1]][[1]])) )

for( R.type in names(RKC)[my.R.types] ){

    #For miRNA detected in ANY condition
    exoR.all <- by.R.type_kits.count[[R.type]]$exoRneasy$all
    miRv.all <- by.R.type_kits.count[[R.type]]$miRvana$all
    exoR.all.names <- rownames(exoR.all)
    miRv.all.names <- rownames(miRv.all)
    all.cmp.names <- intersect(exoR.all.names,miRv.all.names)
    if( length(all.cmp.names) <1 ){
        plot(x = 1, y = 1, col = 'red', pch = 12, cex = 8)
        next
    }
    exoR.cmp <- exoR.all[all.cmp.names,]
    miRv.cmp <- miRv.all[all.cmp.names,]
    for(cond in intersect(colnames(exoR.cmp),colnames(miRv.cmp)) ){
 
        ex <- !is.na(exoR.cmp[,cond])
        mi <- !is.na(miRv.cmp[,cond])
        exoR.cmp.cond <- exoR.cmp[ex,cond,drop = F]
        miRv.cmp.cond <- miRv.cmp[mi,cond, drop =F]    
        
        
        shared.by.cond <- merge(exoR.cmp.cond, miRv.cmp.cond,by = 0)
        colnames(shared.by.cond) <- c('RNA.type',names(RKC[[R.type]]))
        print(cor(shared.by.cond$exoRneasy,shared.by.cond$miRvana))
        
        exl <- length(exoR.cmp.cond[[cond]])
        mvl <- length(miRv.cmp.cond[[cond]])
        shl <- dim(shared.by.cond)[1]
        if( shl < 2 ){
            plot(x = 1, y = 1, col = 'blue', pch = 12, cex = 8)
            next
        }
        
        plt.lim <- max(shared.by.cond$exoRneasy,shared.by.cond$miRvana)
        plot(x = shared.by.cond$exoRneasy, y = shared.by.cond$miRvana,#,
             #xlab = names(RKC[[1]][1]), ylab = names(RKC[[1]][2]),
             xlim = c(0,plt.lim), ylim = c(0,plt.lim) )#,
             #main = paste(cond, ':', R.type, '\n', paste(shl,"shared of",exl,"&",mvl)),
             #sub = paste("Pearson:",round(cor(shared.by.cond$exoRneasy,shared.by.cond$miRvana),2)) )
        ls1 <- lsfit(shared.by.cond$exoRneasy,shared.by.cond$miRvana)
        ls2 <- lsfit(shared.by.cond$miRvana,shared.by.cond$exoRneasy)
        abline(ls1,col="red")
        abline(ls2,col="purple")
        
        #abline(lm(exoRneasy ~ miRvana,data = shared.by.cond), col="blue", lty =2)
        #abline(lm(miRvana ~ exoRneasy,data = shared.by.cond), col="cyan", lty = 2)
        abline(0,1)
        ##POSSIBLY FILL RED IF LEAST SQUARES LINE IS ABOVE AND BLUE IF BELOW (USE TRANSPARENCY)
        
    }
    
    #For miRNA detected in all conditions
    
   
}
# print the overall labels
mtext('x-axis title', side = 1, outer = TRUE, line = 2)
mtext('y-axis title', side = 2, outer = TRUE, line = 2)
dev.off()

#Also try lattic dotplot with mirvana predicting exornaeasy conditioning on e.g. R type and/or condition 
dotplot(Scheme ~ Area_under_ROC | Dataset, data = simulationSummary, layout = c(4,6))

####################################################################################################

if(F){
#Get the correlations for each mRNA type by kit and condition 
n <- by.R.type_kits.count$miRNAmature_sense$exoRneasy$all
m <- by.R.type_kits.count$miRNAmature_sense$miRvana$all
i <- intersect(rownames(n),rownames(m))
n <- n[i,]
m <- m[i,]
n <- n[complete.cases(n$SAH),]
m <- m[complete.cases(m$SAH),]
i <- intersect(rownames(n),rownames(m))
cor(n[i,]$SAH,m[i,]$SAH)
plot(n[i,]$SAH,m[i,]$SAH)

#Subset  by Kit:
tmp <- populate.lst(targ.lst = UH2.CSF.R.K.C.rdct,df.path = df.pathRoot)
jnk <- tmp$miRNAmature_sense$exoRneasy$AD$uniqueReadCount
jnk1 <- tmp$miRNAmature_sense$exoRneasy$Control$uniqueReadCount
names(jnk) <- rownames(tmp$miRNAmature_sense$exoRneasy$AD)
names(jnk1) <- rownames(tmp$miRNAmature_sense$exoRneasy$Control)
bb <-data.frame()
as.data.frame(jnk)
as.data.frame(jnk1)
ccd <- merge.data.frame(bb,jnk1,by = 0,all = T)
head(ccd)
ccd <- merge.data.frame(ccd,jnk,by = 0,all = T)
Exor.miRmat <- tmp$miRNAmature_sense$exoRneasy
Mirv.miRmat <- tmp$miRNAmature_sense$miRvana

#Combine the data to facilitate comparisons
#We start by listing all the unique mature miRNA species identified
#(by rownames for unique reads > 0) 
total.mirna.Exor <- c(rownames(Exor.miRmat$AD),rownames(Exor.miRmat$Control),
                      rownames(Exor.miRmat$GBM),rownames(Exor.miRmat$LGG),
                      rownames(Exor.miRmat$PD),rownames(Exor.miRmat$SAH))
total.mirna.Exor <- unique(total.mirna.Exor)
total.mirna.Mirv <- c(rownames(Mirv.miRmat$AD),rownames(Mirv.miRmat$Control),
                      rownames(Mirv.miRmat$GBM),rownames(Mirv.miRmat$LGG),
                      rownames(Mirv.miRmat$PD),rownames(Mirv.miRmat$SAH))
total.mirna.Mirv <- unique(total.mirna.Mirv)

total.miRna <- unique(c(total.mirna.Exor,total.mirna.Mirv))
total.intersect <- intersect(total.mirna.Mirv,total.mirna.Exor)



#Now that we know how many unique RNA we have, construct a data frame 
my.tst.1 <- cbind(Mirv.miRmat$AD$uniqueReadCount,rownames(Mirv.miRmat$AD))
my.tst.2 <- cbind(Mirv.miRmat$GBM$uniqueReadCount,rownames(Mirv.miRmat$GBM))
head(my.tst.2)

tst <- merge.data.frame(my.tst.1, my.tst.2,by = 2,all = T)
tst[is.na(tst)] <- 0
tail(tst,30)
length(unique(c(rownames(Mirv.miRmat$AD),rownames(Mirv.miRmat$GBM))))
plot(Mirv.miRmat$AD$uniqueReadCount,Mirv.miRmat$AD$uniqueReadCount)


a <- intersect(rownames(Exor.miRmat$AD),rownames(Exor.miRmat$Control))
b <- intersect(rownames(Exor.miRmat$GBM),rownames(Exor.miRmat$LGG))
c <- intersect(rownames(Exor.miRmat$PD),rownames(Exor.miRmat$SAH))
d <- intersect(a,b)
e <- intersect(c,d)



aa <- intersect(rownames(Mirv.miRmat$AD),rownames(Mirv.miRmat$Control))
ba <- intersect(rownames(Mirv.miRmat$GBM),rownames(Mirv.miRmat$LGG))
ca <- intersect(rownames(Mirv.miRmat$PD),rownames(Mirv.miRmat$SAH))
da <- intersect(aa,ba)
ea <- intersect(ca,da)

b1 <- intersect(rownames(Mirv.miRmat$AD),rownames(Exor.miRmat$AD))
b2 <- intersect(rownames(Mirv.miRmat$Control),rownames(Exor.miRmat$Control))
b3 <- intersect(rownames(Mirv.miRmat$GBM),rownames(Exor.miRmat$GBM))
b4 <- intersect(rownames(Mirv.miRmat$LGG),rownames(Exor.miRmat$LGG))
b5 <- intersect(rownames(Mirv.miRmat$PD),rownames(Exor.miRmat$PD))
b6 <- intersect(rownames(Mirv.miRmat$SAH),rownames(Exor.miRmat$SAH))
length(b1)
xyplot(uniqueReadCount ~ multimapAdjustedReadCount, data = Mirv.miRmat$AD)
cor(Mirv.miRmat$AD$uniqueReadCount,Mirv.miRmat$AD$multimapAdjustedReadCount) 

length(e)
length(ea)
length(intersect(e,ea))
####################################################################################################
#To directly compare the split sample by kit
#Create correlation plots (make sure rows are ordered by RNA species name)
#For each type of RNA e.g. circ, mRNA, miRNA

#write a function that imports the file each type of RNA mapping
#and appropriately formats the data
#Create a list structure to organize the samples by: RNA TYPE > KIT > CONDITION




#EXAMPLE FLOW
#Select a set of raw read count DFs and add it to a list for comparison
GBM.miR.mat <- list(exoR = UH2.CSF.readcounts$exoRneasy$GBM$miRNA.mature,
                    miRv = UH2.CSF.readcounts$miRvana$GBM$miRNA.mature)


#get intersects
intersect <- intersect(rownames(GBM.miR.mat$exoR),
                       rownames(GBM.miR.mat$miRv))
pct_a <- length(intersect)/nrow(GBM.miR.mat$exoR)
pct_b <- length(intersect)/nrow(GBM.miR.mat$miRv)

exoR.shared <- GBM.miR.mat$exoR[(rownames(GBM.miR.mat$exoR) %in% intersect),]
miRv.shared <- GBM.miR.mat$miRv[(rownames(GBM.miR.mat$miRv) %in% intersect),]
exoR.shared <- exoR.shared[order(rownames(exoR.shared)),]
miRv.shared <- miRv.shared[order(rownames(miRv.shared)),]

rownames(exoR.shared) == rownames(miRv.shared)
cor(exoR.shared$uniqueReadCount,miRv.shared$uniqueReadCount)

my.sub.titles = c("GBM mature miRNA")
plot(exoR.shared$uniqueReadCount,miRv.shared$uniqueReadCount,sub = my.sub.titles[1])
cor(exoRneasy$uniqueReadCount,miRvana$uniqueReadCount)    

grid.newpage()
draw.pairwise.venn(area1 = , area2 = 20, cross.area = 11, category = c("Dog People", 
                                                                       "Cat People"))




AD1_readCounts_miRNAmature_sense.txt <- read.delim("./input/by_kit/exoRNeasy/AD_1/readCounts_miRNAmature_sense.txt")
AD2_readCounts_miRNAmature_sense.txt <- read.delim("./input/by_kit/miRVana/AD_2/readCounts_miRNAmature_sense.txt")
}