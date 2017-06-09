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
            cond.df <- cond.df$totalReadCount
            names(cond.df) <- rnames
            cond.df <- as.data.frame(cond.df)
            colnames(cond.df) <- cond
            unq <- unique(c(unq,rnames))
            all.df <- merge.data.frame(cond.df,all.df,by = 0,all = T)
            rownames(all.df) <-all.df$Row.names
            all.df <- all.df[,c(2:ncol(all.df)),drop =F]
            if( all(colnames(all.df) %in% names(RKC[[R.type]][[K.type]])) ){
                ubiq.df <- all.df[complete.cases(all.df),]}
            all.df[is.na(all.df)] <- 0
            
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


my.R.types <- c(5:7,10) 
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
        #print(cor(shared.by.cond$exoRneasy,shared.by.cond$miRvana))
        
        exl <- length(exoR.cmp.cond[[cond]])
        mvl <- length(miRv.cmp.cond[[cond]])
        shl <- dim(shared.by.cond)[1]
        if( shl < 2 ){
            plot(x = 1, y = 1, col = 'blue', pch = 12, cex = 8)
            next
        }
        
        plt.lim <- max(shared.by.cond$exoRneasy,shared.by.cond$miRvana)
        #plot(x = shared.by.cond$exoRneasy, y = shared.by.cond$miRvana,#,
             #xlab = names(RKC[[1]][1]), ylab = names(RKC[[1]][2]),
             #xlim = c(-1,plt.lim), ylim = c(-1,plt.lim) )#,
             #main = paste(cond, ':', R.type, '\n', paste(shl,"shared of",exl,"&",mvl)),
             #sub = paste("Pearson:",round(cor(shared.by.cond$exoRneasy,shared.by.cond$miRvana),2)) )
        plot(x = shared.by.cond$exoRneasy, y = shared.by.cond$miRvana,
             col = ifelse( (shared.by.cond$exoRneasy==0|shared.by.cond$miRvana==0) ,"green","black"),
             xlim = c(-1,plt.lim), ylim = c(-1,plt.lim))
        ls1 <- lm(shared.by.cond$exoRneasy,shared.by.cond$miRvana)
        abline(ls1,col="red")
        abline(0,1)
        ##POSSIBLY FILL RED IF LEAST SQUARES LINE IS ABOVE AND BLUE IF BELOW (USE TRANSPARENCY)
        
    }
    
    #For miRNA detected in all conditions
    
   
}
# print the overall labels
mtext('x-axis title', side = 1, outer = TRUE, line = 2)
mtext('y-axis title', side = 2, outer = TRUE, line = 2)
dev.off()


#Also try lattice dotplot with mirvana predicting exornaeasy conditioning on e.g. R type and/or condition 
#dotplot(Scheme ~ Area_under_ROC | Dataset, data = simulationSummary, layout = c(4,6))

####################################################################################################