source('./exrna_atlas_functions.R',echo = F)

options(scipen = 5)
options(digits=4)

All.exps.output <- readRDS('./CellAtlasMod_OutLst10.26.rdat')
All.names <- names(All.exps.output)
expOutData <- names(All.exps.output$exrna1.AD)

query.data <- exRna.Atlas.miR
#query.data <- miRNA.RPM
colors.meta <-Atlas.miR.meta$Condition
shape.meta <- Atlas.miR.meta$Biofluid.Name
#miRNA.RPM <- rbind.data.frame(Tissue.Atlas,SRA.miR,exRna.Atlas.miR)


#Loop throug all experiments eventually
d.do.expts <- grep('Panc',names(All.exps.output),value = T)
for(cur.expt in d.do.expts[1:2] ){
    #for the current experiment get the number of significant modules 
    n.sig.mods<-length(All.exps.output[[cur.expt]][["mods.and.infos"]])
    n.top.mods <- 5
    
    module.members.lst <- list()
    module.info.vec <- c()
    for(j in 1:n.sig.mods){
        modMembs <- unlist(All.exps.output[[cur.expt]][["mods.and.infos"]][[j]][[1]])
        module.members.lst[[j]] <- modMembs
        
        modInfo <- unlist(All.exps.output[[cur.expt]][["mods.and.infos"]][[j]][[2]])
        module.info.vec[j] <- modInfo
    }
    
    #pplot.file <- paste('./output/Cell_Atlas/',cur.expt,'.permodule.cors.hist.pdf',sep = '')
    #pdf(file = pplot.file, onefile = T)
    
    top.mod.im.prs <- list()
    for(m in 1:n.top.mods){
        in.mod <- module.members.lst[[m]]
        if(length(in.mod)<2){next}
       
        all.ex.nodes <- unlist(module.members.lst)
        all.ex.nodes <- unique(all.ex.nodes)
        all.ex.nodes <- setdiff(all.ex.nodes,in.mod)
        
        #Enumerate all intramodule pairs
        im.prs <- NULL
        if(length(in.mod)>1){im.prs <- combn(in.mod,2,simplify = F)
        }else{print('*')}
        #Enumerate all
        
        imex.prs <- list() 
        cnt <-0
        for(in.mod.nd in in.mod){
            for(ex.mod.nd in all.ex.nodes){
                cnt <- cnt+1
                pr <- c(in.mod.nd,ex.mod.nd)
                imex.prs[[cnt]] <- pr
            }
        }
        
        all.ex.prs <- NULL
        #Enumerate all purely extramodule pairs
        all.ex.prs <-combn(all.ex.nodes,2,simplify = F)
        
        im.cors <- c()
        for(aa in 1:length(im.prs)){
            if(!(im.prs[[aa]][1] %in% colnames(query.data))){next}
            if(!(im.prs[[aa]][2] %in% colnames(query.data))){next}   
            ao <- cor(query.data[,im.prs[[aa]][1]],
                      query.data[,im.prs[[aa]][2]])
            im.cors <- c(im.cors,ao)
        }
        
        imex.cors <- c()
        for(aa in 1:length(imex.prs)){
            if(!(imex.prs[[aa]][1] %in% colnames(query.data))){next}
            if(!(imex.prs[[aa]][2] %in% colnames(query.data))){next}
            ao <- cor(query.data[,imex.prs[[aa]][1]],
                      query.data[,imex.prs[[aa]][2]])
            imex.cors <- c(imex.cors,ao)
        }
        
        all.ex.cors <- c()
        for(aa in 1:length(all.ex.prs)){
            if(!(all.ex.prs[[aa]][1] %in% colnames(query.data))){next}
            if(!(all.ex.prs[[aa]][2] %in% colnames(query.data))){next}
            ao <- cor(query.data[,all.ex.prs[[aa]][1]],
                      query.data[,all.ex.prs[[aa]][2]])
            all.ex.cors <- c(all.ex.cors,ao)
        }
        if(F){
        txp <- .3
        hist(imex.cors,breaks = 20, col=rgb(1,0,0,txp), freq = F,xlim = c(-.5,1))
        hist(im.cors,add =T,breaks = 20, col=rgb(0,1,0,txp), alpha =.5, freq = F )
        hist(all.ex.cors,add =T, col=rgb(0,0,1,txp) ,freq = F)
        }
###################################################################################################
        ###################################################################################################                
        pca.data <-(query.data[,in.mod])
        
        non.mod.nds <- colnames(query.data)[!(colnames(query.data) %in% unique(unlist(module.members.lst)))]
        pca.data.rnd <-(query.data[sample(non.mod.nds,length(in.mod))])
    
        #pca.data <-(query.data[,sample(colnames(query.data),length(in.mod))])
        #pca.data <-(query.data[,sample(all.ex.nodes,length(in.mod))])

        ##Look at overview of PCA results
        #res.pca <- prcomp(my.cell.atlas[,-2141])
        
        ########################################################################
        res.pca <- prcomp(pca.data)
        eig <- (res.pca$sdev)^2
        component.variance <- eig/sum(eig)*100
        cumvar <- cumsum(component.variance)
        
        pdf(paste('./output/Cell_Atlas/',cur.expt,'-mod',m,'.PCA.pdf',sep = ''),
            onefile = T,width = 16, height = 9)
        if(T){
        PC4scree <- 12 #number of principla components to plot
        barplot(component.variance[1:PC4scree], names.arg=1:length(component.variance[1:PC4scree]), 
                main = paste(cur.expt,'true module',m,'\n% Variance of PCs 1-12'),
                sub = "% variance explained",
                xlab = "Principal Components",
                ylab = "Percentage of variances",
                col ="skyblue")
        }

        var_cor_func <- function(var.loadings, comp.sdev){
            var.loadings*comp.sdev
        }
        # Variable correlation/coordinates
        loadings <- res.pca$rotation
        sdev <- res.pca$sdev
        var.coord <- var.cor <- t(apply(loadings, 1, var_cor_func, sdev))
        
        tmp1 <- predict(res.pca)
        tmp1 <- as.data.frame(tmp1)

        print(qplot(PC1, PC2, data=tmp1, colour=colors.meta,asp = 9/16, shape = shape.meta,
                    main = paste('Top module',m, 'for', cur.expt),
              xlab = paste('PC1 (', round(component.variance[1],1),'% variance)'),
              ylab = paste('PC2 (', round(component.variance[2],1),'% variance)') ) )
        
        if(ncol(tmp1)>2){
        print( qplot(PC2, PC3, data=tmp1, colour=colors.meta,asp = 9/16, shape = shape.meta,
                     main = paste('Top module',m, 'for', cur.expt),
                     xlab = paste('PC2 (', round(component.variance[2],1),'% variance)'),
                     ylab = paste('PC3 (', round(component.variance[3],1),'% variance)') ) )
        
        print( qplot(PC1, PC3, data=tmp1, colour=colors.meta,asp = 9/16, shape = shape.meta,
                     main = paste('Top module',m, 'for', cur.expt),
                     xlab = paste('PC1 (', round(component.variance[1],1),'% variance)'),
                     ylab = paste('PC3 (', round(component.variance[3],1),'% variance)') ) )
        }
        
        
####################        ####################        ####################        ####################
    
        res.pca <- prcomp(pca.data.rnd)
        eig <- (res.pca$sdev)^2
        component.variance <- eig/sum(eig)*100
        cumvar <- cumsum(component.variance)
        
        if(T){
        PC4scree <- 12 #number of principla components to plot
        barplot(component.variance[1:PC4scree], names.arg=1:length(component.variance[1:PC4scree]), 
                main = paste(cur.expt,'random module',m,'\n% Variance of PCs 1-12'),
                sub = "% variance explained",
                xlab = "Principal Components",
                ylab = "Percentage of variances",
                col ="darkblue")
        }
        
        var_cor_func <- function(var.loadings, comp.sdev){
            var.loadings*comp.sdev
        }
        # Variable correlation/coordinates
        loadings <- res.pca$rotation
        sdev <- res.pca$sdev
        var.coord <- var.cor <- t(apply(loadings, 1, var_cor_func, sdev))
      
        tmp1 <- predict(res.pca)
        tmp1 <- as.data.frame(tmp1)
        
        print(qplot(PC1, PC2, data=tmp1, colour=colors.meta,asp = 9/16, shape = shape.meta,
                    main = paste('"Random" module',m, 'for', cur.expt),
                    xlab = paste('PC1 (', round(component.variance[1],1),'% variance)'),
                    ylab = paste('PC2 (', round(component.variance[2],1),'% variance)') ) )
        
        if(ncol(tmp1)>2){
            print( qplot(PC2, PC3, data=tmp1, colour=colors.meta,asp = 9/16, shape = shape.meta,
                         main = paste('"Random" module',m, 'for', cur.expt),
                         xlab = paste('PC2 (', round(component.variance[2],1),'% variance)'),
                         ylab = paste('PC3 (', round(component.variance[3],1),'% variance)') ) )
            
            print( qplot(PC1, PC3, data=tmp1, colour=colors.meta,asp = 9/16, shape = shape.meta,
                         main = paste('"Random" module',m, 'for', cur.expt),
                         xlab = paste('PC1 (', round(component.variance[1],1),'% variance)'),
                         ylab = paste('PC3 (', round(component.variance[3],1),'% variance)') ) )
        }
        dev.off()
        
    }
    #dev.off()
}