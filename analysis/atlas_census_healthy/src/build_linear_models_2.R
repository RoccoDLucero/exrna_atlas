

################################################################################
#source("./exrna_atlas_functions.R")
library(tsne)
library(nnet)
library(igraph)
library(ggplot2)
library(ggfortify)
library(car)

#Define useful functions:
my.get.meta <-function(df,ptt = "hsa-"){return(df[,-grep(ptt,colnames(df)),])}

################################################################################
my.make.models <- function(cluster.set, atlas.data,predictors = NULL,
                           age = T, sex = T, kit = T){
    fit.age <- NULL
    fit.kit <- NULL
    fit.sex <-NULL
    preds <- predictors
    options(warn = -1)
    atlas.data <- droplevels(atlas.data)
    #cluster.set <- my.cluster.sets
    modeling.results <- vector(mode = 'list', length(cluster.set))
    i <- 0
    for(mds in cluster.set){
        i <- i+1
        module.meta.pred <- vector(mode = 'list', length(mds))
        j <- 0
        for(mdl in mds){
            j <- j+1
            fit.models.aic <- list("sex" = NULL, "age" = NULL, "kit" = NULL, "module" = NULL) 
            md <-  names(V(mdl)) #Get the list of miRNA in current module
            my.data.sex <- atlas.data[,c(md,"Donor.Sex", preds)]
            if(sex){fit.sex <- glm(Donor.Sex ~. ,family=binomial(link='logit'),data = my.data.sex)}
            
            my.data.age <- atlas.data[,c(md,"Donor.Age", preds)]
            if(age){fit.age <- lm(Donor.Age ~. ,data = my.data.age ) 
                s <-summary(fit.age)
            }
     
            my.data.kit <- atlas.data[,c(md,"RNA_isolation_kit",preds)]
            
            if(kit){fit.kit <- multinom(RNA_isolation_kit ~. , data=my.data.kit) }
            
            fit.models.aic$sex <- list('AIC' = fit.sex$aic)
            fit.models.aic$age <- list('AIC' = AIC(fit.age), 'RSQ' = s$r.squared)
            fit.models.aic$kit <- list('AIC' = fit.kit$AIC)
            fit.models.aic$module <- md
            
            module.meta.pred[[j]] <- fit.models.aic 
        }
        
        modeling.results[[i]] <- module.meta.pred
        names(modeling.results)[i] <- names(cluster.set)[i]
    }
    return(modeling.results)
}    
################################################################################
my.compile.results <- function(model.res = modeling.results){
    #Turn the modeling results into a more useful form
    #Should be a list of vectors for, e.g. AIC, from each response metadata variable
    
    #Define Sub-routines
    return.modules <- function(t){
        t[grep("module",names(t))]
    }
    
    #For each threshold, for each module get data:
    tmp <- unlist(model.res)
    sex.AIC <- tmp[grep("sex",names(tmp))]
    age.rsq <- tmp[grep("age.RSQ",names(tmp))]
    age.AIC <- tmp[grep("age.AIC",names(tmp))]
    kit.AIC <- tmp[grep("kit",names(tmp))]
    
    modules <- lapply(model.res,unlist,recursive = F)
    modules <- unlist(sapply(modules,return.modules),recursive = F)
    
    
    names(sex.AIC) <- make.names(names(sex.AIC),unique = T)
    names(age.AIC) <- make.names(names(age.AIC),unique = T)
    names(age.rsq) <- make.names(names(age.rsq),unique = T)
    names(kit.AIC) <- make.names(names(kit.AIC),unique = T)
    names(modules) <- make.names(names(modules),unique = T)
    
    tmp.stats <- cbind(sex.AIC, age.AIC, age.rsq, kit.AIC)
    tmp.stats <- as.data.frame(apply(tmp.stats,2,as.numeric))
    rownames(tmp.stats) <- gsub(".sex.AIC","",names(sex.AIC))
    names(modules) <- gsub(".module","",names(modules))
    #all(rownames(tmp.stats) == names(modules))
    
    return(list('model.scores' = tmp.stats, 'modules' = modules))
    
}
################################################################################
my.make.comp.plots <- function(mod.score.lst, sources.list, show.plots = T, return.combined = F){
    dat <-NULL
    mods <- vector("list",sum(sapply(mod.score.lst,function(x){length(x$modules)})))
    m <- 0
    for(l in 1:length(mod.score.lst)){
        dat.source <- sources.list[[l]]
        scores <- mod.score.lst[[l]]
        md.size <- sapply(scores$modules,length)
        sc <- cbind(scores$model.scores, md.size) #,md.content)
        sc[,'source'] <- dat.source
        dat <- rbind(dat,sc)
        for(i in 1:length(scores$modules)){
            m <- m+1
            #print(m)
            mods[m] <- scores$modules[i]
        }
        
    }
    
    #Plotting ##SEPARATE THIAS FUNCTIONALITY### #
    p1 <- qplot(dat$age.rsq, dat$age.AIC,main = "AIC vs. R-squared for linear model of age",
                xlab = "R-squared", ylab = "AIC", color = dat$source)
    p2 <- qplot(dat$md.size, dat$age.AIC,
                main = "AIC vs. module size for linear model of age",
                xlab = "module size", ylab = "AIC", color = dat$source)
    p3 <- qplot(dat$md.size, dat$age.rsq,
                main = "R-squared vs. module size for linear model of age",
                xlab = "module size", ylab = "R-squared",color = dat$source)
    p4 <- qplot(dat$md.size, dat$sex.AIC,
                main = "AIC vs. module size for logistic model of gender",
                xlab = "module size", ylab = "AIC",color = dat$source)
    p5 <-qplot(dat$md.size, dat$kit.AIC,
               main = "AIC vs. module size for multinomial model of kit",
               xlab = "module size", ylab = "AIC",color = dat$source)
    
    plts <- list(p1,p2,p3,p4,p5)
    if(show.plots){sapply(plts,print)}
    
    
    
    if(return.combined){return(list(data=dat, modules = mods, plots=plts))}
}
################################################################################
################################################################################
################################################################################
#ad.z comes from metadataprocessing.R 
#Subset only healthy controls
ad.z.healthy <- ad.z[ad.z$condition == "Healthy Control",]
#Remove problematic records:
my.data <- ad.z.healthy[grep("Male|Female",ad.z.healthy$Donor.Sex),]
my.data <- my.data[which(my.data$Donor.Age>=10),]
#Some possible transforms on the data
#npn.dat <- huge.npn(  as.matrix( my.data[,grep("hsa-",colnames(my.data))] ) )
qn.dat <- normalize.quantiles(  as.matrix( my.data[,grep("hsa-",colnames(my.data))] ) )
qn.npn.dat <- huge.npn(qn.dat)
my.data[,grep("hsa-",colnames(my.data))] <- qn.npn.dat   

#my.cluster.sets comes from AtlasSecondOrderModuleAnalysis.R
##It is a list of lists of connected components generated at several edgeweight thresholds
my.cluster.sets <- readRDS("./output/test_clusters_30_thresholds")
my.cluster.sets.cmp1 <- readRDS("./output/test_clusters_cmplmt_30_thresholds")
my.cluster.sets.cmp2 <- readRDS("./output/test_clusters_cmplmt_30_thresholds_2")

#Model metadata with the miRNA expression data
#my.data <- my.data[my.data$biofluid_name == "Plasma" ,]
modeling.results <- my.make.models(cluster.set = my.cluster.sets, atlas.data = my.data)
modeling.results.cmp1 <- my.make.models(cluster.set = my.cluster.sets.cmp1, atlas.data = my.data)
modeling.results.cmp2 <- my.make.models(cluster.set = my.cluster.sets.cmp2, atlas.data = my.data)
    
modeling.scores <- my.compile.results(modeling.results)
modeling.scores.cmp1 <- my.compile.results(modeling.results.cmp1)
modeling.scores.cmp2 <- my.compile.results(modeling.results.cmp2)

#Plot overall modeling performance & select most predictive modules
m.sc.lst <- list(modeling.scores,modeling.scores.cmp1,modeling.scores.cmp2)
src.lst <- list("101@15rpm","101random-a","101random-b")
comp.plots <- my.make.comp.plots(m.sc.lst, src.lst,show.plots = T,return.combined = T)

comp.dat <- comp.plots$data[!duplicated(comp.plots$modules),]
comp.mods <- comp.plots$modules[!duplicated(comp.plots$modules)]

sex.modules <- list(dat = comp.dat[order(comp.dat$sex.AIC,decreasing = F),colnames(comp.dat) %in% c("sex.AIC","md.size","source")],
                    mods = comp.mods[order(comp.dat$sex.AIC,decreasing = F)])
age.modules <- list(dat = comp.dat[order(comp.dat$age.AIC,decreasing = F),colnames(comp.dat) %in% c("age.AIC","age.rsq","md.size","source")],
                    mods = comp.mods[order(comp.dat$age.AIC,decreasing = F)])
kit.modules <- list(dat = comp.dat[order(comp.dat$kit.AIC,decreasing = F),colnames(comp.dat) %in% c("kit.AIC","md.size","source")],
                    mods = comp.mods[order(comp.dat$kit.AIC,decreasing = F)])

g <- age.modules$mods[[1]]
prp = 50
ts <- tsne(X = my.data[,g],k = 5,perplexity = prp)
qplot(ts[,3],ts[,2],color = my.data$Donor.Age, shape = my.data$Donor.Sex)
pca <- prcomp(x = my.data[,g])
qplot(pca$x[,3],pca$x[,2],color = my.data$Donor.Age, shape = my.data$Donor.Sex)
ts2 <- tsne(X = my.data[,sample(colnames(my.data[,2:1600]),length(g),replace = F)],k = 5,perplexity = prp)
qplot(ts2[,1],ts2[,2],color = my.data$RNA_isolation_kit, shape = my.data$Donor.Sex)

cbind(age.modules$dat,seq(1:nrow(age.modules$dat) ))
preds <- c(NULL)
g <- age.modules$mods[[27]]  #c(24,27,37,40,42,44,48:54)
my.data.age <- my.data[,c(g,"Donor.Age",preds)]
my.data.age <- my.data[,c(zzz.intersect,"Donor.Age" )]
fit.age <- lm(Donor.Age ~. ,data = my.data.age ) 
#fit.age <- glm(Donor.Age ~. ,data = my.data.age,family = poisson() ) 
summary(fit.age)
#plot(fit.age)
av<- avPlots(model = fit.age,intercept = T)

 g <- sex.modules$mods[[1]]
my.data.sex <- my.data[,c(g,"Donor.Sex")]
fit.sex <- glm(Donor.Sex ~. ,family=binomial(link='logit'),data = my.data.sex)
summary(fit.sex)
avPlots(model = fit.sex,intercept = T)

par(mfrow = c(3,3))
apply(z.healthy[,3:40],2,hist, breaks = 30)
x <-  apply(z.healthy[,3:40],2,log10)
x <-  normalize.quantiles(as.matrix(z.healthy[,3:40]))
apply(x,2,hist, breaks = 30)

