
#For the significant modules...get the largest set of connected nodes
# where any subset creates a significant module

#Later make this for only the most important subset of modules....???
###################################################################################
merge.modules <- function(mds.lst,min.ovlp = 1){
    ##THIS FUNCTION NEEDS TO BE IMPROVED:
    #it must handle cases where a module cannot be merged without discarding such a module
    lmm <- length(mds.lst)
    merged <- vector('list',lmm)
    for(i in 1:(lmm-1)){
        merg.mod <-c() #This will store the updated merged based on current module
        a <- mds.lst[[i]] 
        for(j in (i+1):lmm){
            b <- mds.lst[[j]]
            x <- intersect(a,b)
            if(length(x) >= min.ovlp){
                c <- union(unlist(a),unlist(b))
                merg.mod <- union(d,c)
            }
               
        }
        merged[[i]] <- naturalsort(merg.mod)
        #print(length(d))
    }
    return(unique(merged))
}
##################################################################################

expt.name
tmp <- rmv.sub.modules(module.members )
tmp1 = 0
tmp2 = 4
while(tmp1 != tmp2){
    tmp1 <- length(tmp)
    tmp <- rmv.sub.modules(tmp)
    tmp2 <- length(tmp)
    print(c(tmp1,tmp2))
}
tmp
tmp <- merge.modules(tmp,1)
tmp


################################################################################
mods.by.expt <- out.put.lst
identify.recurrent <- function(mods.by.expt){
    flt.lst <- unlist(mods.by.expt,recursive = F,use.names = T)
    jnk <- vector("list",length(flt.lst))
    idx <-1
    for(md1 in 1:length(flt.lst)){
        exps <- c()
        for(md2 in 1:length(flt.lst)){
            if( identical(flt.lst[[md1]],flt.lst[[md2]]) && (md1!=md2) ){
                exps <- c(exps,names(flt.lst[md2])) 
            }
        }
        jnk[[idx]] <- exps
        idx <- idx+1
    }
    jnk <- jnk[!sapply(jnk,is.null)]
    flt.lst[ jnk[47][[1]] ]
}
flt.lst["kras1.Pancreatic4"]
head(sort(unlist(jnk),decreasing = T),10)
head(sort(sapply(jnk,length),decreasing = F))
jnk$exrna2.AD.CSF81
flt.lst$exrna2.AD.CSF81

#Does the complement of modules allow for clustering of experiments
mod.mirna <- naturalsort(unique(unlist(flt.lst)))

make.modnames <- function(x){cat(unique(x),sep = '.')}
lapply(module.members,make.modnames)


module.members[[1]]
rng <- 1247:1465
cor(atlas.qn[(combn(module.members[[1]],2))[1,1],rng],atlas.qn[(combn(module.members[[1]],2))[2,1],rng])
