#Based on the simpler 'modules as connected nodes' definition, we will look for module recurrence.
#Using the covariance structure, we keep all of the edges in the graph above some threshold for seq(25,550,25)
#  according to the edge weights between a pair of nodes.
#All of the  

source('./exrna_atlas_functions.R', echo=TRUE)

#First load 'Patel-Jensen_recurrentModules.Rdata' into the environment
#load("./input/Lilli/non_clique/recurrentModulez.RData")
#load('./input/Lilli/non_clique/recurrentModules-10comparisons.RData')
#load('./input/Lilli/non_clique/recurrentModules-15comparisons.RData')
load('./input/Lilli/non_clique/recurrentModules.fixed.10comparisons.RData')


my.dfs <- ls()[grep("summary",ls())]
my.dfs
#Prepare lists to receive outputs
out.lst.raw <- vector('list',length(my.dfs))
out.lst.shared.contrib <- vector('list',length(my.dfs))

#Name the output lists appropriately
comp.names <- gsub("summary.|.df","",my.dfs)
names(out.lst.raw) <- comp.names
names(out.lst.shared.contrib) <- comp.names


for(j in 1:length(my.dfs)){
    df <- get(my.dfs[j])
    
    #Rename the data frame columns:
    cur.graph.names <- my.dfs[j]
    cur.graph.names <- gsub("summary.","",cur.graph.names)
    cur.graph.names <- gsub(".df","",cur.graph.names)
    cur.graph.names <- unlist(strsplit(cur.graph.names,"vs") )
    colnames(df)[1] <- 'Merged.Graph.module'
    colnames(df)[2] <- 'retained.edges'
    colnames(df)[3] <- paste(cur.graph.names[1],'.contrib',sep = '')
    colnames(df)[4] <- paste(cur.graph.names[2],'.contrib',sep = '')
    #head(df) 
    
    #Split the summary frame into a list of subframes 
    #based on the number of retained edges
    df.lst <- split(df,df[,2])
     df.lst <- df.lst[naturalsort(names(df.lst))]
    
    out.dfs.lst.raw <- vector('list',length(df.lst))
    names(out.dfs.lst.raw) <- names(df.lst)
    
    out.dfs.lst.shr <- vector('list',length(df.lst))
    names(out.dfs.lst.shr) <- names(df.lst)
    for(k in 1: length(df.lst)){
        df.cur <- df.lst[[k]]
        delim <- '\\.'
        df.cur[,5]<- sapply( (strsplit(df.cur[,1],delim)),length)
        df.cur[,6] <- sapply( (strsplit(df.cur[,3],delim)),length)
        df.cur[,7] <- sapply( (strsplit(df.cur[,4],delim)),length)
        
        #Now that df has been converted to module lengths, we can
        #  get some sense of whether the counts make sense.
        df.cur[,8] <- df.cur[,5] - rowSums(df.cur[,6:7]) 
        colnames(df.cur)[5:8] <- c('Merged.Graph.module.size',
                                   paste(cur.graph.names[1],'.contrib.count',sep = ''),
                                   paste(cur.graph.names[2],'.contrib.count',sep = ''),
                                   "surplus.node.count")
        
        df.cur <- df.cur[,c(5:8,1,3:4)] #Reorder and discard the retained edges info
        out.dfs.lst.raw[[k]] <- df.cur
        
        df.shr <- df.cur[(df.cur[,2] > 0 & df.cur[,3] > 0 ),]
        out.dfs.lst.shr[[k]] <- df.shr
        
    }
    
    out.lst.raw[[j]] <- out.dfs.lst.raw
    out.lst.shared.contrib[[j]] <- out.dfs.lst.shr
    
}

#This provides a quick summary of the module overlap for two graphs
for(n in 1:length(my.dfs)){
    print(my.dfs[n])
    df <- get(my.dfs[n])
    df.T <- df[df$inGraph1 != '' & df$inGraph2 != '' ,]
    print(dim(df))
    print(dim(df.T))
    print(cat(rep('=',30),sep = "=") )
    
    
} 


shared.df <- as.data.frame(matrix(data = NA, nrow = length(out.lst.shared.contrib),
                    ncol = length(names(out.lst.shared.contrib[[1]]))))
colnames(shared.df) <- names(out.lst.shared.contrib[[1]])
rownames(shared.df) <- names(out.lst.shared.contrib)

for(cmp in 1:length(out.lst.shared.contrib)){
    shared.df[cmp,] <- sapply(out.lst.shared.contrib[[cmp]],nrow)
}


raw.df <- as.data.frame(matrix(data = NA, nrow = length(out.lst.raw),
                                  ncol = length(names(out.lst.raw[[1]]))))
colnames(raw.df) <- names(out.lst.raw[[1]])
rownames(raw.df) <- names(out.lst.raw)

for(cmp in 1:length(out.lst.raw)){
    raw.df[cmp,] <- sapply(out.lst.raw[[cmp]],nrow)
}

#write.table(x = shared.df, file = './output/shared.modules.counts10.fix.txt',quote = F, sep = '\t',
#            row.names = T,col.names = T)

#write.table(x = raw.df, file = './output/combined.modules.counts10.fix.txt',quote = F, sep = '\t',
#            row.names = T,col.names = T)



#d <- out.lst.shared.contrib$JensenCSFvsPatelPLASMA$`325`$Merged.Graph.module.size

#Show the dependence of module size and retained edges
#Determine module similarity as retained edges change

#module.sizes.summary <- function(out.lst){
    out.lst <- out.lst.shared.contrib
    ll <- vector('list',length = length(out.lst))
    names(ll) <- names(out.lst)
    tmp1 <- 1
    for(cmp in out.lst){
        df.ms <- matrix(data = NA,nrow = length(out.lst[[tmp1]]),
                        ncol = 3)
        rownames(df.ms) <- names(cmp)
        colnames(df.ms) <- c('A','B','C')
        tmp <- 1
        for(ms in cmp){
            #Put the mode value for module size for each  
            if(nrow(ms)!=0){
                #df.ms[tmp,1] <- my.get.mode(ms[,1])
                #df.ms[tmp,2] <- my.get.mode(ms[,2])
                #df.ms[tmp,3] <- ms[,3]
                df.ms[tmp,1:3] <- sapply(ms[,1:3],my.get.mode.info)
            }else{
                df.ms[tmp,1] <- 0
                df.ms[tmp,2] <- 0
                df.ms[tmp,3] <- 0
            }
            tmp <- tmp + 1
        }
        #a <-hist(ms1,breaks = seq(min(ms1),max(ms1),1))
        #c <-hist(ms1,breaks = seq(100))
        #b <-hist(ms1,breaks = seq(100))
        ll[[tmp1]] <- df.ms
        tmp1 <- tmp1 + 1
    }
m <- 11    
names(ll)[m]    
ll[[m]]
dim(out.lst[[m]][[1]])


intersect(out.lst$JensenCSFvsPatelPLASMA$`75`$Merged.Graph.module,
      out.lst$JensenCSFvsPatelPLASMA$`100`$Merged.Graph.module)




#    a$breaks[which(a$counts == max(a$counts))]
#}
