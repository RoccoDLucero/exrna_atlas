
source('./exrna_atlas_functions.R', echo=F)
load('./input/Lilli/edge_only/7graphs_forRocco.RData')

################################################################################
##This section should be used to replicate the "first order" results Joel
##presented in the Decembe 8 RNA-seq call. We want to know the degree of
##unique and shared miRNA for various sets of bio fluids represented in the Atlas.

##Select the miRNA from a given experiment ('PI.biofluid."data"') which occur
## with a mean value above a given threshold
my.thresh.mirna <- function(df,thresh){x <- rowMeans(df,na.rm = T);x[x>thresh] }

serum.data <- cbind(jensen.serum.data,laurent.serum.data) # Combine Serum experiments to approximated Joel's analysis
mirna.dfs <- ls()[grep(pattern = 'data',ls())]
mirna.dfs <- sapply(mirna.dfs,get)
mirna.dfs <- mirna.dfs[c(1,4:6)] #Select only the 4 dataframes of interest

#Set mean readcound threshold to 100 reads for each data set 
t100.mirna.counts <- lapply(mirna.dfs,my.thresh.mirna,thresh = 100)
#Get the names of the miRNA in each experiment so we can compare in a VENN DIAGRAM
t100.mirna.names <- sapply(t100.mirna.counts,names)

#Plot the Venn Diagram and also get the species present in each overlap
venn.diagram(t100.mirna.names,filename = './output/test.venn.png',
             imagetype = 'png',height = 2000,width = 5000)

cn <- calculate.overlap(t100.mirna.names) #This contains the miRNA species for each section of the Venn Diagram

################################################################################

################################################################################
##This section will be for our 'second order' analysis, which describes the 
## most highly co-varying pairs of miRNA shared by sets of biofluids
##Lilli proveded the igraph structures

E(ig.jensen.serum)[(E(ig.jensen.serum) %in% E(ig.jensen.serum))]
E(ig.jensen.serum)[(E(ig.jensen.serum) %in% E(ig.jensen.csf))]

E(ig.jensen.serum)[(E(ig.jensen.serum) %in% E(ig.laurent.serum))]
E(ig.jensen.serum)[(E(ig.jensen.serum) %in% E(ig.patel.plasma))]

#Wrie a function to remove single-node edges
edg <- E(ig.jensen.csf)
as_ids(edg[1])
ea <- unlist(edge.attributes(ig.jensen.csf)) #command to get edge weights
str(ea)
ea[ea>0] #select all edges with positive weights
head(ea[order(ea,decreasing = T)])

as_ids(edg[order(ea,decreasing = T)][1:100],)
(edg[order(ea,decreasing = T)][1:100])

