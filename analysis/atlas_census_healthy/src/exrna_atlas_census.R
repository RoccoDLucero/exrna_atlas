
source('./exrna_atlas_functions.R', echo=F)

head(all_mirna_rpm[,1:5])

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
##For May 16 (ISEV 2017)
load("./docs/Oscar_5_16_2017/inputVenn.txt")
str(VENN.LIST)

library(VennDiagram)
venn.title <- "15_Threshold_miRNA"    
venn.plot = venn.diagram(VENN.LIST, NULL, cex = 1.25, alpha = rep(.25,5),
                         fill=c("darkmagenta", "darkblue", "yellow", "orange", "gray"),
                         cat.fontface = 1, category.names = names(VENN.LIST),
                         main = venn.title )

pdf(paste(venn.title,"_Mean_Healthy_Biofluids_Final_2.pdf", sep = ''))

grid.draw(venn.plot)

dev.off()
