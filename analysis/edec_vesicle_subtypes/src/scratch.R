select.best.probes <- function(informative_probes, probes_per_comparison, n_classes)

    tmp1 <- lapply(c(1:n_classes),
                   function(x){
                       head(sort(informative_probes_obj$one.vs.rest$t_test$p.values[,x]), probes_ovr)})
    
    tmp2 <- unique(names(unlist(tmp1)))
    
    informative_probes1 <- tmp2
    
    tmp11 <- lapply(c(1:n_pairs),
                    function(x){head(sort(informative_probes_obj$each.pair$t_test$p.values[,x]), probes_ep)})
    
    tmp22 <- unique(names(unlist(tmp11)))
    
    informative_probes2 <- tmp22
    chosen_probes <- unlist(union(informative_probes1, informative_probes2))






##
## FUNCTION TO GET BEST PROBES BY RNA TYPE

perform.probe.select.by.probeset <- function(ref_df, ref_classes,
                                             probe_sets_for_selection,
                                             max_p_val,
                                             n_markers,
                                             method = list("one.vs.rest",
                                                           "each.pair" )){
    
        ref_sets <- lapply(X = probe_sets_for_selection,
                           FUN = function(set){ return(ref_df[set,]) })
            
        sel_probes <- lapply(X = ref_sets,
                             FUN = function(rfs){
                                 perform.probe.selection(ref_probes = rfs,
                                 ref_classes = ref_classes,
                                 max_p_val = max_p_val,
                                 n_markers = n_markers)
                             }) 
                                 
        
        return(list(probe_sets = probe_sets_for_selection,
                    selected_probes = sel_probes) )
        
        
}



"Originally Posted by alexdobin View Post
At the 3' ends of your reads you should see the reverse complementary to the A-Tail RT Primer sequence: 5’-TCTCGGCATTCCTGCTGAACCGCTCTTCCGATCTTTTTTTTTTTTVN 
i.e. 


We trimmed any sequence that contained AAAAAA (6As). This is a very aggressive trimming - we hoped to get rid of all the genomic A-homopolymer priming sites. STAR can do it for you with:
--clip3pAdapterSeq AAAAAA --clip3pAdapterMMp 0
'"

library('ShortRead')

View(gingeras_meta)

fq_dir <- "./input/gingeras_gse24565/cyt_nuc_all"
fq_files <- list.files(fq_dir,full.names = T)

fq <- fq_files[grep("SRR527587",fq_files)]
smp <- FastqSampler(con = fq, n = 10000)
yld <-yield(smp)
tbls <- tables(yld)

fq1 <- fq_files[grep("SRR527585",fq_files)]
smp1 <- FastqSampler(con = fq1, n = 10000)
yld1 <-yield(smp1)
tbls1 <- tables(yld1)
tbls1$top
################################################################################
###make this a function!!!!!!!
################################################################################
## CREATE REFERENCE SAMPLE LABELS FOR EDEC
## MAKE THIS EXPANDABLE TO MULTPLE METADATA VARIABLES IE A FUNCTION ##
if(T){
    meta_field <-gingeras_meta$sample_name
    lbls_frac <- c("cytosol","nucleus")#,"cell")
    use_lbls <- lbls_frac
    
    ref_lbl <- sapply(X = use_lbls,
                      FUN =  function(frac){
                          frac <- gsub(' ','',frac)
                          pos <- grep(frac, meta_field)
                          
                          clss <- rep(frac, length(pos))
                          
                          
                          nm <- gingeras_meta$run_accession[pos]
                          
                          
                          lbl_map <- cbind(nm,clss)
                          
                          return(lbl_map)},
                      simplify = T)
    
    ref_lbl <- Reduce(f = rbind, ref_lbl)
}

if(T){
    meta_field <- gingeras_meta$sample_attribute
    lbls_cell_ty <- unique(gsub("^source_name:|\\|\\|.*$","", meta_field))
    use_lbls <- lbls_cell_ty
    
    alt_lbl <- sapply(X = use_lbls,
                      FUN =  function(frac){
                          frac <- gsub(' ','',frac)
                          pos <- grep(frac, meta_field)
                          
                          clss <- rep(frac, length(pos))
                          
                          
                          nm <- gingeras_meta$run_accession[pos]
                          
                          
                          lbl_map <- cbind(nm,clss)
                          
                          return(lbl_map)},
                      simplify = T)
    
    alt_lbl <- Reduce(f = rbind, alt_lbl)
}

if(F){
    meta_field <- gingeras_meta$instrument_model
    lbls_instr_ty <- unique(meta_field)
    use_lbls <- lbls_instr_ty
    
    alt_lbl <- sapply(X = use_lbls,
                      FUN =  function(frac){
                          
                          pos <- grep(paste("^",frac,"$",sep = ''), meta_field)
                          
                          clss <- rep(frac, length(pos))
                          
                          
                          nm <- gingeras_meta$run_accession[pos]
                          
                          
                          lbl_map <- cbind(nm,clss)
                          
                          return(lbl_map)},
                      simplify = T)
    
    alt_lbl <- Reduce(f = rbind, alt_lbl)
}

#####################################################
## ENSURE THAT OTHER METADATA COME FROM
## SAMPLES SELECTED ON THE PRIMARY METATDATA PROPERTY
meta_reorder <- match(x = ref_lbl[,'nm'], table = alt_lbl)
alt_lbl <- alt_lbl[meta_reorder,]


# check_sra_dnlds.R

##Compare the number of bases read based on the file uploaded to SRA vs.
## the file downloaded from SRA

##ADD THIS TO THE SRA GET RAW DATA SCRIPT.

## LOAD METADATA from SRA_DB SRA table
## LOAD EXCERPT OUTPUT "READ MAPPING SUMMARY"

cel_nuc_mapp_summary <- read_delim("P:/brl/proj/EXRNA_~1/analysis/EDEC_V~1/input/GINGER~1/EXCERP~1/GINGER~1/RDLUCE~1/EXRNAA~1/EXCERP~1.2/GINGER~1/POSTPR~1.3/GI93E2~1.TXT",
                                   "\t", escape_double = FALSE, trim_ws = TRUE)


cyt_mapp_summary <- read_delim("P:/brl/proj/EXRNA_~1/analysis/EDEC_V~1/input/GINGER~1/EXCERP~1/EXCERP~1/RDLUCE~1/EXRNAA~1/EXCERP~1.2/GINGER~1/POSTPR~1.3/GI1D14~1.TXT",
                               "\t", escape_double = FALSE, trim_ws = TRUE)

ref_path  <- "../edec_vesicle_subtypes/input/lasser_atlas_gingeras/ref/"
in_file <- "gingeras_subcell_meta.RDS"
gingeras_meta <- readRDS(paste(ref_path, in_file,sep = ''))


combd_summary <- rbind(cel_nuc_mapp_summary,cyt_mapp_summary)
combd_summary$X1 <- gsub('sample_|_fastq','', combd_summary$X1)

xx <- match(x = combd_summary$X1, table = gingeras_meta$run_accession)

yy <- gingeras_meta[xx,]
View(yy)

yy$run_accession == combd_summary$X1

zz <- cbind(yy$spots,combd_summary$input )
zzz <- yy$spots == combd_summary$input
all(zzz)


x <- my.run.edec.stage.0(version = 'one.vs.rest', reference_meth = probe_sel_refs,
                    reference_classes = ref_classes,
                    max_p_value = probe_sel_max_p, num_markers = n_mrkrs)

