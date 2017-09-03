{
    "collab_server" : "",
    "contents" : "###############################\n## EDec_OscarMurillo\n###############################\nrm(list=ls())\nlibrary(gplots)\nlibrary(class)\nlibrary(RColorBrewer)\nlibrary(stringr)\nlibrary(ggplot2)\nlibrary(reshape2)\nlibrary(survival)\nlibrary(rms)\nlibrary(devtools)\nlibrary(clue)\n\nquantile_normalisation <- function(df){\n  df_rank <- apply(df,2,rank,ties.method=\"min\")\n  df_sorted <- data.frame(apply(df, 2, sort))\n  df_mean <- apply(df_sorted, 1, mean)\n  \n  index_to_mean <- function(my_index, my_mean){\n    return(my_mean[my_index])\n  }\n  \n  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)\n  rownames(df_final) <- rownames(df)\n  return(df_final)\n}\n\n#setwd(\"~/Documents/BCM/Lab/Glioblastoma/EDec_Analysis_12\")\n#source(\"~/Documents/BCM/Lab/Glioblastoma/sources/ComBat.R\")\n#source(\"~/Documents/BCM/Lab/Glioblastoma/sources/PAM50ExerciseFunctions.R\")\n#source(\"~/Documents/BCM/Lab/Glioblastoma/sources/EDec.R\")\n#source(\"~/Documents/BCM/Lab/Glioblastoma/sources/K_cell_number.R\")\n#source(\"~/Documents/BCM/Lab/Glioblastoma/sources/diffExpEDec.R\")\n\ndir.create(\"./Glioblastoma/EDec_Analysis_12/Step1_celltype_results\", showWarnings = TRUE, recursive = FALSE)\ndir.create(\"./Documents/BCM/Lab/Glioblastoma/EDec_Analysis_12/Step3_TCGA_refSamples_results\", showWarnings = TRUE, recursive = FALSE)\n\n############################\n## Read in TCGA data\n## Read in subset beta table\n## Read in subset meta table\n## Read in combat files\n############################\nTCGA.original = read.table(\"~/Documents/BCM/Lab/Glioblastoma/EDec_Analysis_12/input_EDec_Analysis_12/TCGA_Analysis_12.txt\", header=TRUE, check.names = FALSE, sep = \"\\t\", row.names = NULL)\nTCGA.table = read.table(\"~/Documents/BCM/Lab/Glioblastoma/EDec_Analysis_12/input_EDec_Analysis_12/adjusted_TCGA_Analysis_12.txt\", header=TRUE, check.names = FALSE, sep = \"\\t\", row.names = 1)\nTCGA.table.combat = read.table(\"~/Documents/BCM/Lab/Glioblastoma/EDec_Analysis_12/input_EDec_Analysis_12/combat_TCGA_Analysis_12.txt\", header = TRUE, sep = \"\\t\", row.names = NULL)\nTCGA.table.clinic = read.table(\"~/Documents/BCM/Lab/Glioblastoma/EDec_Analysis_12/input_EDec_Analysis_12/clinical_TCGA_Analysis_12.txt\", header = TRUE, sep = \"\\t\", row.names = NULL)\n\nbeta.table = read.table(\"~/Documents/BCM/Lab/Glioblastoma/EDec_Analysis_12/input_EDec_Analysis_12/beta_subsamples_Analysis_12.txt\", header = TRUE, sep = \"\\t\", row.names = 1)\nmeta.table = read.table(\"~/Documents/BCM/Lab/Glioblastoma/EDec_Analysis_12/input_EDec_Analysis_12/meta_subsamples_Analysis_12.txt\", header = TRUE, sep = \"\\t\", row.names = NULL)\nmeta.table.combat = read.table(\"~/Documents/BCM/Lab/Glioblastoma/EDec_Analysis_12/input_EDec_Analysis_12/combat_subsamples_Analysis_12.txt\", header = TRUE, sep = \"\\t\", row.names = NULL)\noverlapProbes = as.character(intersect(row.names(TCGA.table), row.names(beta.table)))\nbeta.table = beta.table[overlapProbes,]\nTCGA.table = TCGA.table[overlapProbes,]\n\n##################################################\n## Subset data to only use XXXXXXXXXXXX subtypes\n##################################################\n## Column 7 contains cell subtype information\nmeta.table.studies = as.character(meta.table[,1])\nmeta.table.subtypes = as.character(meta.table[,7])\nnames(beta.table) = meta.table.subtypes\n\n## Use table(meta.table.subtypes) to print subtypes\nmeta.table.colors = rep(\"white\",length(meta.table.subtypes))\nmeta.table.colors[meta.table.subtypes==\"endothelial_GBM\"] = \"red\"\nmeta.table.colors[meta.table.subtypes==\"glia\"] = \"aquamarine\"\nmeta.table.colors[meta.table.subtypes==\"glioblastoma\"] = \"green\"\nmeta.table.colors[meta.table.subtypes==\"immune\"] = \"navy\"\nmeta.table.colors[meta.table.subtypes==\"neural_progenitor\"] = \"yellow\"\nmeta.table.colors[meta.table.subtypes==\"neuron\"] = \"darkorange\"\n\nTCGA.table.clinic.subtype = as.character(TCGA.table.clinic$GeneExp_Subtype)\nTCGA.table.clinic.colors = rep(\"white\",length(TCGA.table.clinic.subtype))\nTCGA.tablfo[,c(1:2,covariates)]}
	design <- design.mat(saminfo)	


	batches <- list.batch(saminfo)
	n.batch <- length(batches)
	n.batches <- sapply(batches, length)
	n.array <- sum(n.batches)
	
	## Check for missing values
	NAs = any(is.na(dat))
	if(NAs){cat(c('Found',sum(is.na(dat)),'Missing Data Values\n'),sep=' ')}
        #print(dat[1:2,])
	##Standardize Data across genes
	cat('Standardizing Data across genes\n')
	if (!NAs){B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat))}else{B.hat=apply(dat,1,Beta.NA,design)} #Standarization Model
	grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
	if (!NAs){var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)}else{var.pooled <- apply(dat-t(design%*%B.hat),1,var,na.rm=T)}

	stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
	if(!is.null(design)){tmp <- design;tmp[,c(1:n.batch)] <- 0;stand.mean <- stand.mean+t(tmp%*%B.hat)}	
	s.data <- (dat-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array)))

	##Get regression batch effect parameters
	cat("Fitting L/S model and finding priors\n")
	batch.design <- design[,1:n.batch]
	if (!NAs){gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))}else{gamma.hat=apply(s.data,1,Beta.NA,batch.design)}
	delta.hat <- NULL
	for (i in batches){
		delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=T))
		}

	##Find Priors
	gamma.bar <- apply(gamma.hat, 1, mean)
	t2 <- apply(gamma.hat, 1, var)
	a.prior <- apply(delta.hat, 1, aprior)
	b.prior <- apply(delta.hat, 1, bprior)

	
	##Plot empirical and parametric priors

	if (prior.plots & par.prior){
		par(mfrow=c(2,2))
		tmp <- density(gamma.hat[1,])
		plot(tmp,  type='l', main="Density Plot")
		xx <- seq(min(tmp$x), max(tmp$x), length=100)
		lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
		qqnorm(gamma.hat[1,])	
		qqline(gamma.hat[1,], col=2)	
	
		tmp <- density(delta.hat[1,])
		invgam <- 1/rgamma(ncol(delta.hat),a.prior[1],b.prior[1])
		tmp1 <- density(invgam)
		plot(tmp,  typ='l', main="Density Plot", ylim=c(0,max(tmp$y,tmp1$y)))
		lines(tmp1, col=2)
		qqplot(delta.hat[1,], invgam, xlab="Sample Quantiles", ylab='Theoretical Quantiles')	
		lines(c(0,max(invgam)),c(0,max(invgam)),col=2)	
		title('Q-Q Plot')
	}
	
	##Find EB batch adjustments

	gamma.star <- delta.star <- NULL
	if(par.prior){
		cat("Finding parametric adjustments\n")
		for (i in 1:n.batch){
			temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
			gamma.star <- rbind(gamma.star,temp[1,])
			delta.star <- rbind(delta.star,temp[2,])
			}
	}else{
		cat("Finding nonparametric adjustments\n")
		for (i in 1:n.batch){
			temp <- int.eprior(as.matrix(s.data[,batches[[i]]]),gamma.hat[i,],delta.hat[i,])
			gamma.star <- rbind(gamma.star,temp[1,])
			delta.star <- rbind(delta.star,temp[2,])
			}
		}


	### Normalize the Data ###
	cat("Adjusting the Data\n")

	bayesdata <- s.data
	j <- 1
	for (i in batches){
		bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
		j <- j+1
		}

	bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean
	if(write){
		output_file <- paste('Adjusted',expression_xls,'.xls',sep='_')
                 #print(geneinfo[1:2])
                 #print(bayesdata[1:2,1:4])
		 #cat(c(colnames(geneinfo),colnames(dat),'\n'),file=output_file,sep='\t')
		#suppressWarnings(write.table(cbind(geneinfo,formatC(as.matrix(bayesdata), format = "f")), file=output_file, sep="\t", quote=F,row.names=F,col.names=F,append=T))
                outdata <- cbind(ProbeID=geneinfo, bayesdata); write.table(outdata, file=output_file, sep="\t")
		cat("Adjusted data saved in file:",output_file,"\n")
		}else{return(cbind(geneinfo,bayesdata))}
	}

# filters data based on presence/absence call
filter.absent <- function(x,pct){
	present <- T
	col <- length(x)/2
	pct.absent <- (sum(x[2*(1:col)]=="A") + sum(x[2*(1:col)]=="M"))/col
	if(pct.absent > pct){present <- F}
	present
	}

# Next two functions make the design matrix (X) from the sample info file 
build.design <- function(vec, des=NULL, start=2){
	tmp <- matrix(0,length(vec),nlevels(vec)-start+1)
	for (i in 1:ncol(tmp)){tmp[,i] <- vec==levels(vec)[i+start-1]}
	cbind(des,tmp)
	}

design.mat <- function(saminfo){
	tmp <- which(colnames(saminfo) == 'Batch')
	tmp1 <- as.factor(saminfo[,tmp])
	cat("Found",nlevels(tmp1),'batches\n')
	design <- build.design(tmp1,start=1)
	ncov <- ncol(as.matrix(saminfo[,-c(1:2,tmp)]))
	cat("Found",ncov,'covariate(s)\n')
	if(ncov>0){
		for (j in 1:ncov){
			tmp1 <- as.factor(as.matrix(saminfo[,-c(1:2,tmp)])[,j])
			design <- build.design(tmp1,des=design)
			}
		}
	design
	}

# Makes a list with elements pointing to which array belongs to which batch
list.batch <- function(saminfo){
	tmp1 <- as.factor(saminfo[,which(colnames(saminfo) == 'Batch')])
	batches <- NULL
	for (i in 1:nlevels(tmp1)){batches <- append(batches, list((1:length(tmp1))[tmp1==levels(tmp1)[i]]))}
	batches
	}

# Trims the data of extra columns, note your array names cannot be named 'X' or start with 'X.'
trim.dat <- function(dat){
	tmp <- strsplit(colnames(dat),'\\.')
	tr <- NULL
	for (i in 1:length(tmp)){tr <- c(tr,tmp[[i]][1]!='X')}
	tr
	}

# Following four find empirical hyper-prior values
aprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (2*s2+m^2)/s2}
bprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (m*s2+m^3)/s2}
postmean <- function(g.hat,g.bar,n,d.star,t2){(t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)}
postvar <- function(sum2,n,a,b){(.5*sum2+b)/(n/2+a-1)}


# Pass in entire data set, the design matrix for the entire data, the batch means, the batch variances, priors (m, t2, a, b), columns of the data  matrix for the batch. Uses the EM to find the parametric batch adjustments

it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
	n <- apply(!is.na(sdat),1,sum)
	g.old <- g.hat
	d.old <- d.hat
	change <- 1
	count <- 0
	while(change>conv){
		g.new <- postmean(g.hat,g.bar,n,d.old,t2)
		sum2 <- apply((sdat-g.new%*%t(rep(1,ncol(sdat))))^2, 1, sum,na.rm=T)
		d.new <- postvar(sum2,n,a,b)
		change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old)
		g.old <- g.new
		d.old <- d.new
		count <- count+1
		}
	#cat("This batch took", count, "iterations until convergence\n")
	adjust <- rbind(g.new, d.new)
	rownames(adjust) <- c("g.star","d.star")
	adjust
	}

#likelihood function used below
L <- function(x,g.hat,d.hat){prod(dnorm(x,g.hat,sqrt(d.hat)))}

# Monte Carlo integration function to find the nonparametric adjustments
int.eprior <- function(sdat,g.hat,d.hat){
	g.star <- d.star <- NULL
	r <- nrow(sdat)
	for(i in 1:r){
		g <- g.hat[-i]
		d <- d.hat[-i]		
		x <- sdat[i,!is.na(sdat[i,])]
		n <- length(x)
		j <- numeric(n)+1
		dat <- matrix(as.numeric(x),length(g),n,byrow=T)
		resid2 <- (dat-g)^2
		sum2 <- resid2%*%j
		LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
		LH[LH=="NaN"]=0
		g.star <- c(g.star,sum(g*LH)/sum(LH))
		d.star <- c(d.star,sum(d*LH)/sum(LH))
		#if(i%%1000==0){cat(i,'\n')}
		}
	adjust <- rbind(g.star,d.star)
	rownames(adjust) <- c("g.star","d.star")
	adjust	
	} 

#fits the L/S model in the presence of missing data values

Beta.NA = function(y,X){
	des=X[!is.na(y),]
	y1=y[!is.na(y)]
	B <- solve(t(des)%*%des)%*%t(des)%*%y1
	B
	}
