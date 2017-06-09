require(gtools)

simMixes = function(pureCts,methVarMatrix,propsMeans,propsVarFactor,nSamples){
	allNoisyMethProfs = list()
	#For each cell type add Beta distributed noise to each locus to create nSamples
	#simulated methylation profiles
	for(j in 1:ncol(pureCts)){
		noisyMethProfs = matrix(0,nrow=nrow(pureCts),ncol=nSamples)
		#for each locus
		for(i in 1:nrow(pureCts)){
		    #For each locus, get cell-type specific average methylation value
			avg = pureCts[i,j] 
			if(avg == 1){
				avg = 0.999999
			}
			if(avg == 0){
				avg = 0.000001
			}
            #Compute maximal variance at given locus scaled by
			#values in methVarMatrix (5% for normal 10% for cancer)
			methVar = (avg*(1-avg))*methVarMatrix[i,j]
			
			#Shape parameters in Beta distribution 
			a = avg*(((avg*(1-avg))/methVar)-1)
			b = (1-avg)*(((avg*(1-avg))/methVar)-1)
			#Get a vector of Beta distributed samples from the estimated
			#average methylation at current locus
			noisyMethProfs[i,] = rbeta(nSamples,a,b) 
		}
		allNoisyMethProfs[[j]] = noisyMethProfs
	}
	#Create random mixtures from the previously generated noisy methylation profiles
	#by adding noise to the proportions of each cell type in the simulated mix
	# around the average proportions desired for the simulated mix
	props = rdirichlet(nSamples,propsVarFactor*propsMeans)
	mixes = matrix(0,nrow=nrow(pureCts),ncol=nSamples)
	for(i in 1:nSamples){
		mix = rep(0,nrow(pureCts))
		for(j in 1:ncol(pureCts)){
		    #create a mixed methylation profile for each locus based on
		    #pure cell type noisy values and simulated proportions
			mix = mix + props[i,j]*allNoisyMethProfs[[j]][,i]
		}
		mixes[,i] = mix 
	}
	result = list(mixes,props,allNoisyMethProfs)	
	names(result) = c("mixes","proportions","noisyMethProfs")
	return(result)
}

#Calls 'simMixes' to generate mixed cell-type methylation profiles
#and returns the simulated methyaltion profiles plus
#the true proportions of pure cell types that made up each simulated profile 
makeSimMixes = function(simPureCts,labels,possiblePropMeans,propVarFactor=100,methVarMatrix,nSamples=50){
	mixes = NULL
	props = matrix(0,nrow=nSamples,ncol=ncol(simPureCts))
	colnames(props) = colnames(simPureCts)
	for(i in 1:nSamples){
		propMeans = possiblePropMeans[sample(1:nrow(possiblePropMeans),1),]
		pureCts = NULL
		pureCtsNames = NULL
		#This is meant to ensure only one reference sample is used
		# for each cell type in the simulted mix  
		for(lab in unique(labels)){
			pureCt = sample(colnames(simPureCts)[labels==lab],1)
			pureCtsNames = c(pureCtsNames,pureCt)
			pureCts = cbind(pureCts,simPureCts[,pureCt])
		}
		sim = simMixes(pureCts,methVarMatrix,propMeans,propVarFactor,1)
		mixes = cbind(mixes,sim$mix)
		props[i,pureCtsNames] = sim$prop
	}
	colnames(mixes) = paste("Mix",1:nSamples,sep=".")
	rownames(mixes) = rownames(simPureCts)
  rownames(props) = paste("Mix",1:nSamples,sep=".")
	result = list(mixes,props)
	names(result) = c("mixes","proportions")
	return(result)
}