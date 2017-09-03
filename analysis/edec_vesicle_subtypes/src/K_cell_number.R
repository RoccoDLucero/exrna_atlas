#functions

pMatrix.min <- function(A, B) {
  # finds the permutation P of A such that ||PA - B|| is minimum
  # in Frobenius norm
  # Uses the linear-sum assignment problem (LSAP) solver
  # in the "clue" package
  # Returns P%*%A and the permutation vector `pvec' such that
  # A[pvec, ] is the permutation of A closest to B
  n <- nrow(A)
  D <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      #           D[j, i] <- sqrt(sum((B[j, ] - A[i, ])^2))
      D[j, i] <- (sum((B[j, ] - A[i, ])^2))  # this is better
    } }
  vec <- c(solve_LSAP(D))
  list(A=A[vec,], pvec=vec)
}


meth.cor.min <- function(E1, E2,i) {
  d = diag(1, i, i)
  c = cor(E1,E2)
  f= pMatrix.min(c,d)
  min(diag(f$A)) }


prop.cor.min <- function(E1, E2,i) {
  d = diag(1, i, i)
  com = intersect(row.names(E1),row.names(E2))
  c = cor(E1[com,],E2[com,])
  f= pMatrix.min(c,d) 
  min(diag(f$A))}

#find cell number

find.cell.number <- function(lower,upper,reps,probes,tum_df){
measure = data.frame(0,0,0)
for(i in lower:upper){ 
  prop = list()
  meth = list()
  for(k in 1:reps){
  num <- sample(1:ncol(tum_df), round((0.8*ncol(tum_df)),0), replace=F)
  tum = tum_df[,num]
  tum = as.matrix(tum)
  x = EDecStage1(methMixtureSamples = tum,cellTypeSpecificLoci = probes,nCts = i)
  prop[[k]] <- as.data.frame(x$proportions)
  meth[[k]] <- as.data.frame(x$methylation)}
  pro = c()
  met = c()
  comb = combn( 1:length(prop), 2)
  for(col in 1:ncol(comb)){
    met[[length(met)+1]] <- meth.cor.min(meth[[comb[1,col]]],meth[[comb[2,col]]],i)  
    pro[[length(pro)+1]] <-prop.cor.min(prop[[comb[1,col]]],prop[[comb[2,col]]],i) 
  }
  
  results = c(i,round(min(met),3),round(min(pro),3))
  measure = rbind(measure,results)
  }
measure = measure[2:nrow(measure),]
colnames(measure) = c("Cell_Number","Meth_Cor","Prop_Cor")
return(measure)
}
  

