perform_t_tests_all_rows = function(dataGroup1,dataGroup2){
  nGroup1 = ncol(dataGroup1)
  nGroup2 = ncol(dataGroup2)
  dataAll = cbind(dataGroup1,dataGroup2)
  tTestWithErrorHandling = function(x){
    testResult = try(t.test(x[1:nGroup1],x[(nGroup1+1):(nGroup1+nGroup2)]),silent=TRUE);
    if(is.character(testResult)){
      warning(testResult)
      c(NA,NA,NA)
    }else{
      c(testResult$p.value,testResult$estimate)
    }
  }
  results = matrix(unlist(apply(dataAll,1,tTestWithErrorHandling)),ncol=3,byrow=TRUE)
  colnames(results) = c("P.value","Mean.group.1","Mean.group.2")
  rownames(results) = rownames(dataGroup1)
  results
}

perform_t_tests_all_classes_one_vs_rest = function(dataMatrix,classVector){
  if(ncol(dataMatrix)!=length(classVector)){
    stop("Number of columns of data matrix must be equal to the length of the class vector")
  }
  possibleClasses = unique(classVector)
  nClasses = length(possibleClasses)
  
  allPvalues = matrix(NA,nrow=nrow(dataMatrix),ncol=nClasses)
  allDiffMeans = matrix(NA,nrow=nrow(dataMatrix),ncol=nClasses)
  colnames(allPvalues) = possibleClasses
  rownames(allPvalues) = rownames(dataMatrix)
  colnames(allDiffMeans) = possibleClasses
  rownames(allDiffMeans) = rownames(dataMatrix)
  
  for(i in 1:nClasses){
    class = possibleClasses[i]
    resultTest = perform_t_tests_all_rows(dataMatrix[,classVector==class],dataMatrix[,classVector!=class])
    allPvalues[,i] = resultTest[,1]
    allDiffMeans[,i] = resultTest[,2]-resultTest[,3]
  }
  result = list(allPvalues,allDiffMeans)
  names(result) = c("P.Values","Difference.Between.Means")
  return(result)
}

perform_t_tests_all_classes_each_pair = function(dataMatrix,classVector){
  if(ncol(dataMatrix)!=length(classVector)){
    stop("Number of columns of data matrix must be equal to the length of the class vector")
  }
  possibleClasses = unique(classVector)
  nClasses = length(possibleClasses)
  
  allPValues = NULL
  allDiffMeans = NULL
  names = NULL 
  for(i in 1:(nClasses-1)){
    for(j in (i+1):nClasses){
      class1 = possibleClasses[i]
      class2 = possibleClasses[j]
      names = c(names,paste(class1,class2,sep="."))
      result = perform_t_tests_all_rows(dataMatrix[,classVector==class1],dataMatrix[,classVector==class2])
      allPValues = cbind(allPValues,result[,1])
      allDiffMeans = cbind(allDiffMeans,result[,2]/result[,3])
    }
  }
  colnames(allPValues) = names
  colnames(allDiffMeans) = names
  result = list(allPValues,allDiffMeans)
  names(result) = c("P.Values","Difference.Between.Means")
  return(result)
}