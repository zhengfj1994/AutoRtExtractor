ResultSummary <- function(AutomatedRtExtractorList){
  ## 获取第一次出现真值的位置
  for (firstT in c(1:length(AutomatedRtExtractorList))){
    if (nrow(AutomatedRtExtractorList[[firstT]]$MSMSmatchedResult) > 0){
      colnames_MSMSmatchedResult <- colnames(AutomatedRtExtractorList[[firstT]]$MSMSmatchedResult)
      break
    }
  }
  tempMatrix <- as.data.frame(matrix(rep(NA,length(colnames_MSMSmatchedResult)),nrow = 1))
  colnames(tempMatrix) <- colnames_MSMSmatchedResult
  SummaryResult <- cbind(AutomatedRtExtractorList[[firstT]]$CompoundInformation[0,],AutomatedRtExtractorList[[firstT]]$MSMSmatchedResult[0,])

  pb <- txtProgressBar(style=3)
  for (i in c(1:length(AutomatedRtExtractorList))){
    setTxtProgressBar(pb, i/length(AutomatedRtExtractorList))
    ithMSMSmatchedResult <- AutomatedRtExtractorList[[i]]$MSMSmatchedResult
    ithCompoundInformation <- AutomatedRtExtractorList[[i]]$CompoundInformation
    if (is.null(ithMSMSmatchedResult)){
      ithRow <- cbind(ithCompoundInformation,tempMatrix)
    } else if (nrow(ithMSMSmatchedResult) == 0){
      ithRow <- cbind(ithCompoundInformation,tempMatrix)
    } else {
      ithRow <-cbind(ithCompoundInformation[rep(1,nrow(ithMSMSmatchedResult)),],ithMSMSmatchedResult)
    }
    SummaryResult <- rbind(SummaryResult,ithRow)
  }
  close(pb)
  return(SummaryResult)
}
