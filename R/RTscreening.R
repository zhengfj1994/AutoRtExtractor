## 定义RTscreening用于筛选保留时间
##------------------------------------------------------------------------------
RTscreening <- function(AutomatedRtExtractorList){
  ## Only keep the true positive peaks.
  for (ith in c(1:length(AutomatedRtExtractorList))){
    AutomatedRtExtractorList[[ith]]$AddcutsRT <- AutomatedRtExtractorList[[ith]]$Addcuts
    ithAdductList <- AutomatedRtExtractorList[[ith]]$Addcuts
    for (jth in c(1:length(ithAdductList))){
      jthAdductList <- ithAdductList[[jth]]
      jthAdductList <- jthAdductList$RT[which(jthAdductList$RT$RTjudgment == "2_Automated True" | jthAdductList$RT$RTjudgment == "5_Manually True") ,]
      AutomatedRtExtractorList[[ith]]$AddcutsRT[[jth]] <- jthAdductList
    }
  }
  return(AutomatedRtExtractorList)
}
