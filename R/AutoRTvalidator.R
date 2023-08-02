##------------------------------------------------------------------------------
## 定义RTvalidator用于核对保留时间
##------------------------------------------------------------------------------
AutoRTvalidator <- function(AutomatedRtExtractorList, Intensity1 = 2000, Intensity2 = 1000, IE1 = 0.5){
  library(patchwork)
  library(ggplot2)
  library(ggplotify)
  library(ggbump)
  library(dplyr)
  library(gWidgets2)
  library(ggrepel)

  ## 循环list的每一个元素
  for (ith in c(1:length(AutomatedRtExtractorList))){
    ## 获取Addcut信息
    ithAdductList <- AutomatedRtExtractorList[[ith]]$Addcuts
    ## 循环每一个Addcut
    for (jth in c(1:length(ithAdductList))){
      ## 获取第j个Addcut的信息
      jthAdductList <- ithAdductList[[jth]]
      ## 获取TIC
      TIC <- jthAdductList$TIC
      ## 获取DetectedPeaks
      DetectedPeaks <- jthAdductList$DetectedPeaks
      ## 如果DetectedPeaks为空，RTMatrix也为空
      if (nrow(DetectedPeaks) == 0){
        RTMatrix <- data.frame(RT = NA,
                               Intensity = NA,
                               EIC = NA,
                               IE = NA,
                               predicttedCurve = NA,
                               residual = NA,
                               isotopeSimilarity = NA,
                               isotopeAlignment = NA,
                               RTjudgment = NA)
        ## 将RTMatrix写入list
        AutomatedRtExtractorList[[ith]]$Addcuts[[jth]]$RT <- RTMatrix
      }else{
        ## 生成和DetectedPeaks相同行数的RTMatrix
        RTMatrix <- data.frame(RT = rep(NA,nrow(DetectedPeaks)),
                               Intensity = NA,
                               EIC = NA,
                               IE = NA,
                               predicttedCurve = NA,
                               residual = NA,
                               isotopeSimilarity = NA,
                               isotopeAlignment = NA,
                               RTjudgment = NA)
        ## 循环DetectedPeaks的每一行
        for (kth in c(1:nrow(DetectedPeaks))){
          if (kth == 1){ ## 对于DetectedPeaks的第一行，执行判断
            ## 获取EIC矩阵
            EICmatrix <- do.call("rbind", strsplit(unlist(strsplit(DetectedPeaks$EIC[kth], split = ";")), split = " "))
            ## 获取保留时间
            RetentionTime <- DetectedPeaks$RetentionTime[kth]
            LeftRetentionTime <- DetectedPeaks$LeftRetentionTime[kth]
            ## 如果LeftRetentionTime为NA，分配其为EIC最小保留时间
            if (is.na(LeftRetentionTime)){
              LeftRetentionTime <- as.numeric(EICmatrix[1,2])
            }
            RightRetentionTime <- DetectedPeaks$RightRetentionTime[kth]
            ## 如果RightRetentionTime为NA，分配其为EIC最大保留时间
            if (is.na(RightRetentionTime)){
              RightRetentionTime <- as.numeric(EICmatrix[nrow(EICmatrix),2])
            }
            ## 获取强度
            Intensity <- DetectedPeaks$Intensity[kth]
            LeftIntensity <- DetectedPeaks$LeftIntensity[kth]
            ## 如果LeftIntensity为NA，分配其为EIC最小保留时间下的强度
            if (is.na(LeftIntensity)){
              LeftIntensity <- as.numeric(EICmatrix[1,3])
            }
            RightIntensity <- DetectedPeaks$RightIntensity[kth]
            ## 如果RightIntensity为NA，分配其为EIC最大保留时间下的强度
            if (is.na(RightIntensity)){
              RightIntensity <- as.numeric(EICmatrix[nrow(EICmatrix),3])
            }
            ## 获取IE
            IE <- DetectedPeaks$IE[kth]
            ## 执行判断
            if (Intensity > Intensity1){
              if (IE > 2){
                RTjudgment <- "3_Automated False"
              }else if (IE > 1){
                if (Intensity/LeftIntensity > 100 &
                    Intensity/RightIntensity > 100 &
                    RightRetentionTime - LeftRetentionTime < 25){
                  RTjudgment <- "1_Not yet judged"
                }else {
                  RTjudgment <- "3_Automated False"
                }
              }else if (IE > IE1){
                if (Intensity/LeftIntensity > 100 &
                    Intensity/RightIntensity > 100 &
                    RightRetentionTime - LeftRetentionTime < 25){
                  RTjudgment <- "2_Automated True"
                }else {
                  RTjudgment <- "1_Not yet judged"
                }
              }else {
                RTjudgment <- "2_Automated True"
              }
            }else if (Intensity > Intensity2){
              if (IE > 0.75){
                RTjudgment <- "3_Automated False"
              }else if (IE > IE1){
                if (Intensity/LeftIntensity > 100 &
                    Intensity/RightIntensity > 100 &
                    RightRetentionTime - LeftRetentionTime < 25){
                  RTjudgment <- "1_Not yet judged"
                }else {
                  RTjudgment <- "3_Automated False"
                }
              }else if (IE > 0.2){
                if (Intensity/LeftIntensity > 100 &
                    Intensity/RightIntensity > 100 &
                    RightRetentionTime - LeftRetentionTime < 25){
                  RTjudgment <- "2_Automated True"
                }else {
                  RTjudgment <- "1_Not yet judged"
                }
              }else {
                RTjudgment <- "2_Automated True"
              }
            }else {
              RTjudgment <- "3_Automated False"
            }

            ## 信息写入RTMatrix
            RTMatrix$RT[kth] <- RetentionTime
            RTMatrix$Intensity[kth] <- Intensity
            RTMatrix$EIC[kth] <- DetectedPeaks$EIC[kth]
            RTMatrix$IE[kth] <- IE
            RTMatrix$predicttedCurve[kth] <- DetectedPeaks$predicttedCurve[kth]
            RTMatrix$residual[kth] <- DetectedPeaks$residual[kth]
            RTMatrix$isotopeSimilarity[kth] <- DetectedPeaks$isotopeSimilarity[kth]
            RTMatrix$isotopeAlignment[kth] <- DetectedPeaks$isotopeAlignment[kth]
            RTMatrix$RTjudgment[kth] <- RTjudgment
          }
          else { ## 对于DetectedPeaks的第二行及以下，执行判断
            ## 获取EIC矩阵
            EICmatrix <- do.call("rbind", strsplit(unlist(strsplit(DetectedPeaks$EIC[kth], split = ";")), split = " "))
            ## 获取保留时间
            RetentionTime <- DetectedPeaks$RetentionTime[kth]
            LeftRetentionTime <- DetectedPeaks$LeftRetentionTime[kth]
            ## 如果LeftRetentionTime为NA，分配其为EIC最小保留时间
            if (is.na(LeftRetentionTime)){
              LeftRetentionTime <- as.numeric(EICmatrix[1,2])
            }
            RightRetentionTime <- DetectedPeaks$RightRetentionTime[kth]
            ## 如果RightRetentionTime为NA，分配其为EIC最大保留时间
            if (is.na(RightRetentionTime)){
              RightRetentionTime <- as.numeric(EICmatrix[nrow(EICmatrix),2])
            }
            ## 获取强度
            Intensity <- DetectedPeaks$Intensity[kth]
            LeftIntensity <- DetectedPeaks$LeftIntensity[kth]
            ## 如果LeftIntensity为NA，分配其为EIC最小保留时间下的强度
            if (is.na(LeftIntensity)){
              LeftIntensity <- as.numeric(EICmatrix[1,3])
            }
            RightIntensity <- DetectedPeaks$RightIntensity[kth]
            ## 如果RightIntensity为NA，分配其为EIC最大保留时间下的强度
            if (is.na(RightIntensity)){
              RightIntensity <- as.numeric(EICmatrix[nrow(EICmatrix),3])
            }
            ## 获取信息熵
            IE <- DetectedPeaks$IE[kth]
            ## 如果强度大于最大峰强度的指定比例（十分之一）或者小于指定比例，但是最大强度峰判定为假，执行下述语句
            if (Intensity/max(DetectedPeaks$Intensity) > 0.1 |
                (Intensity/max(DetectedPeaks$Intensity) <= 0.1 &
                 RTMatrix$RTjudgment[1] != "2_Automated True")){
              ## 执行判断
              if (Intensity > Intensity1){
                if (IE > 2){
                  RTjudgment <- "3_Automated False"
                }else if (IE > 1){
                  if (Intensity/LeftIntensity > 100 &
                      Intensity/RightIntensity > 100 &
                      RightRetentionTime - LeftRetentionTime < 25){
                    RTjudgment <- "1_Not yet judged"
                  }else {
                    RTjudgment <- "3_Automated False"
                  }
                }else if (IE > IE1){
                  if (Intensity/LeftIntensity > 100 &
                      Intensity/RightIntensity > 100 &
                      RightRetentionTime - LeftRetentionTime < 25){
                    RTjudgment <- "2_Automated True"
                  }else {
                    RTjudgment <- "1_Not yet judged"
                  }
                }else {
                  RTjudgment <- "2_Automated True"
                }
              }else if (Intensity > Intensity2){
                if (IE > 0.75){
                  RTjudgment <- "3_Automated False"
                }else if (IE > IE1){
                  if (Intensity/LeftIntensity > 100 &
                      Intensity/RightIntensity > 100 &
                      RightRetentionTime - LeftRetentionTime < 25){
                    RTjudgment <- "1_Not yet judged"
                  }else {
                    RTjudgment <- "3_Automated False"
                  }
                }else if (IE > 0.2){
                  if (Intensity/LeftIntensity > 100 &
                      Intensity/RightIntensity > 100 &
                      RightRetentionTime - LeftRetentionTime < 25){
                    RTjudgment <- "2_Automated True"
                  }else {
                    RTjudgment <- "1_Not yet judged"
                  }
                }else {
                  RTjudgment <- "2_Automated True"
                }
              }else {
                RTjudgment <- "3_Automated False"
              }
            }else {
              ## 如果强度小于最大峰强度的十分之一，并且最大强度峰已被判定为真，则不对此峰进行判断
              RTjudgment <- "4_There is already high intensity true peak"
            }
            ## 信息写入RTMatrix
            RTMatrix$RT[kth] <- RetentionTime
            RTMatrix$Intensity[kth] <- Intensity
            RTMatrix$EIC[kth] <- DetectedPeaks$EIC[kth]
            RTMatrix$IE[kth] <- IE
            RTMatrix$predicttedCurve[kth] <- DetectedPeaks$predicttedCurve[kth]
            RTMatrix$residual[kth] <- DetectedPeaks$residual[kth]
            RTMatrix$isotopeSimilarity[kth] <- DetectedPeaks$isotopeSimilarity[kth]
            RTMatrix$isotopeAlignment[kth] <- DetectedPeaks$isotopeAlignment[kth]
            RTMatrix$RTjudgment[kth] <- RTjudgment
          }
        }
        ## RTMatrix写入结果list
        AutomatedRtExtractorList[[ith]]$Addcuts[[jth]]$RT <- RTMatrix
      }
    }
  }
  return(AutomatedRtExtractorList)
}
##------------------------------------------------------------------------------
