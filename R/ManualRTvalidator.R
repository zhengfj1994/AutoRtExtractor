##------------------------------------------------------------------------------
## 定义RTvalidator用于核对保留时间
##------------------------------------------------------------------------------
ManualRTvalidator <- function(AutomatedRtExtractorList, Intensity2 = 1000){
  library(patchwork)
  library(ggplot2)
  library(ggplotify)
  library(ggbump)
  library(dplyr)
  library(gWidgets2)
  library(ggrepel)
  ## Remove breakCommend
  if (exists("breakCommend")){
    remove(breakCommend)
  }

  startPoint <- 1
  while (length(grep("1_Not yet judged", AutomatedRtExtractorList[[startPoint]]$Addcuts[[1]]$RT$RTjudgment)) == 0){
    if (startPoint < length(AutomatedRtExtractorList)){
      startPoint <- startPoint + 1
    }
    else{
      break
    }
  }
  if (startPoint > 1){
    startPoint <- startPoint - 1
  }
  ## 获取所有的分子式
  AllMF <- as.data.frame(matrix(data = NA, nrow = length(AutomatedRtExtractorList), ncol = 2))
  colnames(AllMF) <- c("Mix", "MF")
  for (ith in c(1:length(AutomatedRtExtractorList))){
    AllMF$Mix[ith] <- unlist(strsplit(names(AutomatedRtExtractorList)[ith], split = "_"))[1]
    AllMF$MF[ith] <- AutomatedRtExtractorList[[ith]]$CompoundInformation$MolecularFormula
  }
  ## 循环list的每一个元素
  for (ith in c(startPoint:length(AutomatedRtExtractorList))){
    ## 判断breakCommend是否存在，如果存在，终止循环
    tryResult <- try(breakCommend == "Break", silent = T)
    if (!'try-error' %in% class(tryResult)){break}
    ## 判断是否为第一个出现的分子式
    IDwithSameMF <- which(AllMF$MF == AutomatedRtExtractorList[[ith]]$CompoundInformation$MolecularFormula & AllMF$Mix == unlist(strsplit(names(AutomatedRtExtractorList)[ith], split = "_"))[1])
    ## 如果是第一个出现的分子式，执行判断。如果不是，直接赋予第一个分子式的结果
    if (which(IDwithSameMF == ith) == 1){
      ## 获取Addcut信息
      ithAdductList <- AutomatedRtExtractorList[[ith]]$Addcuts
      ## 循环每一个Addcut
      for (jth in c(1:length(ithAdductList))){
        ## 判断breakCommend是否存在，如果存在，终止循环
        tryResult <- try(breakCommend == "Break", silent = T)
        if (!'try-error' %in% class(tryResult)){break}
        ## 获取第j个Addcut的信息
        jthAdductList <- ithAdductList[[jth]]
        ## 获取TIC
        TIC <- jthAdductList$TIC
        ## 获取DetectedPeaks
        DetectedPeaks <- jthAdductList$DetectedPeaks
        ## 获取RTmatrix
        RTmatrix <- jthAdductList$RT
        ## 如果全部为自动判断，跳过
        if (length(which(RTmatrix$RTjudgment == "1_Not yet judged")) == 0){
          next()
        }
        ## 如果只匹配得到两个峰，并且两个峰的保留时间偏差小于15s, 加起来的IE小于1.75，说明这两个峰是单独的两个峰，自动为真
        else if (length(which(RTmatrix$RTjudgment == "1_Not yet judged")) == 2 &
                 nrow(RTmatrix) == 2 &
                 max(RTmatrix$RT)-min(RTmatrix$RT) < 15 &
                 sum(RTmatrix$IE) < 2){
          RTmatrix$RTjudgment <- "2_Automated True"
        }
        ## 否则执行判断
        else{
          NotYetJudged <- which(RTmatrix$RTjudgment == "1_Not yet judged")
          ## 获取保留时间
          RetentionTime <- DetectedPeaks$RetentionTime[NotYetJudged[1]]
          ## 获取强度
          Intensity <- DetectedPeaks$Intensity[NotYetJudged[1]]
          ## 获取EIC
          EIC <- DetectedPeaks$EIC[NotYetJudged[1]]
          ## 拆分EIC为矩阵，用于画图
          EICmatrix <- strsplit(strsplit(EIC, split = ";")[[1]], split = " ")
          EICmatrix <- data.frame(matrix(as.numeric(unlist(EICmatrix)), byrow = T, ncol = 3))
          colnames(EICmatrix) <- c("scan", "scantime", "intensity")
          ## 获取IE
          IE <- DetectedPeaks$IE[NotYetJudged[1]]
          ## 获取预测曲线
          predicttedCurve <- DetectedPeaks$predicttedCurve[NotYetJudged[1]]
          ## 拆分预测曲线，用于画图
          predicttedCurveMatrix <- strsplit(strsplit(predicttedCurve, split = ";")[[1]], split = " ")
          predicttedCurveMatrix <- data.frame(matrix(as.numeric(unlist(predicttedCurveMatrix)), byrow = T, ncol = 2))
          colnames(predicttedCurveMatrix) <- c("scantime", "intensity")
          ## 获取残差
          residual <- DetectedPeaks$residual[NotYetJudged[1]]
          ## 获取同位素相似性
          isotopeSimilarity <- DetectedPeaks$isotopeSimilarity[NotYetJudged[1]]
          ## 获取同位素对齐矩阵
          isotopeAlignment <- DetectedPeaks$isotopeAlignment[NotYetJudged[1]]
          ## 如果加和离子形式为M+H或者M-H或者第一个为M+或M-，执行判断
          if ((ithAdductList[[jth]]$AdductInformation$Name == "M+H" |
               ithAdductList[[jth]]$AdductInformation$Name == "M-H") |
              (jth == 1 & (ithAdductList[[jth]]$AdductInformation$Name == "M+" |
                           ithAdductList[[jth]]$AdductInformation$Name == "M++" |
                           ithAdductList[[jth]]$AdductInformation$Name == "M-"))){
            ## 绘制TIC
            TICforPlot <- TIC
            TICforPlot$AutoOrManually <- NA
            TICforPlot$AutoOrManually[DetectedPeaks$Scan] <- RTmatrix$RTjudgment
            TICforPlot$order <- NA
            TICforPlot$order[DetectedPeaks$Scan] <- c(1:nrow(DetectedPeaks))
            TICforPlot <- na.omit(TICforPlot)

            p0 <- ggplot() +
              ggtitle(paste(ith, " | ", AutomatedRtExtractorList[[ith]]$CompoundInformation$English.Name, " | ", AutomatedRtExtractorList[[ith]]$CompoundInformation$MolecularFormula, " | ", ithAdductList[[jth]]$AdductInformation$Name)) +
              theme(plot.title = element_text(hjust = 0.5),
                    legend.title = element_blank(),
                    legend.position = "top") +
              geom_line(data = TIC, aes(x=scantime,y=intensity)) +
              geom_hline(yintercept = max(TIC$intensity)*0.1, color = "blue",linetype = 2) +
              geom_hline(yintercept = Intensity2, color = "red",linetype = 2) +
              geom_point(data = TICforPlot, aes(x = scantime, y = intensity, color = AutoOrManually)) +
              geom_text_repel(data = TICforPlot, aes(x = scantime, y = intensity, label = order)) +
              scale_fill_brewer(palette = "Set1")

            ## 绘制EIC
            p1 <- ggplot() +
              geom_line(data = TIC[which(TIC$scantime-RetentionTime+50 >= 0)[1]:rev(which((TIC$scantime-RetentionTime-50) <= 0))[1],], aes(x=scantime,y=intensity)) +
              geom_point(data = EICmatrix, aes(x=scantime[which.max(intensity)],y=intensity[which.max(intensity)])) +
              scale_x_continuous(guide = guide_axis(position = "bottom")) +
              annotate("text",
                       x = EICmatrix$scantime[which.max(EICmatrix$intensity)]-10,
                       y = EICmatrix$intensity[which.max(EICmatrix$intensity)],
                       label = paste("IE =",round(IE,2))) +
              geom_line(data = predicttedCurveMatrix, aes(x=scantime,y=intensity), color = "red",linetype = 2) +
              geom_vline(xintercept = EICmatrix$scantime[1], color = "blue",linetype = 2) +
              geom_vline(xintercept = rev(EICmatrix$scantime)[1], color = "blue",linetype = 2) +
              geom_point(data = predicttedCurveMatrix, aes(x=scantime[which.max(intensity)],y=intensity[which.max(intensity)]), color = "red") +
              annotate("text",
                       x = predicttedCurveMatrix$scantime[which.max(predicttedCurveMatrix$intensity)]+15,
                       y = predicttedCurveMatrix$intensity[which.max(predicttedCurveMatrix$intensity)]*1.05,
                       label = paste("residual =",round(residual,2)),
                       color = "red")
            ## 绘制同位素对比图
            if (!is.na(isotopeAlignment)){
              isotopeAlignmentMatrix <- strsplit(strsplit(isotopeAlignment, split = ";")[[1]], split = " ")
              isotopeAlignmentMatrix <- data.frame(matrix(as.numeric(unlist(isotopeAlignmentMatrix)), byrow = T, ncol = 3))
              colnames(isotopeAlignmentMatrix) <- c("mz","experimental","theoretical")
              plotCommand = "p2 <- ggplot(data = isotopeAlignmentMatrix) + labs(x='mz',y='intensity') + scale_x_continuous(breaks = pretty(range(isotopeAlignmentMatrix$mz), 5)) + scale_y_continuous(breaks = c(-100,-50,0,50,100), labels = c(100,50,0,50,100)) + geom_abline(slope = 0, intercept = 0) + annotate('text', x = (min(pretty(range(isotopeAlignmentMatrix$mz), 5)) + max(pretty(range(isotopeAlignmentMatrix$mz), 5)))/2, y = 95, label = 'Experiment', color = 'blue') + annotate('text', x = (min(pretty(range(isotopeAlignmentMatrix$mz), 5)) + max(pretty(range(isotopeAlignmentMatrix$mz), 5)))/2, y = -95, label = 'Theory', color = 'red') + annotate('text', x = (min(pretty(range(isotopeAlignmentMatrix$mz), 5)) + max(pretty(range(isotopeAlignmentMatrix$mz), 5)))/2, y = 75, label = paste('Similarity =', round(isotopeSimilarity, 2)))"
              for (i in c(1:nrow(isotopeAlignmentMatrix))){
                ithIons <- paste0(" + geom_segment(aes(x = mz[",i,"], y = 0, xend = mz[",i,"], yend = experimental[",i,"]), color = 'blue') + geom_segment(aes(x = mz[",i,"], y = 0, xend = mz[",i,"], yend = -1*theoretical[",i,"]), color = 'red')")
                plotCommand <- paste0(plotCommand,ithIons)
              }
              eval(parse(text = plotCommand))
            } else{
              p2 = ggplot()
            }
            ## 将三个图排版显示
            plot(p0 / (p1 | p2))
            ## 暂停一秒钟再做判断
            Sys.sleep(0.5)
            ## 生成判断对话框
            print(paste(ith, "in", length(AutomatedRtExtractorList)))
            RTjudgment <- ginput(msg = "If there is x true peaks, please input x.",
                                 text = "1", title = "RTjudgment", icon = "question",
                                 parent = c(100,100), toolkit = guiToolkit())
            ## 信息写入RTmatrix
            ## 如果对话框输入为e，结束整个循环
            if(RTjudgment == "e"){
              breakCommend <- "Break"
              break
            }else if (RTjudgment == "0"){
              RTmatrix$RTjudgment[NotYetJudged] <- "6_Manually False"
            }else if (grepl("\\?", RTjudgment) | grepl("？", RTjudgment)){
              RTmatrix$RTjudgment[as.numeric(unlist(strsplit(RTjudgment, split = " "))[-1])] <- "?_DoubtfulTrue"
            }else {
              ManuallyTrueID <- as.numeric(unlist(strsplit(RTjudgment, split = " ")))
              if (max(ManuallyTrueID) > nrow(RTmatrix)){
                return(AutomatedRtExtractorList)
                stop("The input is wrong!")
              }
              else{
                RTmatrix$RTjudgment[ManuallyTrueID[which(ManuallyTrueID > 0)]] <- "5_Manually True"
                RTmatrix$RTjudgment[abs(ManuallyTrueID[which(ManuallyTrueID < 0)])] <- "6_Manually False"
                RTmatrix$RTjudgment[setdiff(NotYetJudged, abs(ManuallyTrueID))] <- "7_No Manually Judge"
              }
            }
          } else {
            next()
          }
          ## RTmatrix写入结果list
          AutomatedRtExtractorList[[ith]]$Addcuts[[jth]]$RT <- RTmatrix
        }
      }
    } else {
      ## 获取Addcut信息
      ithAdductList <- AutomatedRtExtractorList[[ith]]$Addcuts
      for (jth in c(1:length(ithAdductList))){
        ## 循环每一个Addcut
        AutomatedRtExtractorList[[ith]]$Addcuts[[jth]]$RT <- AutomatedRtExtractorList[[IDwithSameMF[1]]]$Addcuts[[jth]]$RT
      }
    }
  }
  return(AutomatedRtExtractorList)
}
##------------------------------------------------------------------------------
