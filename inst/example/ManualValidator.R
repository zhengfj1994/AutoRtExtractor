##------------------------------------------------------------------------------
## 定义RTvalidator用于核对保留时间
##------------------------------------------------------------------------------
ManualRTvalidator <- function(AutomatedRtExtractorList){
  library(patchwork)
  library(ggplot2)
  library(ggplotify)
  library(ggbump)
  library(dplyr)
  library(gWidgets2)
  library(ggrepel)
  ## Remove breakCommend
  remove(breakCommend)

  ## 循环list的每一个元素
  for (ith in c(1:length(AutomatedRtExtractorList))){
    ## 判断breakCommend是否存在，如果存在，终止循环
    tryResult <- try(breakCommend == "Break", silent = T)
    if (!'try-error' %in% class(tryResult)){break}
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

      ## 获取保留时间
      RetentionTime <- DetectedPeaks$RetentionTime[1]
      ## 获取强度
      Intensity <- DetectedPeaks$Intensity[1]
      ## 获取EIC
      EIC <- DetectedPeaks$EIC[1]
      ## 拆分EIC为矩阵，用于画图
      EICmatrix <- strsplit(strsplit(EIC, split = ";")[[1]], split = " ")
      EICmatrix <- data.frame(matrix(as.numeric(unlist(EICmatrix)), byrow = T, ncol = 3))
      colnames(EICmatrix) <- c("scan", "scantime", "intensity")
      ## 获取IE
      IE <- DetectedPeaks$IE[1]
      ## 获取预测曲线
      predicttedCurve <- DetectedPeaks$predicttedCurve[1]
      ## 拆分预测曲线，用于画图
      predicttedCurveMatrix <- strsplit(strsplit(predicttedCurve, split = ";")[[1]], split = " ")
      predicttedCurveMatrix <- data.frame(matrix(as.numeric(unlist(predicttedCurveMatrix)), byrow = T, ncol = 2))
      colnames(predicttedCurveMatrix) <- c("scan", "intensity")
      ## 获取残差
      residual <- DetectedPeaks$residual[1]
      ## 获取同位素相似性
      isotopeSimilarity <- DetectedPeaks$isotopeSimilarity[1]
      ## 获取同位素对齐矩阵
      isotopeAlignment <- DetectedPeaks$isotopeAlignment[1]

      ## 绘制TIC
      TICforPlot <- TIC[,c(1,3)]
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
        geom_line(data = TIC, aes(x=scan,y=intensity)) +
        geom_hline(yintercept = max(TIC$intensity)*0.1, color = "blue",linetype = 2) +
        geom_hline(yintercept = 1000, color = "red",linetype = 2) +
        geom_point(data = TICforPlot, aes(x = scan, y = intensity, color = AutoOrManually)) +
        geom_text_repel(data = TICforPlot, aes(x = scan, y = intensity, label = order)) +
        scale_fill_brewer(palette = "Set1")

      ## 绘制EIC
      p1 <- ggplot() +
        geom_line(data = TIC[(which(TIC$scantime == RetentionTime) - 50):(which(TIC$scantime == RetentionTime) + 50),], aes(x=scan,y=intensity)) +
        geom_point(data = EICmatrix, aes(x=scan[which.max(intensity)],y=intensity[which.max(intensity)])) +
        scale_x_continuous(guide = guide_axis(position = "bottom")) +
        annotate("text",
                 x = EICmatrix$scan[which.max(EICmatrix$intensity)]-10,
                 y = EICmatrix$intensity[which.max(EICmatrix$intensity)],
                 label = paste("IE =",round(IE,2))) +
        geom_line(data = predicttedCurveMatrix, aes(x=scan,y=intensity), color = "red",linetype = 2) +
        geom_vline(xintercept = EICmatrix$scan[1], color = "blue",linetype = 2) +
        geom_vline(xintercept = rev(EICmatrix$scan)[1], color = "blue",linetype = 2) +
        geom_point(data = predicttedCurveMatrix, aes(x=scan[which.max(intensity)],y=intensity[which.max(intensity)]), color = "red") +
        annotate("text",
                 x = predicttedCurveMatrix$scan[which.max(predicttedCurveMatrix$intensity)]+15,
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
      Sys.sleep(1)
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
      }else {
        next()
      }
    }
  }
}
##------------------------------------------------------------------------------
