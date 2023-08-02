##------------------------------------------------------------------------------
## 定义RTvalidator用于核对保留时间
##------------------------------------------------------------------------------
RTvalidator <- function(AutomatedRtExtractorList){
  library(patchwork)
  library(ggplot2)
  library(ggplotify)
  library(ggbump)
  library(dplyr)
  library(gWidgets2)
  library(ggrepel)
  ## Remove breakCommend
  remove(breakCommend)

  startPoint <- 1
  while (!is.null(AutomatedRtExtractorList[[startPoint]]$Addcuts[[1]]$RT)){
    startPoint <- startPoint + 1
  }
  if (startPoint > 1){
    startPoint <- startPoint - 1
  }

  ## 循环list的每一个元素
  for (ith in c(startPoint:length(AutomatedRtExtractorList))){
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
            ## 获取保留时间
            RetentionTime <- DetectedPeaks$RetentionTime[kth]
            ## 获取强度
            Intensity <- DetectedPeaks$Intensity[kth]
            ## 获取EIC
            EIC <- DetectedPeaks$EIC[kth]
            ## 拆分EIC为矩阵，用于画图
            EICmatrix <- strsplit(strsplit(EIC, split = ";")[[1]], split = " ")
            EICmatrix <- data.frame(matrix(as.numeric(unlist(EICmatrix)), byrow = T, ncol = 3))
            colnames(EICmatrix) <- c("scan", "scantime", "intensity")
            ## 获取IE
            IE <- DetectedPeaks$IE[kth]
            ## 获取预测曲线
            predicttedCurve <- DetectedPeaks$predicttedCurve[kth]
            ## 拆分预测曲线，用于画图
            predicttedCurveMatrix <- strsplit(strsplit(predicttedCurve, split = ";")[[1]], split = " ")
            predicttedCurveMatrix <- data.frame(matrix(as.numeric(unlist(predicttedCurveMatrix)), byrow = T, ncol = 2))
            colnames(predicttedCurveMatrix) <- c("scan", "intensity")
            ## 获取残差
            residual <- DetectedPeaks$residual[kth]
            ## 获取同位素相似性
            isotopeSimilarity <- DetectedPeaks$isotopeSimilarity[kth]
            ## 获取同位素对齐矩阵
            isotopeAlignment <- DetectedPeaks$isotopeAlignment[kth]
            ## 如果加和离子形式为M+H或者M-H或者第一个为M+或M-，执行判断
            if ((ithAdductList[[jth]]$AdductInformation$Name == "M+H" |
                ithAdductList[[jth]]$AdductInformation$Name == "M-H") |
                (jth == 1 & (ithAdductList[[jth]]$AdductInformation$Name == "M+" |
                             ithAdductList[[jth]]$AdductInformation$Name == "M++" |
                             ithAdductList[[jth]]$AdductInformation$Name == "M-"))){
              ## 强度大于阈值，同时信息熵小于阈值的时候，自动判断为真实信号（Automated True）
              if (Intensity > 1000 & IE < 0.5){
                RTjudgment <- "Automated True"
              }else if (Intensity < 1000 | IE > 2){
                ## 强度小于阈值，或者信息熵大于阈值的时候，自动判断为假信号（Automated False）
                RTjudgment <- "Automated False"
              }else { ## 都不满足，进行判断
                ## 绘制TIC
                TICforPlot <- TIC[,c(1,3)]
                TICforPlot$color <- NA
                TICforPlot$color[DetectedPeaks$Scan] <- "blue"
                TICforPlot$color[DetectedPeaks$Scan[kth]] <- "red"
                TICforPlot$color <- factor(TICforPlot$color)
                TICforPlot$order <- NA
                TICforPlot$order[DetectedPeaks$Scan] <- c(1:nrow(DetectedPeaks))
                TICforPlot <- na.omit(TICforPlot)

                p0 <- ggplot() + ggtitle(paste(ith, " | ",
                                               AutomatedRtExtractorList[[ith]]$CompoundInformation$English.Name, " | ",
                                               AutomatedRtExtractorList[[ith]]$CompoundInformation$MolecularFormula, " | ",
                                               ithAdductList[[jth]]$AdductInformation$Name)) +
                  theme(plot.title = element_text(hjust = 0.5)) +
                  geom_line(data = TIC, aes(x=scan,y=intensity)) +
                  geom_hline(yintercept = max(TIC$intensity)*0.1, color = "blue",linetype = 2) +
                  geom_hline(yintercept = 1000, color = "red",linetype = 2) +
                  geom_point(data = TICforPlot, aes(x = scan, y = intensity, color = color), show.legend = F) +
                  geom_text_repel(data = TICforPlot, aes(x = scan, y = intensity, label = order))

                ## 绘制EIC
                p1 <- ggplot() +
                  geom_line(data = EICmatrix, aes(x=scan,y=intensity)) +
                  geom_point(data = EICmatrix, aes(x=scan[which.max(intensity)],y=intensity[which.max(intensity)])) +
                  scale_x_continuous(guide = guide_axis(position = "bottom")) +
                  annotate("text",
                           x = EICmatrix$scan[which.max(EICmatrix$intensity)-10],
                           y = EICmatrix$intensity[which.max(EICmatrix$intensity)],
                           label = paste("IE =",round(IE,2))) +
                  geom_line(data = predicttedCurveMatrix, aes(x=scan,y=intensity), color = "red",linetype = 2) +
                  geom_point(data = predicttedCurveMatrix, aes(x=scan[which.max(intensity)],y=intensity[which.max(intensity)]), color = "red") +
                  annotate("text",
                           x = predicttedCurveMatrix$scan[which.max(predicttedCurveMatrix$intensity)+15],
                           y = predicttedCurveMatrix$intensity[which.max(predicttedCurveMatrix$intensity)],
                           label = paste("residual =",round(residual,2)),
                           color = "red")
                ## 绘制同位素对比图
                if (!is.na(isotopeAlignment)){
                  isotopeAlignmentMatrix <- strsplit(strsplit(isotopeAlignment, split = ";")[[1]], split = " ")
                  isotopeAlignmentMatrix <- data.frame(matrix(as.numeric(unlist(isotopeAlignmentMatrix)), byrow = T, ncol = 3))
                  colnames(isotopeAlignmentMatrix) <- c("mz","experimental","theoretical")
                  plotCommand =
                    "p2 <- ggplot(data = isotopeAlignmentMatrix) +
  labs(x='mz',y='intensity') +
  scale_x_continuous(breaks = pretty(range(isotopeAlignmentMatrix$mz), 5)) +
  scale_y_continuous(breaks = c(-100,-50,0,50,100), labels = c(100,50,0,50,100)) +
  geom_abline(slope = 0, intercept = 0) +
  annotate('text',
           x = (min(pretty(range(isotopeAlignmentMatrix$mz), 5)) + max(pretty(range(isotopeAlignmentMatrix$mz), 5)))/2,
           y = 95,
           label = 'Experiment',
           color = 'blue') +
  annotate('text',
           x = (min(pretty(range(isotopeAlignmentMatrix$mz), 5)) + max(pretty(range(isotopeAlignmentMatrix$mz), 5)))/2,
           y = -95,
           label = 'Theory',
           color = 'red') +
  annotate('text',
           x = (min(pretty(range(isotopeAlignmentMatrix$mz), 5)) + max(pretty(range(isotopeAlignmentMatrix$mz), 5)))/2,
           y = 75,
           label = paste('Similarity =', round(isotopeSimilarity, 2)))"
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
                ## 如果对话框输入为e，结束整个循环
                if(RTjudgment == "e"){
                  breakCommend <- "Break"
                  break
                }
              }
            }else { ## 如果加和离子不是M+H或者M-H，只进行自动判断，不进行人工判断
              if (Intensity > 1000 & IE < 0.5){
                ## 强度大于阈值，同时信息熵小于阈值的时候，自动判断为真实信号（Automated True）
                RTjudgment <- "Automated True"
              }else if (Intensity < 1000 | IE > 2){
                ## 强度小于阈值，或者信息熵大于阈值的时候，自动判断为假信号（Automated False）
                RTjudgment <- "Automated False"
              }else {
                ## 都不满足，进行判断
                RTjudgment <- "Not yet judged"
              }
            }
            ## 信息写入RTMatrix
            RTMatrix$RT[kth] <- RetentionTime
            RTMatrix$Intensity[kth] <- Intensity
            RTMatrix$EIC[kth] <- EIC
            RTMatrix$IE[kth] <- IE
            RTMatrix$predicttedCurve[kth] <- predicttedCurve
            RTMatrix$residual[kth] <- residual
            RTMatrix$isotopeSimilarity[kth] <- isotopeSimilarity
            RTMatrix$isotopeAlignment[kth] <- isotopeAlignment
            RTMatrix$RTjudgment[as.numeric(unlist(strsplit(RTjudgment, split = " ")))] <- "Manually True"
          }
          else { ## 对于DetectedPeaks的第二行及以下，执行判断
            ## 获取保留时间
            RetentionTime <- DetectedPeaks$RetentionTime[kth]
            ## 获取强度
            Intensity <- DetectedPeaks$Intensity[kth]
            ## 获取EIC
            EIC <- DetectedPeaks$EIC[kth]
            ## 拆分EIC为矩阵
            EICmatrix <- strsplit(strsplit(EIC, split = ";")[[1]], split = " ")
            EICmatrix <- data.frame(matrix(as.numeric(unlist(EICmatrix)), byrow = T, ncol = 3))
            colnames(EICmatrix) <- c("scan", "scantime", "intensity")
            ## 获取信息熵
            IE <- DetectedPeaks$IE[kth]
            ## 获取预测曲线
            predicttedCurve <- DetectedPeaks$predicttedCurve[kth]
            ## 拆分预测曲线为矩阵
            predicttedCurveMatrix <- strsplit(strsplit(predicttedCurve, split = ";")[[1]], split = " ")
            predicttedCurveMatrix <- data.frame(matrix(as.numeric(unlist(predicttedCurveMatrix)), byrow = T, ncol = 2))
            colnames(predicttedCurveMatrix) <- c("scan", "intensity")
            ## 获取残差
            residual <- DetectedPeaks$residual[kth]
            ## 获取同位素相似性
            isotopeSimilarity <- DetectedPeaks$isotopeSimilarity[kth]
            ## 获取同位素对齐矩阵
            isotopeAlignment <- DetectedPeaks$isotopeAlignment[kth]
            ## 如果强度大于最大峰强度的指定比例（十分之一）或者小于指定比例，但是最大强度峰判定为假，执行下述语句
            if (Intensity/max(DetectedPeaks$Intensity) > 0.1 |
                (Intensity/max(DetectedPeaks$Intensity) <= 0.1 &
                 RTMatrix$RTjudgment[1] != "Automated True" &
                 RTMatrix$RTjudgment[1] != "1" &
                 RTMatrix$RTjudgment[1] != "2" &
                 RTMatrix$RTjudgment[1] != "3")){
              ## 如果加和离子形式为M+H或者M-H或者第一个为M+或M-，执行判断
              if ((ithAdductList[[jth]]$AdductInformation$Name == "M+H" |
                   ithAdductList[[jth]]$AdductInformation$Name == "M-H") |
                  (jth == 1 & (ithAdductList[[jth]]$AdductInformation$Name == "M+" |
                               ithAdductList[[jth]]$AdductInformation$Name == "M++" |
                               ithAdductList[[jth]]$AdductInformation$Name == "M-"))){
                if (Intensity > 1000 & IE < 0.5){
                  ## 强度大于阈值，同时信息熵小于阈值的时候，自动判断为真实信号（Automated True）
                  RTjudgment <- "Automated True"
                }else if (Intensity < 1000 | IE > 2){
                  ## 强度小于阈值，或者信息熵大于阈值的时候，自动判断为假信号（Automated False）
                  RTjudgment <- "Automated False"
                }else { ## 都不满足，进行判断
                  ## 绘制TIC
                  TICforPlot <- TIC[,c(1,3)]
                  TICforPlot$color <- NA
                  TICforPlot$color[DetectedPeaks$Scan] <- "blue"
                  TICforPlot$color[DetectedPeaks$Scan[kth]] <- "red"
                  TICforPlot$color <- factor(TICforPlot$color)
                  TICforPlot$order <- NA
                  TICforPlot$order[DetectedPeaks$Scan] <- c(1:nrow(DetectedPeaks))
                  TICforPlot <- na.omit(TICforPlot)

                  p0 <- ggplot() + ggtitle(paste(ith, " | ",
                                                 AutomatedRtExtractorList[[ith]]$CompoundInformation$English.Name, " | ",
                                                 AutomatedRtExtractorList[[ith]]$CompoundInformation$MolecularFormula, " | ",
                                                 ithAdductList[[jth]]$AdductInformation$Name)) +
                    theme(plot.title = element_text(hjust = 0.5)) +
                    geom_line(data = TIC, aes(x=scan,y=intensity)) +
                    geom_hline(yintercept = max(TIC$intensity)*0.1, color = "blue",linetype = 2) +
                    geom_hline(yintercept = 1000, color = "red",linetype = 2) +
                    geom_point(data = TICforPlot, aes(x = scan, y = intensity, color = color), show.legend = F) +
                    geom_text_repel(data = TICforPlot, aes(x = scan, y = intensity, label = order))
                  ## 绘制EIC
                  p1 <- ggplot() +
                    geom_line(data = EICmatrix, aes(x=scan,y=intensity)) +
                    geom_point(data = EICmatrix, aes(x=scan[which.max(intensity)],y=intensity[which.max(intensity)])) +
                    scale_x_continuous(guide = guide_axis(position = "bottom")) +
                    annotate("text",
                             x = EICmatrix$scan[which.max(EICmatrix$intensity)-10],
                             y = EICmatrix$intensity[which.max(EICmatrix$intensity)],
                             label = paste("IE =",round(IE,2))) +
                    geom_line(data = predicttedCurveMatrix, aes(x=scan,y=intensity), color = "red",linetype = 2) +
                    geom_point(data = predicttedCurveMatrix, aes(x=scan[which.max(intensity)],y=intensity[which.max(intensity)]), color = "red") +
                    annotate("text",
                             x = predicttedCurveMatrix$scan[which.max(predicttedCurveMatrix$intensity)+15],
                             y = predicttedCurveMatrix$intensity[which.max(predicttedCurveMatrix$intensity)],
                             label = paste("residual =",round(residual,2)),
                             color = "red")
                  ## 绘制同位素对比图
                  if (!is.na(isotopeAlignment)){
                    isotopeAlignmentMatrix <- strsplit(strsplit(isotopeAlignment, split = ";")[[1]], split = " ")
                    isotopeAlignmentMatrix <- data.frame(matrix(as.numeric(unlist(isotopeAlignmentMatrix)), byrow = T, ncol = 3))
                    colnames(isotopeAlignmentMatrix) <- c("mz","experimental","theoretical")
                    plotCommand =
                      "p2 <- ggplot(data = isotopeAlignmentMatrix) +
  labs(x='mz',y='intensity') +
  scale_x_continuous(breaks = pretty(range(isotopeAlignmentMatrix$mz), 5)) +
  scale_y_continuous(breaks = c(-100,-50,0,50,100), labels = c(100,50,0,50,100)) +
  geom_abline(slope = 0, intercept = 0) +
  annotate('text',
           x = (min(pretty(range(isotopeAlignmentMatrix$mz), 5)) + max(pretty(range(isotopeAlignmentMatrix$mz), 5)))/2,
           y = 95,
           label = 'Experiment',
           color = 'blue') +
  annotate('text',
           x = (min(pretty(range(isotopeAlignmentMatrix$mz), 5)) + max(pretty(range(isotopeAlignmentMatrix$mz), 5)))/2,
           y = -95,
           label = 'Theory',
           color = 'red') +
  annotate('text',
           x = (min(pretty(range(isotopeAlignmentMatrix$mz), 5)) + max(pretty(range(isotopeAlignmentMatrix$mz), 5)))/2,
           y = 75,
           label = paste('Similarity =', round(isotopeSimilarity, 2)))"
                    for (i in c(1:nrow(isotopeAlignmentMatrix))){
                      ithIons <- paste0(" + geom_segment(aes(x = mz[",i,"], y = 0, xend = mz[",i,"], yend = experimental[",i,"]), color = 'blue') + geom_segment(aes(x = mz[",i,"], y = 0, xend = mz[",i,"], yend = -1*theoretical[",i,"]), color = 'red')")
                      plotCommand <- paste0(plotCommand,ithIons)
                    }
                    eval(parse(text = plotCommand))
                  } else{
                    p2 = ggplot()
                  }
                  ## 图片排版显示
                  plot(p0 / (p1 | p2))
                  ## 暂停一秒再执行判断
                  Sys.sleep(1)
                  ## 生成对话框执行判断
                  print(paste(ith, "in", length(AutomatedRtExtractorList)))
                  RTjudgment <- ginput(msg = "If there is x true peaks, please input x.",
                                       text = "1", title = "RTjudgment", icon = "question",
                                       parent = c(100,100), toolkit = guiToolkit())
                  ## 输入e终止所有循环
                  if(RTjudgment == "e"){
                    breakCommend <- "Break"
                    break
                  }
                }
              }else { ## 如果加和离子形式不是M+H或者M-H，不进行人工判断
                if (Intensity > 1000 & IE < 0.5){
                  ## 强度大于阈值，同时信息熵小于阈值的时候，自动判断为真实信号（Automated True）
                  RTjudgment <- "Automated True"
                }else if (Intensity < 1000 | IE > 2){
                  ## 强度小于阈值，或者信息熵大于阈值的时候，自动判断为假信号（Automated False）
                  RTjudgment <- "Automated False"
                }else {
                  ## 都不满足，进行判断
                  RTjudgment <- "Not yet judged"
                }
              }
            }else {
              ## 如果强度小于最大峰强度的十分之一，并且最大强度峰已被判定为真，则不对此峰进行判断
              RTjudgment <- "There is already high intensity true peak"
            }
            ## 信息写入RTMatrix
            RTMatrix$RT[kth] <- RetentionTime
            RTMatrix$Intensity[kth] <- Intensity
            RTMatrix$EIC[kth] <- EIC
            RTMatrix$IE[kth] <- IE
            RTMatrix$predicttedCurve[kth] <- predicttedCurve
            RTMatrix$residual[kth] <- residual
            RTMatrix$isotopeSimilarity[kth] <- isotopeSimilarity
            RTMatrix$isotopeAlignment[kth] <- isotopeAlignment
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
