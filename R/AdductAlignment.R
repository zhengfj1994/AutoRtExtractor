## 定义AdductAlignment
##------------------------------------------------------------------------------
AdductAlignment <- function(AutomatedRtExtractorList, deltaRT = 12, intensity1 = 2000){
  ## 加载颜色的R包
  library(ggplot2)
  library(ggrepel)
  library(ggplotify)
  library(patchwork)
  library(dplyr)
  library(RColorBrewer)
  library(gWidgets2)
  library(ChemmineR)
  ## Remove breakCommend
  if (exists("breakCommend")){
    remove(breakCommend)
  }

  startPoint <- 1
  while (!is.null(AutomatedRtExtractorList[[startPoint]]$AdductAlignmentResult)){
    startPoint <- startPoint + 1
  }
  if (startPoint > 1){
    startPoint <- startPoint - 1
  }

  ## 生成9种颜色
  colorVector <- brewer.pal(12,"Paired")
  ## 循环列表
  for (ith in c(startPoint:length(AutomatedRtExtractorList))){
    ## 判断breakCommend是否存在，如果存在，终止循环
    tryResult <- try(breakCommend == "Break", silent = T)
    if (!'try-error' %in% class(tryResult)){break}
    ## 获取列表的Adduct
    ithAdductList <- AutomatedRtExtractorList[[ith]]$AddcutsRT
    ## 合并列表为矩阵
    ithAdductDF <- do.call("rbind", ithAdductList)
    ## 如果矩阵为0行，跳过
    if (nrow(ithAdductDF) == 0){
      ## 赋予矩阵一些信息
      ithAdductDF$addcut <- NULL
      ithAdductDF$addcutID <- NULL
      ithAdductDF$addcut_addcutID <- NULL
      ## 赋予矩阵M+H峰对齐列
      ithAdductDF$AlignNum <- NULL
      ithAdductDF$AlignID <- NULL
      ithAdductDF$AdductAlignmentJudge <- NULL
    }
    else { ## 如果矩阵不为0行
      ## 赋予矩阵一些信息
      ithAdductDF$addcut <- row.names(ithAdductDF)
      ithAdductDF$addcutID <- c(1:nrow(ithAdductDF))
      ithAdductDF$addcut_addcutID <- paste0(ithAdductDF$addcutID, ": ",ithAdductDF$addcut)
      ## 赋予矩阵M+H峰对齐列
      ithAdductDF$AlignNum <- NA
      ithAdductDF$AlignID <- NA
      for (ithRow_ithAdductDF in c(1:nrow(ithAdductDF))){
        ithAdductDF$AlignNum[ithRow_ithAdductDF] <- length(which(abs(ithAdductDF$RT - ithAdductDF$RT[ithRow_ithAdductDF]) < deltaRT))
        ithAdductDF$AlignID[ithRow_ithAdductDF] <- paste(which(abs(ithAdductDF$RT - ithAdductDF$RT[ithRow_ithAdductDF]) < deltaRT), collapse = " ")
      }
      ## 赋予矩阵新的结果列
      ithAdductDF$AdductAlignmentJudge <- NA
      ## 构建绘图用EIC矩阵
      EICmatrixAll <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 4))
      for (jth in c(1:nrow(ithAdductDF))){
        EICmatrix <- strsplit(strsplit(ithAdductDF$EIC[jth], split = ";")[[1]], split = " ")
        EICmatrix <- data.frame(matrix(as.numeric(unlist(EICmatrix)), byrow = T, ncol = 3))
        colnames(EICmatrix) <- c("scan", "scantime", "intensity")
        EICmatrix$adduct <- ithAdductDF$addcut_addcutID[jth]
        EICmatrixAll <- rbind(EICmatrixAll,EICmatrix)
      }
      colnames(EICmatrixAll) <- c("scan", "scantime", "intensity", "adduct")
      ## 构建绘图用EIC矩阵的最大值子矩阵
      maxEICmatrixAll <- aggregate(EICmatrixAll,by=list(adduct=EICmatrixAll$adduct),FUN=max)[,-1]
      maxEICmatrixAll$scan <- aggregate(EICmatrixAll$scan,by=list(adduct=EICmatrixAll$adduct),FUN=mean)$x
      maxEICmatrixAll$scantime <- aggregate(EICmatrixAll$scantime,by=list(adduct=EICmatrixAll$adduct),FUN=mean)$x
      ## 绘图
      # ithSmiles <- AutomatedRtExtractorList[[ith]]$CompoundInformation$CanonicalSMILES
      # ithSDF <- smiles2sdf(as(ithSmiles, "SMIset"))
      # p_mol <- as.ggplot(~openBabelPlot(ithSDF))
      ## 绘图
      p_line <- ggplot(EICmatrixAll) +
        ggtitle(label = paste(ith, "|", AutomatedRtExtractorList[[ith]]$CompoundInformation$MolecularFormula,"|", AutomatedRtExtractorList[[ith]]$CompoundInformation$CanonicalSMILES)) +
        geom_line(aes(x = EICmatrixAll[,2],
                      y = EICmatrixAll[,3],
                      color = EICmatrixAll[,4]),
                  size=0.75) +
        scale_color_manual(values = colorVector[1:nrow(ithAdductDF)]) +
        #坐标间距调整
        scale_x_continuous(breaks = pretty(range(EICmatrixAll[,2]), 5)) +
        scale_y_continuous(labels = scales::scientific,
                           breaks = pretty(range(EICmatrixAll[,3]), 5)) +
        #标注坐标名称、legend图例名称
        labs(x="scantime",y="intensity",color="adduct") +
        geom_text_repel(aes(label=adduct),maxEICmatrixAll,
                        x = maxEICmatrixAll$scantime,
                        y = maxEICmatrixAll$intensity,
                        size = 5, #注释文本的字体大小
                        box.padding = 0.5, #字到点的距离
                        point.padding = 0, #字到点的距离，点周围的空白宽度
                        min.segment.length = 0, #短线段可以省略
                        segment.color = "black", #segment.colour = NA, 不显示线段
                        show.legend = F
        ) +
        #改图片其他参数
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.title = element_text(size = 12,face = "italic"),
          legend.text = element_text(size = 12, face = "italic"),
          legend.position = 'right',
          legend.key.size = unit(0.5,'cm'),
          axis.text.x = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5),
          axis.text.y = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5),
          axis.title.x = element_text(size = 12,face = "bold", vjust = -0.5, hjust = 0.5),
          axis.title.y = element_text(size = 12,face = "bold", vjust = 1.2, hjust = 0.5),
          panel.background = element_rect(fill = "transparent",colour = "black"),
          plot.background = element_rect(fill = "transparent",colour = "white"))

      ## 判断部分---------------------------------------------------------------
      TemporaryJudgment <- NA
      ## 获取Adduct名字
      AdductNames <- do.call("rbind", strsplit(ithAdductDF$addcut, split = "\\."))[,1]
      ## 如果只有一个Adduct，并且是M+H/M+/M++/M-H/M-, 结果自动为真
      if (length(unique(AdductNames)) == 1 &
          (unique(AdductNames)[1] ==  "M+H" | unique(AdductNames)[1] == "M+" | unique(AdductNames)[1] == "M++" | unique(AdductNames)[1] == "M-H" | unique(AdductNames)[1] == "M-")){
        ithAdductDF$AdductAlignmentJudge[which.max(ithAdductDF$Intensity)] <- "AutoTrue"
      }
      ##  如果只有一个Adduct，并且不是M+H/M+/M++/M-H/M-, 结果自动为疑惑真
      else if (length(unique(AdductNames)) == 1){
        ithAdductDF$AdductAlignmentJudge[which.max(ithAdductDF$Intensity)] <- "DoubtfulTrue"
      }
      else { ## 如果有多个Adduct
        ## 获取最多的对齐结果
        MaxAlignNumID <- which(ithAdductDF$AlignNum == max(ithAdductDF$AlignNum))
        ##获取M+H所在ID
        MplusHID <- which(grepl("^M\\+H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M+H")
        ##获取M-H所在ID
        MreduceHID <- which(grepl("^M\\-H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M-H")
        ## 如果多个Adduct没有一个超过两个的分组
        if (max(ithAdductDF$AlignNum) == 1){
          ## 如果M+H只有一个,并且没有M-H
          if (length(which(grepl("^M\\+H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M+H")) == 1 &
              length(which(grepl("^M\\-H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M-H")) == 0){
            ## 如果唯一的M+H强度大于intensity1
            if(ithAdductDF$Intensity[which(grepl("^M\\+H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M+H")] > intensity1 |
               ithAdductDF$Intensity[which(grepl("^M\\+H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M+H")] > 0.5 * max(ithAdductDF$Intensity)){
              ithAdductDF$AdductAlignmentJudge[which(grepl("^M\\+H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M+H")] <- "AutoTrue"
            }
            ## 如果唯一的M+H强度小于等于intensity1
            else {
              TemporaryJudgment <- "NotSure"
            }
          }
          ## 如果M-H只有一个，并且没有M+H
          else if (length(which(grepl("^M\\+H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M+H")) == 0 &
                   length(which(grepl("^M\\-H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M-H")) == 1){
            ## 如果唯一的M-H强度大于intensity1
            if(ithAdductDF$Intensity[which(grepl("^M\\-H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M-H")] > intensity1 |
               ithAdductDF$Intensity[which(grepl("^M\\-H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M-H")] > 0.5 * max(ithAdductDF$Intensity)){
              ithAdductDF$AdductAlignmentJudge[which(grepl("^M\\-H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M-H")] <- "AutoTrue"
            }
            ## 如果唯一的M-H强度小于等于intensity1
            else {
              TemporaryJudgment <- "NotSure"
            }
          }
          ## 如果M+H没有，M-H也没有
          else if (length(which(grepl("^M\\+H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M+H")) == 0 &
                   length(which(grepl("^M\\-H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M-H")) == 0){
            ## 获得最大强度的Adduct
            AdductWithMaxIntensity <- ithAdductDF$addcut[which.max(ithAdductDF$Intensity)]
            ## 最大强度自动为怀疑真
            ithAdductDF$AdductAlignmentJudge[which.max(ithAdductDF$Intensity)] <- "DoubtfulTrue"
          }
          ## 如果M+H和M-H都有一个
          else if (length(which(grepl("^M\\+H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M+H")) == 1 &
                   length(which(grepl("^M\\-H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M-H")) == 1){
            ## 如果M+H峰强度大于10倍的M-H峰，M+H自动为真
            if (ithAdductDF$Intensity[which(grepl("^M\\+H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M+H")] >
                10 *ithAdductDF$Intensity[which(grepl("^M\\-H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M-H")]){
              ithAdductDF$AdductAlignmentJudge[which(grepl("^M\\+H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M+H")] <- "AutoTrue"
            }
            ## 如果M-H峰大于10倍M+H峰，M-H峰自动为真
            else if (ithAdductDF$Intensity[which(grepl("^M\\-H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M-H")] >
                     10 *ithAdductDF$Intensity[which(grepl("^M\\+H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M+H")]){
              ithAdductDF$AdductAlignmentJudge[which(grepl("^M\\-H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M-H")] <- "AutoTrue"
            }
            ## 如果M+H和M-H强度相差在10倍以内，都为真
            else {
              ithAdductDF$AdductAlignmentJudge[which(grepl("^M\\+H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M+H")] <- "DoubtfulTrue"
              ithAdductDF$AdductAlignmentJudge[which(grepl("^M\\-H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M-H")] <- "DoubtfulTrue"
            }
          }
          ## 如果有多个M+H而没有M-H
          else if (length(which(grepl("^M\\+H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M+H")) > 1 &
                   length(which(grepl("^M\\-H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M-H")) == 0){
            ## 获取M+H的ID
            tempMplusHID <- which(grepl("^M\\+H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M+H")
            ## 获取M+H的强度
            tempMplusHInt <- ithAdductDF$Intensity[tempMplusHID]
            ## 如果强度大于0.1倍的最大强度的M+H峰只有一个，这个最大峰自动为真
            if (length(which(tempMplusHInt > 0.1 * max(tempMplusHInt))) == 1){
              ithAdductDF$AdductAlignmentJudge[tempMplusHID[which(tempMplusHInt > 0.1 * max(tempMplusHInt))]] <- "AutoTrue"
            }
            ## 如果强度大于0.1倍的最大强度的M+H峰有多个，都为疑惑真
            else {
              ithAdductDF$AdductAlignmentJudge[tempMplusHID[which(tempMplusHInt > 0.1 * max(tempMplusHInt))]] <- "DoubtfulTrue"
            }
          }
          ## 如果有多个M-H而没有M+H
          else if (length(which(grepl("^M\\+H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M+H")) == 0 &
                   length(which(grepl("^M\\-H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M-H")) > 1){
            ## 获取M-H的ID
            tempMreduceHID <- which(grepl("^M\\-H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M-H")
            ## 获取M-H的强度
            tempMreduceHInt <- ithAdductDF$Intensity[tempMreduceHID]
            ## 如果强度大于0.1倍的最大强度的M-H峰只有一个，这个最大峰自动为真
            if (length(which(tempMreduceHInt > 0.1 * max(tempMreduceHInt))) == 1){
              ithAdductDF$AdductAlignmentJudge[tempMreduceHID[which(tempMreduceHInt > 0.1 * max(tempMreduceHInt))]] <- "AutoTrue"
            }
            ## 如果强度大于0.1倍的最大强度的M-H峰有多个，都为疑惑真
            else {
              ithAdductDF$AdductAlignmentJudge[tempMreduceHID[which(tempMreduceHInt > 0.1 * max(tempMreduceHInt))]] <- "DoubtfulTrue"
            }
          }
          ## 其他情况
          else {
            ## 获取M+H和M-H的ID
            tempMHID <- which(grepl("^M\\+H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M+H" | grepl("^M\\-H\\.",ithAdductDF$addcut) | ithAdductDF$addcut == "M-H")
            ## 获取M+H/M-H的强度
            tempMHInt <- ithAdductDF$Intensity[tempMHID]
            ## 如果强度大于0.1倍的最大强度的M+H/M-H峰只有一个，这个最大峰自动为真
            if (length(which(tempMHInt > 0.1 * max(tempMHInt))) == 1){
              ithAdductDF$AdductAlignmentJudge[tempMHID[which(tempMHInt > 0.1 * max(tempMHInt))]] <- "AutoTrue"
            }
            ## 如果强度大于0.1倍的最大强度的M+H/M-H峰有多个，都为疑惑真
            else {
              ithAdductDF$AdductAlignmentJudge[tempMHID[which(tempMHInt > 0.1 * max(tempMHInt))]] <- "DoubtfulTrue"
            }
          }
        }
        ## 如果最多的对齐结果只有一组
        else if (length(unique(ithAdductDF$AlignID[MaxAlignNumID])) == 1){
          AlignNumID <- as.numeric(unlist(strsplit(ithAdductDF$AlignID[MaxAlignNumID], split = " ")))
          ## 如果M+H和M-H都在这一组
          if (length(intersect(AlignNumID, MplusHID)) > 0 &
              length(intersect(AlignNumID, MreduceHID)) > 0){
            ithAdductDF$AdductAlignmentJudge[AlignNumID] <- "AutoTrue"
          }
          ## 如果有M+H或者M-H
          else if (length(intersect(AlignNumID, MplusHID)) > 0 |
                   length(intersect(AlignNumID, MreduceHID)) > 0){
            ## 最大峰组在3个以上，或者占三分之一以上
            if (length(AlignNumID) > 2 |
                length(AlignNumID)/nrow(ithAdductDF) > 1/3){
              ithAdductDF$AdductAlignmentJudge[AlignNumID] <- "AutoTrue"
            }
            ## 最大峰组少于3个，占比小于三分之一，待定
            else {
              TemporaryJudgment <- "NotSure"
            }
          }
          ## 如果最大峰组没有M+H或者M-H
          else {
            ## 最大峰组在3个以上，或者占二分之一以上
            if (length(AlignNumID) > 2 |
                length(AlignNumID)/nrow(ithAdductDF) > 1/2){
              ithAdductDF$AdductAlignmentJudge[AlignNumID] <- "AutoTrue"
            }
            ## 最大峰组少于4个，占比小于二分之一，待定
            else {
              TemporaryJudgment <- "NotSure"
            }
          }
        }
        ## 如果最多的对齐结果有多组
        else {
          ## 获得最大强度的Adduct
          AdductWithMaxIntensity <- ithAdductDF$addcut[which.max(ithAdductDF$Intensity)]
          ## 如果最大强度的Adduct是M+H或者M-H，并且强度大于其他的10倍，自动为真
          if ((ithAdductDF$Intensity[order(ithAdductDF$Intensity,decreasing = T)]/max(ithAdductDF$Intensity))[2] < 0.1 &
              (grepl("^M\\+H\\.",AdductWithMaxIntensity) | AdductWithMaxIntensity == "M+H" | grepl("^M\\-H\\.",AdductWithMaxIntensity) | AdductWithMaxIntensity == "M-H")){
            ithAdductDF$AdductAlignmentJudge[as.numeric(unlist(strsplit(ithAdductDF$AlignID[which.max(ithAdductDF$Intensity)], split = " ")))] <- "AutoTrue"
          }
          ## 如果最大强度峰不是M+H或者M-H且大于其他的10倍，执行下述判断
          else {
            ## 获取分组信息
            tempGroup <- unique(ithAdductDF$AlignID[MaxAlignNumID])
            ## 创建分组信息矩阵
            tempGroupMatrix <- as.data.frame(matrix(data = NA, nrow = length(tempGroup), ncol = 5))
            colnames(tempGroupMatrix) <- c("Group", "GroupAdduct", "M+H","M-H", "M+-H")
            tempGroupMatrix$Group <- tempGroup
            ## 循环分组
            for (ithGroup in c(1:length(tempGroup))){
              ## 获取分组的ID信息
              ithGroups <- as.numeric(unlist(strsplit(tempGroup[ithGroup], split = " ")))
              ## 获取分组的Adduct信息
              ithAdducts <- ithAdductDF$addcut[ithGroups]
              ## 生成GroupAdduct信息
              tempGroupMatrix$`GroupAdduct`[ithGroup] <- paste(sort(do.call("rbind",strsplit(ithAdducts, split = "\\."))[,1]), collapse = "|")
              ## 查看是否有M+H
              isMplusH <- length(which(grepl("^M\\+H\\.",ithAdducts) | ithAdducts == "M+H")) > 0
              tempGroupMatrix$`M+H`[ithGroup] <- isMplusH
              ## 查看是否有M-H
              isMreduceH <- length(which(grepl("^M\\-H\\.",ithAdducts) | ithAdducts == "M-H")) > 0
              tempGroupMatrix$`M-H`[ithGroup] <- isMreduceH
              ## M+H和M-H
              tempGroupMatrix$`M+-H`[ithGroup] <- isMplusH + isMreduceH
            }
            ## 如果同时有M+H和M-H的分组多于一个，自动为疑惑真
            if (length(which(tempGroupMatrix$`M+-H` == 2)) > 1){
              ithAdductDF$AdductAlignmentJudge[as.numeric(unlist(strsplit(tempGroupMatrix$Group[which(tempGroupMatrix$`M+-H` == 2)], split = " ")))] <- "DoubtfulTrue"
            }
            ## 如果同时有M+H和M-H的分组只有一个，自动为真
            else if (length(which(tempGroupMatrix$`M+-H` == 2)) == 1){
              ithAdductDF$AdductAlignmentJudge[as.numeric(unlist(strsplit(tempGroupMatrix$Group[which(tempGroupMatrix$`M+-H` == 2)], split = " ")))] <- "AutoTrue"
            }
            ## 如果同时有M+H和M-H的分组一个也没有，执行下述判断
            else {
              ## 如果有M+H或者M-H的分组大于1个，自动为疑惑真
              if (length(which(tempGroupMatrix$`M+-H` == 1)) > 1){
                ithAdductDF$AdductAlignmentJudge[as.numeric(unlist(strsplit(tempGroupMatrix$Group[which(tempGroupMatrix$`M+-H` == 1)], split = " ")))] <- "DoubtfulTrue"
              }
              ## 如果有M+H或者M-H的分组大于1个，自动为真
              else if (length(which(tempGroupMatrix$`M+-H` == 1)) == 1){
                ithAdductDF$AdductAlignmentJudge[as.numeric(unlist(strsplit(tempGroupMatrix$Group[which(tempGroupMatrix$`M+-H` == 1)], split = " ")))] <- "AutoTrue"
              }
              ## 如果有M+H或者M-H的分组一个也没有，暂时不确认
              else {
                TemporaryJudgment <- "NotSure"
              }
            }
          }
        }
        ## 对于没有确定的，进行人工判断
        if (!is.na(TemporaryJudgment)){
          ## 绘图
          print(paste0(ith,"/",length(AutomatedRtExtractorList)," | ", AutomatedRtExtractorList[[ith]]$CompoundInformation$CanonicalSMILES))
          plot(p_line)
          Sys.sleep(1)
          AdductComfirmed <- gWidgets2::ginput(msg = "Please input the ID of true adduct.",
                                               text = "1", title = "Adduct Judgment", icon = "question",
                                               parent = c(100,100), toolkit = guiToolkit())
          ## 结束命令
          if (AdductComfirmed == "e"){
            breakCommend <- "Break"
            break
          }
          ## 有问号和没问号
          else if (grepl("\\?", AdductComfirmed) | grepl("？", AdductComfirmed)){
            ithAdductDF$AdductAlignmentJudge[as.vector(unique(as.numeric(unlist(strsplit(ithAdductDF$AlignID[as.vector(as.numeric(unlist(strsplit(AdductComfirmed, split = " "))[-1]))], split = " ")))))] <- "DoubtfulTrue"
          }
          ## 没有任何一个是对的
          else if (AdductComfirmed == "0"){
            ithAdductDF$AdductAlignmentJudge <- "ManualFalse"
          }
          ## 输入不带问号的数字的时候
          else {
            ## 获取原始的输入的ID
            RawAdductComfirmedID <- as.vector(as.numeric(unlist(strsplit(AdductComfirmed, split = " "))))
            ## 获取所有的ID
            AdductComfirmedID <- as.vector(unique(as.numeric(unlist(strsplit(ithAdductDF$AlignID[RawAdductComfirmedID], split = " ")))))
            ## 如果输入的是一个分组里的，为手动真
            if (length(unique(ithAdductDF$AlignID[RawAdductComfirmedID])) == 1){
              ithAdductDF$AdductAlignmentJudge[as.vector(unique(as.numeric(unlist(strsplit(ithAdductDF$AlignID[as.vector(as.numeric(unlist(strsplit(AdductComfirmed, split = " "))))], split = " ")))))] <- "ManualTrue"
            }
            ## 如果输入的是多个分组里的，为手动疑惑真
            else {
              ithAdductDF$AdductAlignmentJudge[as.vector(unique(as.numeric(unlist(strsplit(ithAdductDF$AlignID[as.vector(as.numeric(unlist(strsplit(AdductComfirmed, split = " "))))], split = " ")))))] <- "ManualDoubtfulTrue"
            }

          }
        }
      }
      ## 判断部分---------------------------------------------------------------
    }
    AutomatedRtExtractorList[[ith]]$AdductAlignmentResult <- ithAdductDF
  }
  return(AutomatedRtExtractorList)
}
