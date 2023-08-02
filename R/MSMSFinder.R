MSMSFinder <- function(AutomatedRtExtractorList, mgfFilePath, ppm = 20, deltaRt = 6, MSMSthreshold = 0.001, cores = 2){
  library(doSNOW)
  library(progress)
  # mgfFilePath <- "U:/zfjDB/Data/mgf/CS26"
  ## 定义mgfReader函数用于读取mgf文件
  mgfReader <- function(mgfFile){
    mgfFileContent <- read.csv2(file = mgfFile, header = F)
    beginIndex <- which(mgfFileContent == "BEGIN IONS")
    endIndex <- which(mgfFileContent == "END IONS")
    RTIndex <- which(grepl("RTINSECONDS=",mgfFileContent$V1))
    pepmassIndex <- which(grepl("PEPMASS=",mgfFileContent$V1))
    mgfMatrix <- as.data.frame(matrix(data = NA, nrow = length(beginIndex), ncol = 4))
    colnames(mgfMatrix) <- c("FileName", "RT", "PEPMASS", "MSMS")
    mgfMatrix$FileName <- unlist(strsplit(rev(unlist(strsplit(mgfFile, split = "/")))[1], split = ".mgf"))[1]
    for (i in c(1:length(beginIndex))){
      mgfMatrix$RT[i] <- strsplit(mgfFileContent$V1[RTIndex[i]], split = "=")[[1]][2]
      mgfMatrix$PEPMASS[i] <- strsplit(strsplit(mgfFileContent$V1[pepmassIndex[i]], split = "=")[[1]][2], split = " ")[[1]][1]

      ithMSMS <- mgfFileContent$V1[(beginIndex[i]):(endIndex[i])]
      ithMSMS <- ithMSMS[!grepl("[a-zA-Z]", ithMSMS)]
      if (length(ithMSMS) == 0){
        mgfMatrix$MSMS[i] <- NA
      } else {
        if (length(unlist(strsplit(ithMSMS[1], split = " "))) == 2){
          ithMSMS <- strsplit(ithMSMS, split = " ")
          ithMSMS <- data.frame(matrix(unlist(ithMSMS), byrow = T, ncol = 2), stringsAsFactors = F)
          colnames(ithMSMS) <- c("mz", "intensity")
          ithMSMS$mz <- as.numeric(ithMSMS$mz)
          ithMSMS$intensity <- as.numeric(ithMSMS$intensity)
          ithMSMS <- ithMSMS[which(ithMSMS$intensity/max(ithMSMS$intensity) > MSMSthreshold),]
          mgfMatrix$MSMS[i] <- paste(paste(ithMSMS$mz,ithMSMS$intensity, sep = " "), collapse = ";")
        } else if (length(unlist(strsplit(ithMSMS[1], split = " "))) == 3){
          ithMSMS <- strsplit(ithMSMS, split = " ")
          ithMSMS <- data.frame(matrix(unlist(ithMSMS), byrow = T, ncol = 3), stringsAsFactors = F)
          colnames(ithMSMS) <- c("mz", "intensity", "charge")
          ithMSMS$mz <- as.numeric(ithMSMS$mz)
          ithMSMS$intensity <- as.numeric(ithMSMS$intensity)
          ithMSMS <- ithMSMS[which(ithMSMS$intensity/max(ithMSMS$intensity) > MSMSthreshold),]
          mgfMatrix$MSMS[i] <- paste(paste(ithMSMS$mz,ithMSMS$intensity, sep = " "), collapse = ";")
        } else {
          packageStartupMessage("The format of mgf is wrong!")
        }
      }
      # if (beginIndex[i] + 4 == endIndex[i]){
      #   mgfMatrix$MSMS[i] <- NA
      # }
      # else {
      #   ithMSMS <- mgfFileContent$V1[(beginIndex[i] + 4):(endIndex[i] - 1)]
      #   ithMSMS <- strsplit(ithMSMS, split = " ")
      #   ithMSMS <- data.frame(matrix(unlist(ithMSMS), byrow = T, ncol = 2), stringsAsFactors = F)
      #   colnames(ithMSMS) <- c("mz", "intensity")
      #   ithMSMS$intensity <- as.numeric(ithMSMS$intensity)
      #   ithMSMS <- ithMSMS[which(ithMSMS$intensity/max(ithMSMS$intensity) > 0.005),]
      #   mgfMatrix$MSMS[i] <- paste(paste(ithMSMS$mz,ithMSMS$intensity, sep = " "), collapse = ";")
      # }
    }
    mgfMatrix <- mgfMatrix[which(!is.na(mgfMatrix$MSMS) & mgfMatrix$MSMS != ""),]
    return(mgfMatrix)
  }
  ## 处理mgf文件，生成MSMS信息矩阵
  ## 获取所有的mgf文件的全名
  mgfFiles <- list.files(mgfFilePath, pattern = "mgf", recursive = T, full.names = T)
  ## 读取mgf文件并存放于矩阵中
  cl <- makeSOCKcluster(cores)
  registerDoSNOW(cl)
  # progress bar ------------------------------------------------------------
  iterations <- length(mgfFiles)
  pb <- progress_bar$new(
    format = ":letter [:bar] :elapsed | Remaining time: :eta <br>",
    total = iterations,
    width = 120)
  # allowing progress bar to be used in foreach -----------------------------
  progress <- function(n){
    pb$tick(tokens = list(letter = "Progress of reading mgf files."))
  }
  opts <- list(progress = progress)
  ## 读取所有mgf文件的命令
  mgfMatrix <- foreach(i=mgfFiles, .options.snow=opts, .combine='rbind') %dopar% mgfReader(i)
  ## 将mgfMatrix中的第一列分列
  fileNameDf <- as.data.frame(do.call("rbind", strsplit(mgfMatrix$FileName, split = "_")))
  colnames(fileNameDf) <- c("CSname","IonMode","Mix")
  ## 重整mgfMatrix
  mgfMatrix <- cbind(fileNameDf, mgfMatrix[,-1])

  ## 获取AutomatedRtExtractorList的名字
  AutomatedRtExtractorListNames <- names(AutomatedRtExtractorList)
  ## 创建进度条
  pb <- txtProgressBar(style=3)
  ## 循环AutomatedRtExtractorList
  for (ith in c(1:length(AutomatedRtExtractorList))){
    setTxtProgressBar(pb, ith/length(AutomatedRtExtractorList))
    ## 获取第i个AutomatedRtExtractorListNames
    ithAutomatedRtExtractorListName <- AutomatedRtExtractorListNames[ith]
    ithAutomatedRtExtractorListName <- unlist(strsplit(ithAutomatedRtExtractorListName,split = "_"))[1]
    ## 获取第i个AutomatedRtExtractorList
    ithAutomatedRtExtractorList <- AutomatedRtExtractorList[[ith]]
    ## 获取AdductAlignmentResult
    ithAdductAlignmentResult <- ithAutomatedRtExtractorList$AdductAlignmentResult
    ## 在ithAdductAlignmentResult最前面新增m/z列,charge列和MSMS列
    newMzColumn <- as.data.frame(matrix(data = NA, nrow = nrow(ithAdductAlignmentResult), ncol = 1))
    colnames(newMzColumn) <- "mz"
    newChargeColunm <- as.data.frame(matrix(data = NA, nrow = nrow(ithAdductAlignmentResult), ncol = 1))
    colnames(newChargeColunm) <- "charge"
    newMSMSColunm <- as.data.frame(matrix(data = NA, nrow = nrow(ithAdductAlignmentResult), ncol = 1))
    colnames(newMSMSColunm) <- "MSMS"
    ithAdductAlignmentResult <- cbind(newMzColumn,newChargeColunm,ithAdductAlignmentResult,newMSMSColunm)
    ## 获取Adduct
    ithAddcuts <- ithAutomatedRtExtractorList$Addcuts

    ## 找到每个Adduct的m/z和charge
    ## 如果ithAdductAlignmentResult为0行，跳过
    if (nrow(ithAdductAlignmentResult) == 0){
      ## 删除无关信息节省空间
      AutomatedRtExtractorList[[ith]]$Addcuts <- NULL
      AutomatedRtExtractorList[[ith]]$AddcutsRT <- NULL
      AutomatedRtExtractorList[[ith]]$AdductAlignmentResult <- NULL
      ## 写入结果
      AutomatedRtExtractorList[[ith]]$MSMSmatchedResult <- ithAdductAlignmentResult
      next()
    }
    ## 如果ithAdductAlignmentResult不为0行，执行下述命令
    for (jth in c(1:nrow(ithAdductAlignmentResult))){
      ## 获取第j个Adduct的名字
      jthAdductName <- ithAdductAlignmentResult$addcut[jth]
      ## 只保留Adduct，去除其编号
      jthAdductName <- unlist(strsplit(jthAdductName, split = "\\."))[1]
      ## 根据Adduct的名字，获取其相关信息
      jthAdductInformation <- ithAddcuts[[jthAdductName]]$AdductInformation
      ## 获取m/z
      jthMz <- jthAdductInformation$mzFromFormula[1]
      ## 将mz写入到ithAdductAlignmentResult
      ithAdductAlignmentResult$mz[jth] <- jthMz
      ## 将charge写入到ithAdductAlignmentResult
      ithAdductAlignmentResult$charge[jth] <- jthAdductInformation$Charge[1]
    }

    ## 为每个Adduct匹配二级质谱
    ## 如果ithAdductAlignmentResult为0行，跳过
    if (nrow(ithAdductAlignmentResult) == 0){
      ## 删除无关信息节省空间
      AutomatedRtExtractorList[[ith]]$Addcuts <- NULL
      AutomatedRtExtractorList[[ith]]$AddcutsRT <- NULL
      AutomatedRtExtractorList[[ith]]$AdductAlignmentResult <- NULL
      ## 写入结果
      AutomatedRtExtractorList[[ith]]$MSMSmatchedResult <- ithAdductAlignmentResult
      next()
    }
    ## 如果ithAdductAlignmentResult不为0行，执行下述命令
    else {
      ## 获取判断为真的行的ID
      tureID <- grep("True",ithAdductAlignmentResult$AdductAlignmentJudge)
      ## 如果没有判断为真的结果，跳过
      if (length(tureID) == 0){
        next()
      }
      ## 如果有判断为真的结果，执行下述命令
      else {
        ## 循环tureID
        for (kthTrue in tureID){
          ## 获取mz
          kthMz <- ithAdductAlignmentResult$mz[kthTrue]
          ## 获取Rt
          kthRt <- ithAdductAlignmentResult$RT[kthTrue]
          ## 获取Charge
          kthCharge <- ithAdductAlignmentResult$charge[kthTrue]
          ## 获取IonMode
          if (kthCharge > 0){
            kthIonMode <- "POS"
          }
          else {
            kthIonMode <- "NEG"
          }
          ## 根据离子模式，混标分组，质荷比和保留时间进行匹配
          matchID <- which(mgfMatrix$IonMode == kthIonMode &
                  mgfMatrix$Mix == ithAutomatedRtExtractorListName &
                  abs(as.numeric(mgfMatrix$PEPMASS) - kthMz)/kthMz*1000000 < ppm &
                  abs(as.numeric(mgfMatrix$RT) - kthRt) < deltaRt)
          ## 如果没有匹配结果，MSMS为NA
          if (length(matchID) == 0){
            ithAdductAlignmentResult$MSMS[kthTrue] <- NA
          }
          ## 如果有匹配结果，MSMS写入
          else {
            ithAdductAlignmentResult$MSMS[kthTrue] <- paste(mgfMatrix$MSMS[matchID], collapse = "|||")
          }
        }
      }
    }
    ## 删除无关信息节省空间
    AutomatedRtExtractorList[[ith]]$Addcuts <- NULL
    AutomatedRtExtractorList[[ith]]$AddcutsRT <- NULL
    AutomatedRtExtractorList[[ith]]$AdductAlignmentResult <- NULL
    ## 写入结果
    AutomatedRtExtractorList[[ith]]$MSMSmatchedResult <- ithAdductAlignmentResult
  }
  ## 关闭进度条
  close(pb)
  return(AutomatedRtExtractorList)
}
