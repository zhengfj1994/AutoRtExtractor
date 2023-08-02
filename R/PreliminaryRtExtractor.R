## PreliminaryRtExtractor defined by Zheng Fujian
#' @title Extract retention time from mzML data
#' @author Fujian Zheng
PreliminaryRtExtractor <- function(adductFile = NA, mzMLFile, MixedStandards, ppm = 10, haldWidth = 20, ignore_threshold = 0.1, span = 9){
  library(xcms)
  library(MALDIquant)
  library(MALDIquantForeign)
  library(magrittr)
  library(scorepeak)
  library(OrgMassSpecR)
  library(enviPat)
  library(MetEx)
  library(ggpmisc)
  if (is.na(adductFile)){
    data(AREadducts)
  } else {
    adducts <- read.csv(adductFile)
  }
  data(isotopes) # enviPat

  ## Read mzML file using xcmaRaw
  mxMLRawData <- xcms::xcmsRaw(mzMLFile)

  ## Get the mzML file name
  mzMLFileName <- rev(unlist(strsplit(mzMLFile, split = "/")))[1]
  mzMLFileName <- unlist(strsplit(mzMLFileName, split = ".mz"))[1]

  ## Get the mix number

  mzMLFileName

  # ## Read the information of chemical standards.
  # MixedStandards <- read.csv(file = MixedStandardsFile, fileEncoding = 'GBK')
  ## Screen the chemical standards.
  MixedStandards_1 <- MixedStandards[which(MixedStandards$Group.Name == unlist(strsplit(mzMLFileName, split = "_|-"))[3] | MixedStandards$Group.Name == "IS"),]
  ## Create a list for saving results.
  AutomatedRtExtractorList <- vector("list", length = nrow(MixedStandards_1))
  names(AutomatedRtExtractorList) <- paste0(unlist(strsplit(mzMLFileName, split = "_|-"))[3], "_", c(1:nrow(MixedStandards_1)))

  ## 将化合物信息存到list里.
  for (ithRow in c(1:nrow(MixedStandards_1))){
    AutomatedRtExtractorList[[ithRow]]$CompoundInformation <- MixedStandards_1[ithRow,]
  }

  ## 获取每个加和形式的TIC并保存于list中
  for (ithRow in c(1:length(AutomatedRtExtractorList))){
    # Get the MF and MW in table
    Charge <- AutomatedRtExtractorList[[ithRow]]$CompoundInformation$Charge
    if (Charge == 0){
      MF <- AutomatedRtExtractorList[[ithRow]]$CompoundInformation$MolecularFormula
      MF <- check_chemform(isotopes,MF)$new_formula
      MW <- check_chemform(isotopes,MF)$monoisotopic_mass
    }else {
      MF <- AutomatedRtExtractorList[[ithRow]]$CompoundInformation$MolecularFormula
      MF <- substr(MF,1,nchar(MF)-abs(Charge))
      MF <- check_chemform(isotopes,MF)$new_formula
      MW <- check_chemform(isotopes,MF)$monoisotopic_mass
    }

    # Calculate m/z of different adduct
    ithAdducts <- adducts
    if (Charge == 0){
      ithAdducts <- ithAdducts[which(ithAdducts$Formula_add != "FALSE" | ithAdducts$Formula_ded != "FALSE"),]
    } else if (Charge == 1){
      ithAdducts <- ithAdducts[which(ithAdducts$Name == "M+"),]
    } else if (Charge == 2){
      ithAdducts <- ithAdducts[which(ithAdducts$Name == "M++"),]
    } else if (Charge == -1){
      ithAdducts <- ithAdducts[which(ithAdducts$Name == "M-"),]
    }
    ithAdducts$mzFromcalc <- MW/abs(ithAdducts$Charge)*ithAdducts$Mult + ithAdducts$Mass
    ithAdducts$Formula <- ithAdducts$MF
    for (i in c(1:nrow(ithAdducts))){
      if (ithAdducts$Formula_add[i] != "FALSE" & ithAdducts$Formula_ded[i] != "FALSE"){
        ithAdducts$AdductFormula[i] <- enviPat::subform(
          enviPat::mergeform(
            enviPat::multiform(MF,ithAdducts$Mult[i]),
            ithAdducts$Formula_add[i]),
          ithAdducts$Formula_ded[i])
      } else if (ithAdducts$Formula_add[i] != "FALSE" & ithAdducts$Formula_ded[i] == "FALSE"){
        ithAdducts$AdductFormula[i] <- enviPat::mergeform(
          enviPat::multiform(MF,ithAdducts$Mult[i]),ithAdducts$Formula_add[i])
      } else if (ithAdducts$Formula_add[i] == "FALSE" & ithAdducts$Formula_ded[i] != "FALSE"){
        ithAdducts$AdductFormula[i] <- enviPat::subform(
          enviPat::multiform(MF,ithAdducts$Mult[i]),ithAdducts$Formula_ded[i])
      }else{
        ithAdducts$AdductFormula[i] <- MF
      }
    }
    ithAdducts <- ithAdducts[which(!grepl("formula", ithAdducts$AdductFormula)),]
    ithAdducts$MonoisotopicMass <- check_chemform(isotopes,ithAdducts$AdductFormula)$monoisotopic_mass
    ithAdducts$mzFromFormula <- (ithAdducts$MonoisotopicMass - ithAdducts$Charge * 0.00054857990924)/abs(ithAdducts$Charge)
    ithAdducts$mzmin <- ithAdducts$mzFromFormula - ithAdducts$mzFromFormula*ppm/2000000
    ithAdducts$mzmax <- ithAdducts$mzFromFormula + ithAdducts$mzFromFormula*ppm/2000000
    if (unlist(strsplit(mzMLFileName, split = "_|-"))[2] == "POS"){
      ithAdducts <- ithAdducts[which(ithAdducts$Ion_mode == "positive"),]
    } else {
      ithAdducts <- ithAdducts[which(ithAdducts$Ion_mode == "negative"),]
    }

    if (nrow(ithAdducts) == 0){
      AutomatedRtExtractorList[[ithRow]]$Addcuts <- NA
    } else {
      #  Use the function named "rawEIC" in xcms to get the EIC of specific m/z
      batchRawEIC <- function(mzmin, mzmax){
        EICdata <- xcms::rawEIC(mxMLRawData,
                                mzrange = c(mzmin,mzmax))
        EICdata$scantime <- mxMLRawData@scantime
        return(EICdata)
      }
      mzminList <- split(ithAdducts$mzmin, 1:nrow(ithAdducts))
      mzmaxList <- split(ithAdducts$mzmax, 1:nrow(ithAdducts))
      TIClist <- mapply(batchRawEIC, mzminList, mzmaxList, SIMPLIFY = F)
      names(TIClist) <- ithAdducts$Name

      # Combind the mz and intensity to a data.frame
      TIClistCombiner <- function(ithTIClist){
        TICmatrix <- as.data.frame(cbind(ithTIClist$scan, ithTIClist$scantime, ithTIClist$intensity))
        colnames(TICmatrix) <- c("scan", "scantime", "intensity")
        return(TICmatrix)
      }
      TIClist <- lapply(TIClist, TIClistCombiner)

      AdductList <- vector("list", length = length(TIClist))
      names(AdductList) <- names(TIClist)
      for (i in c(1:length(AdductList))){
        AdductList[[i]]$AdductInformation <- ithAdducts[i,]
        AdductList[[i]]$TIC <- TIClist[[i]]
      }
      AutomatedRtExtractorList[[ithRow]]$Addcuts <- AdductList
    }
  }

  ## 对TIC进行峰检测获取信息后保存于list中
  for (ithRow in c(1:length(AutomatedRtExtractorList))){
    if (!is.na(AutomatedRtExtractorList[[ithRow]]$Addcuts[1])){
      for (jthAdduct in c(1:length(AutomatedRtExtractorList[[ithRow]]$Addcuts))){
        ## 获取指定的TIC
        jthTIC <- AutomatedRtExtractorList[[ithRow]]$Addcuts[[jthAdduct]]$TIC
        # 获取TIC的强度列
        intForPeakDetection <- jthTIC$intensity
        # ## 依据scorepeak中的detect_localmaxima进行峰检测
        # detectedPeaks <- scorepeak::detect_localmaxima(intForPeakDetection, 9, boundary = "discard")
        # ## 对检测到的峰进行打分
        # detectedPeaksScore <- score_type2(intForPeakDetection, 51)
        # ## 筛选打分大于1的峰
        # ## turePeaks <- detectedPeaksScore > 0.3 & detectedPeaks
        # turePeaks <- Reduce(intersect,
        #                     list(which(detectedPeaksScore > 1) + 1,
        #                          which(detectedPeaksScore > 1),
        #                          which(detectedPeaksScore > 1) - 1,
        #                          which(detectedPeaks)))
        turePeaks <- which(ggpmisc:::find_peaks(intForPeakDetection, ignore_threshold = ignore_threshold, span = span, strict = F))
        tureValleys <- which(ggpmisc:::find_peaks(-intForPeakDetection, ignore_threshold = ignore_threshold, span = span, strict = F))
        LeftScan <- turePeaks
        RightScan <- turePeaks
        if (length(turePeaks) > 0){
          for (ithTurePeaks in c(1:length(turePeaks))){
            LeftScan[ithTurePeaks] <- rev(tureValleys[which(tureValleys - turePeaks[ithTurePeaks] < 0)])[1]
            RightScan[ithTurePeaks] <- tureValleys[which(tureValleys - turePeaks[ithTurePeaks] > 0)][1]
          }
        }

        ## 构建检测到的峰的矩阵
        DetectedPeaks <- as.data.frame(cbind(turePeaks,
                                             mxMLRawData@scantime[turePeaks],
                                             intForPeakDetection[turePeaks],
                                             LeftScan,
                                             mxMLRawData@scantime[LeftScan],
                                             intForPeakDetection[LeftScan],
                                             RightScan,
                                             mxMLRawData@scantime[RightScan],
                                             intForPeakDetection[RightScan]))
        colnames(DetectedPeaks) <- c("Scan","RetentionTime","Intensity","LeftScan","LeftRetentionTime", "LeftIntensity", "RightScan", "RightRetentionTime", "RightIntensity")
        ## 按强度从大到小排序
        DetectedPeaks <- DetectedPeaks[order(DetectedPeaks$Intensity, decreasing = T),]

        ## 将峰矩阵存入列表
        AutomatedRtExtractorList[[ithRow]]$Addcuts[[jthAdduct]]$DetectedPeaks <- DetectedPeaks
      }
    }
  }

  ## 获取检测到的峰的多维信息
  for (ithRow in c(1:length(AutomatedRtExtractorList))){
    if (!is.na(AutomatedRtExtractorList[[ithRow]]$Addcuts[1])){
      for (jthAdduct in c(1:length(AutomatedRtExtractorList[[ithRow]]$Addcuts))){
        ## 获取AdductFormula
        AdductFormula <- AutomatedRtExtractorList[[ithRow]]$Addcuts[[jthAdduct]]$AdductInformation$AdductFormula
        ## 获取m/z
        mzFromFormula <- AutomatedRtExtractorList[[ithRow]]$Addcuts[[jthAdduct]]$AdductInformation$mzFromFormula
        ## 获取Charge
        Charge <- AutomatedRtExtractorList[[ithRow]]$Addcuts[[jthAdduct]]$AdductInformation$Charge
        ## 获取TIC
        TIC <- AutomatedRtExtractorList[[ithRow]]$Addcuts[[jthAdduct]]$TIC
        ## 获取检测到的峰
        DetectedPeaks <- AutomatedRtExtractorList[[ithRow]]$Addcuts[[jthAdduct]]$DetectedPeaks
        ## 获取检测到的峰的附近的EIC
        DetectedPeaks$StartRt <- DetectedPeaks$RetentionTime - haldWidth
        DetectedPeaks$EndRt <- DetectedPeaks$RetentionTime + haldWidth
        ## 如果没有检测到峰，所有信息都设置为空
        if (nrow(DetectedPeaks) == 0){
          DetectedPeaks$EIC <- DetectedPeaks$Scan
          DetectedPeaks$IE <- DetectedPeaks$Scan
          DetectedPeaks$amplitude <- DetectedPeaks$Scan
          DetectedPeaks$center <- DetectedPeaks$Scan
          DetectedPeaks$sigma <- DetectedPeaks$Scan
          DetectedPeaks$gamma <- DetectedPeaks$Scan
          DetectedPeaks$fitStatus <- DetectedPeaks$Scan
          DetectedPeaks$curveModel <- DetectedPeaks$Scan
          DetectedPeaks$predicttedCurve <- DetectedPeaks$Scan
          DetectedPeaks$residual <- DetectedPeaks$Scan
          DetectedPeaks$MS1data <- DetectedPeaks$Scan
          DetectedPeaks$experimentalIsotopePattern <- DetectedPeaks$Scan
          DetectedPeaks$theoreticallyIsotopePattern <- DetectedPeaks$Scan
          DetectedPeaks$isotopeSimilarity <- DetectedPeaks$Scan
          DetectedPeaks$isotopeAlignment <- DetectedPeaks$Scan
        }else { ## 如果检测到峰，获取多维信息
          for (i in c(1:nrow(DetectedPeaks))){
            ## 获取EIC
            ithEIC <- TIC[which(TIC$scantime - DetectedPeaks$StartRt[i] >= 0)[1]:rev(which((TIC$scantime-DetectedPeaks$EndRt[i]) <= 0))[1],]
            ithEIC <- na.omit(ithEIC)
            ## 写入EIC
            DetectedPeaks$EIC[i] <- paste(paste(ithEIC$scan,ithEIC$scantime,ithEIC$intensity, sep = " "), collapse = ";")
            ## 写入IE
            DetectedPeaks$IE[i] <- as.numeric(MetEx::entropyCalculator(ithEIC[,c(1,3)])[1,2])
            ## 使用peakPantheR中的fitCurve模拟峰形
            fittedCurve <- fitCurve(ithEIC$scantime, ithEIC$intensity, curveModel = "emgGaussian")
            ## 获取fitCurve的多维参数
            DetectedPeaks$amplitude[i] <- fittedCurve$amplitude
            DetectedPeaks$center[i] <- fittedCurve$center
            DetectedPeaks$sigma[i] <- fittedCurve$sigma
            DetectedPeaks$gamma[i] <- fittedCurve$gamma
            DetectedPeaks$fitStatus[i] <- fittedCurve$fitStatus
            DetectedPeaks$curveModel[i] <- fittedCurve$curveModel
            ## 根据fitCurve预测峰
            predicttedCurve <- as.data.frame(cbind(ithEIC$scantime,predictCurve(fittedCurve,ithEIC$scantime)))
            colnames(predicttedCurve) <- c("scantime", "intensity")
            ## 预测峰写入列表
            DetectedPeaks$predicttedCurve[i] <- paste(paste(predicttedCurve$scantime,predicttedCurve$intensity, sep = " "), collapse = ";")
            ## 计算实际峰和预测峰的残差并写入列表
            DetectedPeaks$residual[i] <- sum(abs(ithEIC$intensity - predicttedCurve$intensity))/sum(ithEIC$intensity)
            ## 根据峰顶点的Scan点获取MS1谱图
            MS1Data <- xcms::getScan(mxMLRawData, scan = DetectedPeaks$Scan[i])
            ## 只保留指定质荷比前后的MS1谱图
            MS1Data <- as.matrix(MS1Data[which(MS1Data[,1] > mzFromFormula - 0.5 & MS1Data[,1] < mzFromFormula + 5.5),])
            ## 如果MS1Data只有一列，需要进行转置
            if (ncol(MS1Data) == 1){
              MS1Data <- t(MS1Data)
            }
            ## 将MS1写入list
            DetectedPeaks$MS1data[i] <- paste(paste(MS1Data[,1],MS1Data[,2], sep = " "), collapse = ";")
            ## Create MassSpectrum by using "createMassSpectrum" in MALDIquant
            MassSpectrum <- MALDIquant::createMassSpectrum(mass = MS1Data[,1], intensity = MS1Data[,2], metaData=list(name="example spectrum"))
            ## Detect peaks in isotope pattern
            peaks = tryCatch({
              MALDIquant::detectPeaks(MassSpectrum, method="MAD",halfWindowSize=3,SNR=0.1)
            }, warning = function(w) {
              MALDIquant::detectPeaks(MassSpectrum, method="MAD",halfWindowSize=3,SNR=0.1)
            }, error = function(e) {
              createMassPeaks(mass=MS1Data[,1], intensity=MS1Data[,2],
                              metaData=list(name="peaks"))
            }, finally = {})
            ## Create a matrix that save mz, intensity and snr
            experimentalIsotopePattern <- as.data.frame(cbind(peaks@mass,peaks@intensity,peaks@snr))
            colnames(experimentalIsotopePattern) <- c("mass", "intensity", "snr")

            ## 保存实验同位素分布到DetectedPeaks
            DetectedPeaks$experimentalIsotopePattern[i] <- paste(paste(experimentalIsotopePattern$mass,experimentalIsotopePattern$intensity,experimentalIsotopePattern$snr, sep = " "), collapse = ";")
            # Calculate the theoretically isotope pattern
            pattern <- enviPat::isopattern(isotopes,
                                           chemforms = AdductFormula,
                                           threshold = 0.1,
                                           plotit = FALSE,
                                           charge = Charge,
                                           emass = 0.00054858,
                                           algo = 2,
                                           verbose = F)
            profiles <- enviPat::envelope(pattern,
                                          ppm = FALSE,
                                          dmz = 0.0001,
                                          frac = 1/4,
                                          env = "Gaussian",
                                          resolution = 30000,
                                          plotit = FALSE,
                                          verbose = F)
            theoreticallyIsotopePattern <- vdetect(profiles,
                                                   detect = "centroid",
                                                   plotit = F,
                                                   verbose = F)
            theoreticallyIsotopePattern <- as.data.frame(theoreticallyIsotopePattern[[1]])
            colnames(theoreticallyIsotopePattern) <- c("mass", "intensity")

            ## 写入理论同位素分布到DetectedPeaks
            DetectedPeaks$theoreticallyIsotopePattern[i] <- paste(paste(theoreticallyIsotopePattern$mass,theoreticallyIsotopePattern$intensity, sep = " "), collapse = ";")
            ## 归一化实验同位素分布
            experimentalIsotopePattern$intensity <- experimentalIsotopePattern$intensity/max(experimentalIsotopePattern$intensity) * 100

            ## 如果最大理论同位素分布小于最小实验同位素分布，同位素相似性没有
            if (max(theoreticallyIsotopePattern$mass) < min(experimentalIsotopePattern$mass)){
              DetectedPeaks$isotopeSimilarity[i] <- NA
              DetectedPeaks$isotopeAlignment[i] <- NA
            } else{ ## 计算同位素相似性
              isotopeSimilarity <- SpectrumSimilarity(experimentalIsotopePattern[,c(1,2)], theoreticallyIsotopePattern,
                                                      top.label = "Experimental Isotope Pattern",
                                                      bottom.label = "Theoretically Isotope Pattern",
                                                      xlim = c(min(theoreticallyIsotopePattern$mass)-1, max(theoreticallyIsotopePattern$mass)+1),
                                                      output.list = TRUE,
                                                      print.graphic = F)
              ## 同位素相似性写入DetectedPeaks
              DetectedPeaks$isotopeSimilarity[i] <- isotopeSimilarity$similarity.score
              ## 同位素对齐结果写入DetectedPeaks
              DetectedPeaks$isotopeAlignment[i] <- paste(paste(isotopeSimilarity$alignment$mz,isotopeSimilarity$alignment$intensity.top,isotopeSimilarity$alignment$intensity.bottom, sep = " "), collapse = ";")
            }
          }
        }
        ## 将DetectedPeaks写入list
        DetectedPeaks[order(DetectedPeaks$Intensity, decreasing = T),] -> AutomatedRtExtractorList[[ithRow]]$Addcuts[[jthAdduct]]$DetectedPeaks
      }
    }
  }
  return(AutomatedRtExtractorList)
}
##------------------------------------------------------------------------------
