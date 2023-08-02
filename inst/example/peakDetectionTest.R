library(scorepeak)
library(ggplot2)
library(ggpmisc)
for (i in c(1:1000)){
  print(i)
  TIC <- AutomatedRtExtractorList[[i]]$Addcuts[[1]]$TIC[,c(1,3)]
  # TIC <- TIC[0:400,]
  intForPeakDetection <- TIC$intensity
  # ## 依据scorepeak中的detect_localmaxima进行峰检测
  # detectedPeaks <- scorepeak::detect_localmaxima(intForPeakDetection, 9, boundary = "discard")
  # ## 对检测到的峰进行打分
  # detectedPeaksScore <- score_type2(intForPeakDetection, 51)
  # ## 筛选打分大于1的峰
  # turePeaks <- Reduce(intersect,
  #                     list(which(detectedPeaksScore > 1) + 1,
  #                          which(detectedPeaksScore > 1),
  #                          which(detectedPeaksScore > 1) - 1,
  #                          which(detectedPeaks)))
  turePeaks <- ggpmisc:::find_peaks(intForPeakDetection, ignore_threshold = 0.1, span = 5, strict = F)
  turePeaks <- ggpmisc:::find_peaks(-intForPeakDetection, ignore_threshold = 0.1, span = 5, strict = F)

  # p0 <- ggplot(TIC, aes(scan, intensity)) + geom_line() +
  #   stat_peaks(colour = "red") +
  #   stat_valleys(colour = "blue")

  p0 <- ggplot() +
    geom_line(data = TIC, aes(x=scan,y=intensity)) +
    geom_point(data = TIC[turePeaks,], aes(x = scan, y = intensity), color = "red")
  plot(p0)
  Sys.sleep(2)
}

#
# local_peaks_screened <- detect_localmaxima(intForPeakDetection, 13)
# score <- score_type2(intForPeakDetection, 51)
# plot(intForPeakDetection, type = "l", main = "Screened Local Peaks")
# points(which(local_peaks), intForPeakDetection[local_peaks], col = "red", pch = 19)
# points(seq(length(score)), score, type = "l", col = "blue")
