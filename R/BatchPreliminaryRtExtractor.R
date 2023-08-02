BatchPreliminaryRtExtractor <- function(mzMLFolder, CSprefix = "CS", MixedStandardsFile, ResultPath, ppm = 10, haldWidth = 20, ignore_threshold = 0.1, span = 9){
  AllmzMLFiles_FullName <- list.files(path = mzMLFolder, full.names = T)
  AllmzMLFiles <- list.files(path = mzMLFolder)
  CSname <- unique(do.call("rbind",strsplit(AllmzMLFiles,split = "_|-"))[,1])
  CSname <- CSname[order(as.numeric(do.call("rbind",strsplit(CSname, CSprefix))[,2]))]

  ## Read the information of chemical standards.
  MixedStandards <- read.csv(file = MixedStandardsFile, fileEncoding = 'GBK')

  for (ithCS in CSname){
    ithCS_mzMLFiles_FullName <- AllmzMLFiles_FullName[grep(paste0(ithCS, "_"), AllmzMLFiles)]
    ithCS_mzMLFiles <- AllmzMLFiles[grep(paste0(ithCS, "_"), AllmzMLFiles)]

    sortMatrix <- as.data.frame(do.call("rbind", strsplit(unlist(strsplit(ithCS_mzMLFiles, split = ".mzML")), split = "_")))
    sortMatrix[,3] <- as.numeric(do.call("rbind", strsplit(sortMatrix[,3],split = "Mix"))[,2])
    ithCS_mzMLFiles_FullName <- ithCS_mzMLFiles_FullName[order(sortMatrix[,2],sortMatrix[,3])]
    ithCS_mzMLFiles_FullName_POS <- ithCS_mzMLFiles_FullName[grep("_POS_", ithCS_mzMLFiles_FullName)]
    ithCS_mzMLFiles_FullName_NEG <- ithCS_mzMLFiles_FullName[grep("_NEG_", ithCS_mzMLFiles_FullName)]

    ## mzML file check
    packageStartupMessage("mzML file checking")
    for (mzMLFile in ithCS_mzMLFiles_FullName){
      print(mzMLFile)
      tryCatch({
        xcms::xcmsRaw(mzMLFile)
      }, warning = function(w){
        # 这里是出现warning状态时，应该怎么做，可以用print打印出来，可以执行其它命令
        packageStartupMessage(paste0("Warning! ", mzMLFile))
      }, error = function(e){
        # 这里时出现Error状态时，应该怎么做，可以用print打印出来，也可以执行其它命令
        packageStartupMessage(paste0("Error! ", mzMLFile))
      },finally = {
        # 这是运行正常时，应该怎么做，可以用print打印出来，也可以执行其它命令
        packageStartupMessage(paste0("Right! ", mzMLFile))
      })
    }
    packageStartupMessage("mzML file check finised")
    ## mzML file check

    AutomatedRtExtractorList <- list()
    for (mzMLFile in ithCS_mzMLFiles_FullName_POS){
      print(mzMLFile)
      ithAutomatedRtExtractorList_POS <- PreliminaryRtExtractor(adductFile = NA,
                                         mzMLFile = mzMLFile,
                                         MixedStandards = MixedStandards,
                                         ppm = ppm,
                                         haldWidth = haldWidth,
                                         ignore_threshold = ignore_threshold,
                                         span = span)
      if (file.exists(gsub("_POS_", "_NEG_", mzMLFile))){
        ithAutomatedRtExtractorList_NEG <- PreliminaryRtExtractor(adductFile = NA,
                                           mzMLFile = gsub("_POS_", "_NEG_", mzMLFile),
                                           MixedStandards = MixedStandards,
                                           ppm = ppm,
                                           haldWidth = haldWidth,
                                           ignore_threshold = ignore_threshold,
                                           span = span)
        ithAutomatedRtExtractorList <- ithAutomatedRtExtractorList_POS
        for (ithList in c(1:length(ithAutomatedRtExtractorList))){
          if (!is.na(ithAutomatedRtExtractorList_POS[[ithList]][["Addcuts"]][1]) &
              !is.na(ithAutomatedRtExtractorList_NEG[[ithList]][["Addcuts"]][1])){
            ithAutomatedRtExtractorList[[ithList]][["Addcuts"]] <- c(ithAutomatedRtExtractorList_POS[[ithList]][["Addcuts"]], ithAutomatedRtExtractorList_NEG[[ithList]][["Addcuts"]])
          } else if(!is.na(ithAutomatedRtExtractorList_POS[[ithList]][["Addcuts"]][1]) &
                    is.na(ithAutomatedRtExtractorList_NEG[[ithList]][["Addcuts"]][1])){
            ithAutomatedRtExtractorList[[ithList]][["Addcuts"]] <- ithAutomatedRtExtractorList_POS[[ithList]][["Addcuts"]]
          } else if(is.na(ithAutomatedRtExtractorList_POS[[ithList]][["Addcuts"]][1]) &
                    !is.na(ithAutomatedRtExtractorList_NEG[[ithList]][["Addcuts"]][1])){
            ithAutomatedRtExtractorList[[ithList]][["Addcuts"]] <-  ithAutomatedRtExtractorList_NEG[[ithList]][["Addcuts"]]
          } else {
            print("Wrong!")
          }

        }
      } else {
        ithAutomatedRtExtractorList <- ithAutomatedRtExtractorList_POS
      }
      AutomatedRtExtractorList <- c(AutomatedRtExtractorList, ithAutomatedRtExtractorList)
    }
    save(AutomatedRtExtractorList, file = paste0(ResultPath, "/1_PreliminaryRtExtractor_", ithCS, ".Rdata"))
  }
}
