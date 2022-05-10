library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/input/"
GTEX_OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/output/"
TCGA_OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/tcga/output/"
DISP_FUNC_SAVE_PATH <- paste(GTEX_OUT_DIR,"sampled_local_dispersion_function.Rds",sep="")

vst2 <- function (object, blind = TRUE, nsub = 1000, fitType = "parametric", dispFuncSavePath = paste("./",fitType,".Rds",sep="")) 
{
  if (nrow(object) < nsub) {
    stop("less than 'nsub' rows,\n  it is recommended to use varianceStabilizingTransformation directly")
  }
  if (is.null(colnames(object))) {
    colnames(object) <- seq_len(ncol(object))
  }
  if (is.matrix(object)) {
    matrixIn <- TRUE
    object <- DESeqDataSetFromMatrix(object, DataFrame(row.names = colnames(object)), 
                                     ~1)
  }
  else {
    if (blind) {
      design(object) <- ~1
    }
    matrixIn <- FALSE
  }
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    object <- estimateSizeFactors(object)
  }
  baseMean <- rowMeans(counts(object, normalized = TRUE))
  if (sum(baseMean > 5) < nsub) {
    stop("less than 'nsub' rows with mean normalized count > 5, \n  it is recommended to use varianceStabilizingTransformation directly")
  }
  object.sub <- object[baseMean > 5, ]
  baseMean <- baseMean[baseMean > 5]
  o <- order(baseMean)
  idx <- o[round(seq(from = 1, to = length(o), length = nsub))]
  object.sub <- object.sub[idx, ]
  object.sub <- estimateDispersionsGeneEst(object.sub, quiet = TRUE)
  object.sub <- estimateDispersionsFit(object.sub, fitType = fitType, 
                                       quiet = TRUE)
  suppressMessages({
    dispersionFunction(object) <- dispersionFunction(object.sub)
    saveRDS(dispersionFunction(object),file=dispFuncSavePath)
  })
  vsd <- varianceStabilizingTransformation(object, blind = FALSE)
  if (matrixIn) {
    return(assay(vsd))
  }
  else {
    return(vsd)
  }
}

gtex_dds <- readRDS(paste(GTEX_OUT_DIR,"dds.Rds",sep=""))
gtex.vst.counts <- vst2(gtex_dds,
                       blind = FALSE,
                       fitType = "local",
                       dispFuncSavePath = DISP_FUNC_SAVE_PATH)

tcga_dds <- readRDS(paste(TCGA_OUT_DIR,"dds.Rds",sep=""))
DESeq2::dispersionFunction(tcga_dds) <- readRDS(DISP_FUNC_SAVE_PATH)

tcga.vst.counts <- varianceStabilizingTransformation(tcga_dds, blind = FALSE)

saveRDS(gtex.vst.counts, file=paste(GTEX_OUT_DIR,"gtex_vst_counts.Rds",sep=""))
saveRDS(tcga.vst.counts, file=paste(TCGA_OUT_DIR,"tcga_vst_counts.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))