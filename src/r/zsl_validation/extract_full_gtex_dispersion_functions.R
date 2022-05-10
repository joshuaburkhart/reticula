library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

MAXIT <- 2
IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/input/"
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/output/"

varStabFunc_w_MAXIT <- function (object, blind = TRUE, fitType = "parametric") 
{
  if (is.null(colnames(object))) {
    colnames(object) <- seq_len(ncol(object))
  }
  if (is.matrix(object)) {
    matrixIn <- TRUE
    object <- DESeqDataSetFromMatrix(object, DataFrame(row.names = colnames(object)), 
                                     ~1)
  }
  else {
    matrixIn <- FALSE
  }
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    object <- estimateSizeFactors(object)
  }
  if (blind) {
    design(object) <- ~1
  }
  if (blind | is.null(attr(dispersionFunction(object), "fitType"))) {
    object <- estimateDispersionsGeneEst(object, quiet = TRUE,
                                         maxit=MAXIT)
    object <- estimateDispersionsFit(object, quiet = TRUE, 
                                     fitType,
                                     maxit=MAXIT)
  }
  vsd <- getVarianceStabilizedData(object)
  if (matrixIn) {
    return(vsd)
  }
  se <- SummarizedExperiment(assays = vsd, colData = colData(object), 
                             rowRanges = rowRanges(object), metadata = metadata(object))
  DESeqTransform(se)
}

dds <- readRDS(paste(OUT_DIR,"dds.Rds",sep=""))

modeled_dds <- DESeq2::estimateSizeFactors(dds)

ggp_fit_modeled_dds <- DESeq2::estimateDispersions(modeled_dds, fitType="glmGamPoi",maxit=MAXIT)
local_fit_modeled_dds <- DESeq2::estimateDispersions(modeled_dds, fitType="local",maxit=MAXIT)

ggp_fit_dispersion_function <- DESeq2::dispersionFunction(ggp_fit_modeled_dds)
local_fit_dispersion_function <- DESeq2::dispersionFunction(local_fit_modeled_dds)

saveRDS(ggp_fit_dispersion_function,
        paste(OUT_DIR,"ggp_fit_dispersion_function.Rds",sep=""))
saveRDS(local_fit_dispersion_function,
        paste(OUT_DIR,"local_fit_dispersion_function.Rds",sep=""))

# # check dispersion function assignment matches parameterized VST output
# # how does parameterized VST pass maxit to 
# 
# vst.counts <- DESeq2::varianceStabilizingTransformation(dds,
#                                                         blind = FALSE,
#                                                         fitType = "local")
# 
# DESeq2::dispersionFunction(dds) <- local_fit_dispersion_function
# 
# vst.counts.doublecheck <- DESeq2::varianceStabilizingTransformation(dds,
#                                                                     blind = FALSE)
# 
# stopifnot(assertthat::are_equal(vst.counts,vst.counts.doublecheck))
# 
# print("Precalculated dispersion function recapitulates local vst().")

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))