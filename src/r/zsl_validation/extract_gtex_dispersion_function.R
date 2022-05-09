library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/input/"
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/output/"

dds <- readRDS(paste(OUT_DIR,"dds.Rds",sep=""))

modeled_dds <- DESeq2::estimateSizeFactors(dds)
modeled_dds <- DESeq2::estimateDispersions(modeled_dds, fitType="local")
dispersion_function <- DESeq2::dispersionFunction(modeled_dds)

vst.counts <- DESeq2::vst(dds,
                          blind = FALSE,
                          fitType = "local")

DESeq2::dispersionFunction(dds) <- dispersion_function

vst.counts.doublecheck <- DESeq2::vst(dds,
                                      blind = FALSE)

stopifnot(assertthat::are_equal(vst.counts,vst.counts.doublecheck))

print("Precalculated dispersion function recapitulates local vst().")

saveRDS(dispersion_function,
        paste(OUT_DIR,"dispersion_function.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))