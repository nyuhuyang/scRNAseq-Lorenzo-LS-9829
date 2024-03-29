#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.0.3 linux
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(magrittr)

source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

# ====== load single cell =============
object = readRDS(file = "data/Lorenzo-LS6_20210411.rds")

sce <- SingleCellExperiment(list(logcounts=object[["SCT"]]@data),
                                colData=DataFrame(object@meta.data))
rm(object);GC()

# ====== load reference =============
ref <- celldex::MouseRNAseqData()
table(ref$label.fine)
system.time(pred <- SingleR(test = sce, ref = ref, assay.type.test=1,
                                 labels = ref$label.fine))
# elapsed 4872.846 sec
saveRDS(object = pred, file = "output/Lorenzo-LS6_20210411_singleR_pred.rds")
