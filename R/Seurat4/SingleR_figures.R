# conda activate r4.0.3
library(Seurat)
library(magrittr)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggpubr)
library(S4Vectors)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================

##############################
# create singleR data frame
###############################
pred = readRDS("output/Lorenzo-LS6_20210411_singleR_pred.rds")
object = readRDS(file = "data/Lorenzo-LS6_20210411.rds")

singlerDF = data.frame("label.fine" = pred$pruned.labels,
                       row.names = rownames(pred))
singlerDF$label.fine[is.na(singlerDF$label.fine)]= "unknown"
singlerDF$label.main = gsub(" activated$","", singlerDF$label.fine)

##############################
# process color scheme
##############################
table(colnames(object) == rownames(singlerDF))
object@meta.data %<>% cbind(singlerDF[,c("label.main","label.fine")])
saveRDS(object, file = "data/Lorenzo-LS6_20210411.rds")
