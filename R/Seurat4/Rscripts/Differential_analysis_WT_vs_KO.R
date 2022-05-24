########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot",
                   "magrittr","data.table","future","ggplot2","tidyr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# Need 64GB
# load files
object  =readRDS(file = "data/Lorenzo-LS6_20220411.rds")
# Need 32GB
#DefaultAssay(object) = "SCT"
#Idents(object) = "Doublets"
#object <- subset(object, idents = "Singlet")
object$label.fine %<>% gsub("Macrophages activated","Macrophages",.)
object$label.fine %<>% gsub("Fibroblasts.*","Fibroblasts",.)

Idents(object) = "label.fine"
cell.types <- c("B cells","Dendritic cells","Endothelial cells",
                "Fibroblasts","Granulocytes","Macrophages","Monocytes","NK cells",
                "T cells")
cell.type = cell.types[args]
object <- subset(object, idents = cell.type)
Idents(object) = "condition"
system.time(markers <- FindAllMarkers(object, 
                                       logfc.threshold = 0.01, 
                                       return.thresh = 1, 
                                       only.pos = F,latent.vars = "nFeature_SCT",
                                        test.use = "MAST"))
markers$cell.type = cell.type
if(args < 10) args = paste0("0", args)
write.csv(markers,paste0(path,args,"_FC0.01_",cell.type,".csv"))