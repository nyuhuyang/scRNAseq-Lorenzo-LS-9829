########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.0.3
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","sctransform"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)


########################################################################
#
#  1 Seurat Alignment
#
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/20210408_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) <- colnames(df_samples) %>% tolower
rm = "LS1"
df_samples = df_samples[!(df_samples$sample %in% rm),]

print(df_samples)
(samples = df_samples$sample)
nrow(df_samples)

#======1.2 load  Seurat =========================
(load(file = "data/Lorenzo-LS6_20210408.Rda"))

meta.data = object@meta.data
for(i in 1:length(samples)){
    cells <- meta.data$orig.ident %in% samples[i]
    meta.data[cells,"conditions"] = df_samples$conditions[i]
}
meta.data$orig.ident %<>% factor(levels = samples)
meta.data$conditions %<>% factor(levels = c("WT", "KO"))
table(meta.data$conditions)

table(rownames(object@meta.data) == rownames(meta.data))
table(colnames(object) == rownames(meta.data))

object@meta.data = meta.data
#======1.6 Performing SCTransform and integration =========================
set.seed(100)
object_list <- SplitObject(object, split.by = "orig.ident")
remove(object);GC()

object_list %<>% lapply(SCTransform,method = "glmGamPoi")
object.features <- SelectIntegrationFeatures(object.list = object_list,
                                             nfeatures = 2000)
options(future.globals.maxSize= object.size(object_list)*1.5)
npcs = 30
anchors <- FindIntegrationAnchors(object.list = object_list, 
                                         anchor.features = object.features)
remove(object_list);GC()
# this command creates an 'integrated' data assay
object <- IntegrateData(anchorset = anchors,normalization.method = "SCT")
remove(anchors);GC()
save(object, file = "data/Lorenzo-LS6_20210408.Rda")

# Perform an integrated analysis
# Now we can run a single integrated analysis on all cells!
    
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(object) <- "integrated"

# Run the standard workflow for visualization and clustering
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(npcs = 30, verbose = FALSE)
object %<>% RunUMAP(reduction = "pca", dims = 1:30)
system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:30))

object %<>% FindNeighbors(reduction = "pca", dims = 1:30)
object %<>% FindClusters(resolution = 0.5)

#=======1.9 summary =======================================
lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun)
    fun(object, group.by="orig.ident",pt.size = 0.1,label = F,
        label.repel = T,alpha = 0.9,cols = c(Singler.colors,Singler.colors),
        no.legend = T,label.size = 4, repel = T, title = "Seurat Integration",
        do.print = T, do.return = F))

g1 <- UMAPPlot.1(object, group.by="orig.ident",pt.size = 0.1,label = F,
                 label.repel = T,alpha = 0.9,cols = c(Singler.colors,Singler.colors),
                 no.legend = T,label.size = 4, repel = T, title = "with data Integration",
                 do.print = F, do.return = T)

object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)

save(object, file = "data/Lorenzo-LS6_20210408.Rda")
format(object.size(object@assays$RNA),unit = "GB")
format(object.size(object@assays$integrated),unit = "GB")

object[['RNA']] <- NULL
object[['integrated']] <- NULL
   
format(object.size(object),unit = "GB")

save(object, file = "data/Lorenzo-LS6_20210408_SCT.Rda")

#= compare data integration==============
rm(object);GC()
load(file = "data/Lorenzo-LS6_20210408_SCT.Rda")

# Run the standard workflow for visualization and clustering
object %<>% RunPCA(npcs = 30, verbose = FALSE)
object %<>% RunUMAP(reduction = "pca", dims = 1:30)

g2 <- UMAPPlot.1(object, group.by="orig.ident",pt.size = 0.1,label = F,
                 label.repel = T,alpha = 0.9,cols = c(Singler.colors,Singler.colors),
                 no.legend = T,label.size = 4, repel = T, title = "without Integration",
                 do.print = F, do.return = T)

jpeg(paste0(path,"UMAP_data_integration.jpeg"), units="in", width=10, height=7,res=600)
patchwork::wrap_plots(list(g2,g1), nrow = 1, ncol = 2)
dev.off()
