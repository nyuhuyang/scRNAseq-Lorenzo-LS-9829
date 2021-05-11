########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
library(Seurat)
library(dplyr)
library(kableExtra)
library(magrittr)
library(ggplot2)
library(cowplot)
library(fgsea)
library(tibble)
library(ggsci)
library(progress)

source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- "Yang/Figure 2/Figure Sources/"
if(!dir.exists(path)) dir.create(path, recursive = T)

csv_list <- list.files(pattern = "FC0.1",path = "output/20210428",full.names = T)
deg_list <- pbapply::pblapply(csv_list, function(x){
    tmp = read.csv(x,row.names = 1)
    tmp = tmp[tmp$cluster == "KO", ]
    tmp$cluster %<>% paste0("/WT")
    tmp = tmp[order(tmp$avg_log2FC,decreasing = T), ]
    #tmp = tmp[tmp$avg_log2FC > 0, ]
    tmp
    
})


head(deg_list[[1]], 20)
# read pathway

hallmark <- fgsea::gmtPathways("../seurat_resources/msigdb/h.all.v7.3.symbols.gmt")
names(hallmark) = gsub("HALLMARK_","",names(hallmark))
names(hallmark) = gsub("\\_"," ",names(hallmark))


hallmark <- fgsea::gmtPathways("../seurat_resources/msigdb/h.all.v7.3.symbols.gmt")
Biocarta <- gmtPathways("../seurat_resources/msigdb/c2.cp.biocarta.v7.3.symbols.gmt")
kegg <- gmtPathways("../seurat_resources/msigdb/c2.cp.kegg.v7.3.symbols.gmt")
pid <- gmtPathways("../seurat_resources/msigdb/c2.cp.pid.v7.3.symbols.gmt")
reactome <- gmtPathways("../seurat_resources/msigdb/c2.cp.reactome.v7.3.symbols.gmt")
tft <- gmtPathways("../seurat_resources/msigdb/c3.tft.v7.3.symbols.gmt")
Development <- gmtPathways("../seurat_resources/msigdb/Development.gmt")
Physiology <- gmtPathways("../seurat_resources/msigdb/Physiology.gmt")
Disease <- gmtPathways("../seurat_resources/msigdb/Disease.gmt")
Cancer <- gmtPathways("../seurat_resources/msigdb/Cancer.gmt")

msigdb_list <- list("hallmark" = hallmark,
                    "Biocarta" = Biocarta,
                    "kegg" = kegg,
                    "pid" = pid,
                    "reactome" = reactome,
                    "transcription_factor_targets" = tft,
                    "Development" = Development,
                    "Physiology" = Physiology,
                    "Disease" = Disease,
                    "Cancer" = Cancer)
msigdb_list %<>% pbapply::pblapply(function(x) {
    x %<>% lapply(function(x1) Hmisc::capitalize(tolower((x1))))
    x
})

hallmark %<>% lapply(function(x) Hmisc::capitalize(tolower((x))))
#=========
res =  bind_rows(deg_list)
res$cluster = res$cell.type
colnames(res)[2] = "avg_logFC"
(clusters = unique(as.character(res$cluster)))

Fgsea_res <- FgseaDotPlot(stats=res, order.yaxis.by = c(1,"avg_logFC"), pathways=hallmark,
                 title = "enriched hallmark pathways in KO",width = 6,do.return = T)
Fgsea_res = Fgsea_res[order(Fgsea_res$pathway),]
write.csv(Fgsea_res, file = paste0(path,"hallmark_gsea.csv"))
