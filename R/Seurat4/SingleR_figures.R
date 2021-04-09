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
pred = readRDS("output/20210408_Lorenzo-LS6_singleR_pred.rds")
load(file = "data/Lorenzo-LS6_20210408_SCT.Rda")

singlerDF = data.frame("label.fine" = pred$pruned.labels,
                       row.names = rownames(pred))
singlerDF$label.fine[is.na(singlerDF$label.fine)]= "unknown"
singlerDF$label.main = gsub(" activated$","", singlerDF$label.fine)

##############################
# process color scheme
##############################
table(colnames(object) == rownames(singlerDF))
object@meta.data %<>% cbind(singlerDF[,c("label.main","label.fine")])

object %<>% AddMetaColor(label = "label.fine",colors = Singler.colors)
lapply(c("label.fine","label.main"), function(group.by)
    TSNEPlot.1(object = object, label = T, label.repel = T,group.by = group.by,
        no.legend = T,
        pt.size = 0.1,label.size = ifelse(group.by=="label.fine",3,5),
        do.print = T,do.return = F,
        title = paste(group.by, "labeling by Mouse RNA-seq")))

lapply(c("label.fine","label.main"), function(group.by)
    UMAPPlot.1(object = object, label = T, label.repel = T,group.by = group.by,
               no.legend = T,
               pt.size = 0.1,label.size = ifelse(group.by=="label.fine",3,5),
               do.print = T,do.return = F,
               title = paste(group.by, "labeling by Mouse RNA-seq")))
save(object, file = "data/Lorenzo-LS6_20210408_SCT.Rda")


# by barplot
cell_Freq <- table(object$label.fine) %>% as.data.frame
cell_Freq$Percent <- prop.table(cell_Freq$Freq) %>% round(digits = 2) %>% scales::percent()
cols = ExtractMetaColor(object)
cell_Freq$cols = cols[cell_Freq$Var1]
cell_Freq = cell_Freq[order(cell_Freq$Var1),]

cell_Freq = cell_Freq[order(cell_Freq$Freq,decreasing = T),]
cell_Freq$Var1 %<>% factor(levels = as.character(cell_Freq$Var1))
colnames(cell_Freq)[1:2] = c("Cell_Type", "Cell_Number")

jpeg(paste0(path,"cell_type_numbers.jpeg"), units="in", width=6, height=6,res=600)
ggbarplot(cell_Freq, "Cell_Type", "Cell_Number",
          fill = "Cell_Type", color = "black",xlab = "",
          palette = cell_Freq$col,x.text.angle = 90,
          ylab = "Cell Number",
          label = cell_Freq$Percent,
          lab.size = 3,
          sort.val = "desc",
          width = 1, size = 0.5,
          title = "Numbers of cell types in total 6 samples")+NoLegend()+
    theme(plot.title = element_text(hjust = 0.5,size=15),
          axis.text.x = element_text(vjust = 0.5))
dev.off()