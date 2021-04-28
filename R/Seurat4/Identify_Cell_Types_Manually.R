library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# ======== 2.1 =========== test with known markers==================
(load(file = "data/MCL_AIM_74_20210311_SCT.Rda"))
DefaultAssay(object) <- "SCT"
Idents()
# global
features <- c("CD19","CCND1","PCNA",
              "CD3D","CD4","CD8A",
              "MS4A7","CD14","FCGR1A",
              "GNLY","KLRC1","NCAM1")
FeaturePlot.1(object,features = features, pt.size = 0.005,
              cols = c("lightgrey", "red"),
              alpha = 1,reduction = "tsne",
              unique.name = "cell.types",label = F,
              threshold = 1, text.size = 20, border = T,
              file.name="TSNE_AIM_74_markers.jpeg",do.print = T, do.return = F,ncol = 3,
              units = "in",width=9, height=12, no.legend = T)
QC <- c("percent.mt","nCount_SCT","nFeature_SCT")
FeaturePlot.1(object,features = QC, pt.size = 0.005,
              cols = c("lightgrey", "red"),
              alpha = 1,reduction = "tsne",
              unique.name = F,label = F,
              threshold = 1, text.size = 20, border = T,do.print = T, do.return = F,ncol = 3,
              units = "in",width=9, height=4, no.legend = T)
cc <- c("CCND1","CDK4","RB1","E2F1","MCM7","CCNB2")
FeaturePlot.1(object,features = cc, pt.size = 0.005,
              cols = c("lightgrey", "red"),
              alpha = 1,reduction = "tsne",
              unique.name = "cell.types",label = F,
              threshold = 1, text.size = 20, border = T,do.print = T, do.return = F,ncol = 2,
              units = "in",width=9, height=9, no.legend = T)

write.csv(as.data.frame.table(table(object$label.fine)),file = paste0(path,"label.fine.csv"))
object$label = object$label.fine
object$label %<>% gsub("T cells, CD4+.*","T cells, CD4+",.)
object$label %<>% gsub("T cells, CD8+.*","T cells, CD8+",.)
object$label %<>% gsub("Monocytes.*","Monocytes",.)
object$label %<>% gsub("B cells,.*","B cells",.)
write.csv(as.data.frame.table(table(object$label)),file = paste0(path,"label.csv"))

labels = c("T cells, CD4+","T cells, CD8+","NK cells","B cells","MCL","Monocytes")
Idents(object) = "label"

exp_list <- vector("list", length = length(labels))
names(exp_list) = labels
for(i in seq_along(labels)){
    sub_object <- subset(object, idents = labels[i])
    Idents(sub_object) = "orig.ident"
    cell.number.df <- as.data.frame.table(table(sub_object$orig.ident))
    cell.number = cell.number.df$Freq
    names(cell.number) = cell.number.df$Var1
    exp = AverageExpression(sub_object, assays = "SCT")
    exp_list[[labels[i]]] = rbind("cell.number" = cell.number, exp$SCT)
    svMisc::progress(i/length(exp_list)*100)

}

openxlsx::write.xlsx(exp_list,
                     file = paste0(path,"Expression_Cell.types.xlsx"),
                     colNames = TRUE, rowNames = TRUE,
                     borders = "surrounding",colWidths = c(NA, "auto", "auto"))

for(i in seq_along(exp_list)){
    write.table(exp_list[[i]],file = paste0(path,"Expression_",names(exp_list)[i],".txt"),
                sep = "\t",row.names = TRUE)
    svMisc::progress(i/length(exp_list)*100)
}

# Extend Data
Idents(object) = "cell.types"
object %<>% sortIdent()

cell_Freq <- table(Idents(object)) %>% as.data.frame
cell_Freq$Percent <- round(prop.table(cell_Freq$Freq),digits = 3) %>% 
    scales::percent()

cell_Freq = cell_Freq[order(cell_Freq$Var1),]
df_samples <- readxl::read_excel("doc/singler.colors.xlsx")
cell_Freq$col = plyr::mapvalues(cell_Freq$Var1,
                                from = na.omit(df_samples$Cell.types),
                                to = na.omit(df_samples$singler.color2))
cell_Freq$col %<>% as.character
cell_Freq = cell_Freq[order(cell_Freq$Freq,decreasing = T),]
cell_Freq$Var1 %<>% factor(levels = as.character(cell_Freq$Var1))
colnames(cell_Freq)[1:2] = c("Cell_Type", "Cell_Number")
cell_Freq$Cell_Type %<>% gsub("_"," ",.)


jpeg(paste0(path,"cell_type_numbers.jpeg"), units="in", width=6, height=6,res=600)
ggbarplot(cell_Freq, "Cell_Type", "Cell_Number",
          fill = "Cell_Type", color = "black",xlab = "",
          palette = cell_Freq$col,x.text.angle = 45,
          ylab = "Cell Number",
          label = cell_Freq$Percent,
          sort.val = "desc",
          width = 1, size = 0.5,
          title = "Numbers of major cell types in total 74 samples")+NoLegend()+
    theme(plot.title = element_text(hjust = 0.5,size=15))+
    scale_y_continuous(expand = c(0, 0), limits = c(0,max(cell_Freq$Cell_Number)+4000))
dev.off()
