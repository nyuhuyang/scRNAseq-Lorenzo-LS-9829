library(Seurat)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# ======== expression
load(file = "data/Lorenzo-LS6_20210408_SCT.Rda")
object$label.fine %<>% gsub("Macrophages activated","Macrophages",.)

Idents(object) = "label.fine"
EC <- subset(object, idents = "Endothelial cells")
Idents(EC) = "conditions"
PrepareGSEA(EC, k=10, file.name = "EC_KO_WT")

MF <- subset(object, idents = "Macrophages")
Idents(MF) = "conditions"
PrepareGSEA(MF, k=10, file.name = "MF_KO_WT")

#========== prerank
human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

csv_list <- list.files(pattern = "FC0.01",path = "output/20210511",full.names = T)
deg_list <- pbapply::pblapply(csv_list, function(x){
    tmp = read.csv(x,row.names = 1)
    tmp = tmp[tmp$cluster == "KO", ]
    tmp$avg_log2FC = tmp$avg_log2FC * -log10(tmp$p_val)
    tmp = tmp[order(tmp$avg_log2FC,decreasing = T), c("gene","avg_log2FC")]
    tmp$gene %<>% toupper()
    
    return(tmp)
})
names(deg_list) = c("Endothelial_cells","Macrophages")
lapply(names(deg_list), function(x) write.table(deg_list[[x]],row.names = F,col.names = F,
                                                sep = "\t",quote = F,
                                              file = paste0(path,x,"_KO_WT.rnk")))

