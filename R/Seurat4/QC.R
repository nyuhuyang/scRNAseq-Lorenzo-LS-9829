########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","ggplot2","scater"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
        }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("data")) dir.create("data")
if(!dir.exists("doc")) dir.create("doc")
########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======1.1 Load the data files and Set up Seurat object =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/20210408_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) <- colnames(df_samples) %>% tolower
rm = "LS1"
df_samples = df_samples[!(df_samples$sample %in% rm),]

print(df_samples)
(samples = df_samples$sample)
nrow(df_samples)

# check missing data
current <- list.files("data/counts")
(current <- current[!grepl(".Rda|RData",current)])
(missing_data <- df_samples$sample.id[!(df_samples$sample.id %in% current)])

#======1.1.2 record data quality before removing low quanlity cells =========================
# if args 2 is passed
message("read metrics_summary")
QC_list <- lapply(df_samples$sample.id, function(x){
        tmp = read.csv(file = paste0("data/counts/",x,
                               "/outs/metrics_summary.csv"))
        t(tmp)
})

QC = bind_cols(QC_list) %>% as.data.frame()
colnames(QC) = df_samples$sample.id
rownames(QC) = rownames(QC_list[[1]])
QC["Sequencing.Saturation",] %>% gsub("%","",.) %>% as.numeric %>% mean
QC["Estimated.Number.of.Cells",] %>% gsub(",","",.) %>% as.numeric %>% sum

write.csv(QC,paste0(path,"metrics_summary.csv"))

message("Loading the datasets")
## Load the dataset
Seurat_raw <- list()
Seurat_list <- list()
for(i in seq_along(df_samples$sample)){
        Seurat_raw[[i]] <- Read10X(data.dir = paste0("data/counts/",
                                                     df_samples$sample.id[i],"/",
                                                     df_samples$read.path[i]))
        colnames(Seurat_raw[[i]]) = paste0(df_samples$sample[i],"-",colnames(Seurat_raw[[i]]))
        colnames(Seurat_raw[[i]]) %<>% gsub("-[0-9+]","",.)
        Seurat_list[[i]] <- CreateSeuratObject(Seurat_raw[[i]],
                                               min.cells = 0,
                                               names.delim = "-",
                                               min.features = 0)
        Progress(i, length(df_samples$sample))
}
remove(Seurat_raw);GC()
#========1.1.3 g1 QC plots before filteration=================================
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), Seurat_list)
remove(Seurat_list);GC()


# read and select mitochondial genes
if(unique(df_samples$organism) == "hg19") mito = "^MT-"
if(unique(df_samples$organism) == "mm10") mito = "^mt-" # not Mt-
if(unique(df_samples$organism) == "danRer10") mito = "^mt-"
message("mito.genes:")

(mito.features <- grep(pattern = mito, x = rownames(object), value = TRUE))
object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = mito)
Idents(object) %<>% factor(levels = df_samples$sample)
g1 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
        VlnPlot(object = object, features = features, ncol = 1, pt.size = 0.01)+
                theme(axis.text.x = element_text(size=15,angle = 0,hjust = 0.5),legend.position="none")
})
save(g1,file= paste0(path,"g1","_",length(df_samples$sample),"_",gsub("-","",Sys.Date()),".Rda"))

#============1.2 scatter ======================
meta_data = object@meta.data
meta_data$discard = FALSE
for(i in 1:length(df_samples$sample)){
        cell = rownames(meta_data)[meta_data$orig.ident %in% df_samples$`sample`[i]]
        high.mito <- isOutlier(meta_data[cell,"percent.mt"], nmads=3, type="higher")
        low.lib <- isOutlier(log10(meta_data[cell,"nCount_RNA"]), type="lower", nmad=3)
        low.genes <- isOutlier(log10(meta_data[cell,"nFeature_RNA"]), type="lower", nmad=3)
        discard <- high.mito | low.lib | low.genes
        print(data.frame(HighMito= sum(high.mito),LowLib=sum(low.lib),
                         LowNgenes=sum(low.genes),Discard=sum(discard)))
        meta_data[cell,"discard"] = discard
}
object@meta.data = meta_data
table(object$orig.ident, object$discard)

object %<>% subset(subset = discard == FALSE & 
                           nFeature_RNA > 700 &
                           nCount_RNA > 1000 & percent.mt < 10)
# FilterCellsgenerate Vlnplot before and after filteration
Idents(object) = "orig.ident"
Idents(object) %<>% factor(levels = df_samples$sample)

g2 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
        VlnPlot(object = object, features = features, ncol = 1, pt.size = 0.01)+
                theme(axis.text.x = element_text(size=15,angle = 0,hjust = 0.5),legend.position="none")
})
save(g2,file= paste0(path,"g2","_",length(df_samples$sample),"_",gsub("-","",Sys.Date()),".Rda"))
#load(paste0(save.path,"g2_58_20201009.Rda"))

jpeg(paste0(path,"S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nFeature_RNA before filteration")+
                        scale_y_log10(limits = c(100,10000))+
                        theme(plot.title = element_text(hjust = 0.5)),
                g2[[1]]+ggtitle("nFeature_RNA after filteration")+
                        scale_y_log10(limits = c(100,10000))+
                        theme(plot.title = element_text(hjust = 0.5))))
dev.off()
jpeg(paste0(path,"S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nCount_RNA before filteration")+
                        scale_y_log10(limits = c(500,100000))+
                        theme(plot.title = element_text(hjust = 0.5)),
                g2[[2]]+ggtitle("nCount_RNA after filteration")+
                        scale_y_log10(limits = c(500,100000))+
                        theme(plot.title = element_text(hjust = 0.5))))
dev.off()
jpeg(paste0(path,"S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % before filteration")+
                        ylim(c(0,50))+
                        theme(plot.title = element_text(hjust = 0.5)),
                g2[[3]]+ggtitle("mito % after filteration")+
                        ylim(c(0,50))+
                        theme(plot.title = element_text(hjust = 0.5))))
dev.off()

#====
format(object.size(object),unit = "GB")
save(object, file = "data/Lorenzo-LS6_20210408.Rda")