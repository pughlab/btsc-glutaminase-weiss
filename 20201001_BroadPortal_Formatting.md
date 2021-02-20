
---
# Reformat data for Broad Portal Upload

L.Richards  
October 2020  

---

Format data used in Ians paper for Broad portal upload.


---
### 1.0 Data File
---

Transferred to H4H.
/cluster/projects/pughlab/projects/Weiss_Astro_scRNAseq/BroadPortal


```R
library("data.table")
library("dplyr")
library("taRifx")
library("Matrix")
#suppressMessages(library("Seurat", 
#                         lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/" ))
library(Seurat)


#SeuratPath <- "/mnt/work1/users/pughlab/projects/Weiss_Astro_scRNAseq/BTSC_Weiss_merged.RData"
SeuratPath <- "/cluster/projects/pughlab/projects/Weiss_Astro_scRNAseq/data/BTSC_Weiss_merged.RData"
sample <- "Restall_BTSC_Glutaminase"



################################################
print("")
print("-------")
print(sample)
start.time <- Sys.time()
print(start.time)
print("-------")

################################################
print("")
print("")
print("Loading Seurat Object")
print(Sys.time())
print(SeuratPath)

loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}
      
dat <- loadRData(SeuratPath)
print(dim(dat@data))


```


```R
################################################
### Read in Ian's metadata
ian.meta <- read.csv("/cluster/projects/pughlab/projects/Weiss_Astro_scRNAseq/data/CanRes2020_scRNAseq_4lines.csv")
dim(ian.meta)
```


```R
################################################
print("")
print("")
print("Extract Normalized Gene Expression File (txt.gz)")
print(Sys.time())

norm.exp <- data.frame(as.matrix(dat@data))
table(colnames(norm.exp) %in% ian.meta$Cell)
norm.exp <- norm.exp[ ,as.character(ian.meta$Cell)] ##subset by cells Ian included in csv
print(norm.exp[1:5, 1:5])

norm.exp <- tibble::rownames_to_column(norm.exp, "GENE")
print(norm.exp[1:5, 1:5])

file1 <- paste0(sample, "_NormalizedExpression.txt")
print(file1)
      
write.table(norm.exp, 
          file=file1,
          row.names = FALSE,
          quote = FALSE,
            sep = "\t"
          )      
```


```R
################################################
print("")
print("")
print("Extract Meta Data (.txt)")
print(Sys.time())

meta <- data.frame(dat@meta.data)
meta <- remove.factors(meta)
meta <- meta[as.character(ian.meta$Cell), ]
meta <- cbind(meta, ian.meta)

include <- c("nGene",
             "nUMI",
             "percent.mito",
             "BTSC_line",
             "EAAT12_sig",
             "Astro_sig"
            )
meta <- meta[,include]

colnames(meta)[4] <- "Sample"
colnames(meta)[5] <- "EAAT12_GeneSignature_Score"
colnames(meta)[6] <- "Astrocyte_GeneSignature_Score"

meta$nGene <- as.numeric(meta$nGene)
meta$nUMI <- as.numeric(meta$nUMI)
meta$EAAT12_GeneSignature_Score <- as.numeric(meta$EAAT12_GeneSignature_Score)
meta$Astrocyte_GeneSignature_Score <- as.numeric(meta$Astrocyte_GeneSignature_Score)
meta$Sample <- as.character(meta$Sample)

classes <- as.vector(unlist(lapply(meta, class)))
TYPE <- c("TYPE", ifelse(classes == "numeric", "numeric", "group"))
meta <- data.frame(setDT(meta, keep.rownames = "NAME"))
meta <- rbind(TYPE, meta)
      
file2 <- paste0(sample, "_MetaData.txt")
print(file2)

write.table(meta, 
          file=file2,
          row.names = FALSE,
          quote = FALSE,
          sep = "\t",
          col.names = TRUE
          ) 
```


```R
################################################
### Recluster data

dat <- UpdateSeuratObject(dat)
sub <- subset(dat, idents = c("BT48", "BT67", "BT89", "BT94"))
all.genes <- rownames(sub)
sub <- ScaleData(sub, features = all.genes)
sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = 2000)
sub <- RunPCA(sub, features = VariableFeatures(object = sub))
sub <- FindNeighbors(sub, dims = 1:10)
sub <- FindClusters(sub, resolution = 0.5)
sub <- RunUMAP(sub, dims = 1:10)
```


```R
################################################
print("")
print("")
print("Extract UMAP Coordinates (.txt)")
print(Sys.time())

cluster.coords <- as.data.frame(sub@reductions$umap@cell.embeddings)[, c(1:2)]
colnames(cluster.coords) <- c("X", "Y")
cluster.coords <- remove.factors(cluster.coords)

classes <- as.vector(unlist(lapply(cluster.coords, class)))
TYPE <- c("TYPE", ifelse(classes == "numeric", "numeric"))
cluster.coords <- data.frame(setDT(cluster.coords, keep.rownames = "NAME"))
cluster.coords  <- rbind(TYPE, cluster.coords)

##subset out some meta data to include as labels
cluster.coords <- cbind(cluster.coords, meta[ ,c("Sample")] )
colnames(cluster.coords)[4] <- "Sample"                                      

file3 <- paste0(sample, "_UMAP_ClusterCoordinates.txt")
print(file3)
write.table(cluster.coords, 
          file=file3,
          row.names = FALSE,
          quote = FALSE,
           sep = "\t",
          col.names = TRUE
          ) 
```
