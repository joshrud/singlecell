# for running Seurat integration in parallel
# Instructions for running: 
# 1. change dir.main so that it it is a directory where a "saves/" folder is, with the un-integrated object we're integrating
# 2. Run on command line via "Rscript parallelIntegration_final.R unintegrated_object.R"

library(Seurat)
library(future)

# check if multicore is supported ----
if (!future::supportsMulticore()) {
  print("doesn't support multicore...")
  quit(save = "no", status = 1)
} else {
  print("supports multicore! running... :)")
}

# initialization ----
# options(future.globals.maxSize = 6815744000) #for doing scale.data with 12 threads
options(future.globals.maxSize = 11*1024^3) #1024 ^3 is GiB
args <- commandArgs(trailingOnly = T)
if (args[1] == "") {
  stop("no argument specified...")
}

# dirs ----
dir.main <-  gsub("[^\\/]+$", "", x=args[1]) # grabs the path from the argument 
# dir.saves <- paste0(dir.main, "saves/")

savefile.name = gsub(".*\\/", "", x=args[1])
tempfile.name <- gsub("\\.rds", "_TEMP.rds", savefile.name)
rdsIn <- paste0(dir.main,  savefile.name)
rdsOut <- paste0(dir.main, gsub("\\.rds", "_OUT.rds", savefile.name))

if (!file.exists(paste0(dir.main, savefile.name))) {
  stop("savefile doesn't exist, please choose another.")
}
if (!file.exists(paste0(dir.main, tempfile.name))) {
  message("haven't created integrated object, loading unintegrated object...")
  all.monomac <- readRDS(rdsIn)
  DefaultAssay(all.monomac) <- "RNA"
  
  if (!is.na(args[2])) { # added to allow splitting by other cols
    splitcol = args[2]
    message(paste0("Split Column given, using: ", splitcol)) 
  } else {
    splitcol = "orig.ident"
    message("Split Column not given, using DEFAULT: orig.ident")
  }
  
  all.monomac.list <- SplitObject(all.monomac, split.by = splitcol)
  message("normalizing...")
  for (i in 1:length(all.monomac.list)) {
    all.monomac.list[[i]] <- NormalizeData(all.monomac.list[[i]], verbose = FALSE)
    all.monomac.list[[i]] <- FindVariableFeatures(all.monomac.list[[i]], selection.method = "vst", 
                                               nfeatures = 3000, verbose = FALSE)
  }
  allcells = sum( sapply(all.monomac.list, ncol))
  if (allcells <= 1e5) {
    message("less than 100,000 cells, performing single parallel integration")
    message("planning...")
    plan("multicore", workers = 12, future.seed=T)
    message("finding integration anchors....")
    all.monomac.anchors <- FindIntegrationAnchors(object.list = all.monomac.list, 
                                                  dims = 1:30,
                                                  anchor.features=3000)
    message("integrating...")
    all.monomac.int <- IntegrateData(anchorset = all.monomac.anchors, dims = 1:30)
  } else {
    
    message("MORE than 100,000 cells, performing iterative parallel integration")

    for (i in seq_along(all.monomac.list)) {
      message(paste0("working on obj: ", i, " out of ", length(all.monomac.list)))
      message("======================================")
      
      if (i == 1) {
        all.monomac.int <- all.monomac.list[[1]]
        next
      }
      
      message("planning...")
      plan("multicore", workers = 12, future.seed=T)
      message("finding integration anchors....")
      all.monomac.anchors <- FindIntegrationAnchors(object.list = list(all.monomac.list[[i]], 
                                                                       all.monomac.int), 
                                                    dims = 1:30,
                                                    anchor.features=3000)
      message("integrating...")
      all.monomac.int <- IntegrateData(anchorset = all.monomac.anchors, dims = 1:30)
    }
    
  }
  DefaultAssay(all.monomac.int) <- "integrated"
  message("saving temp...")
  saveRDS(all.monomac.int, paste0(dir.main, tempfile.name))
  print(head(all.monomac.int@meta.data))
  print(paste0("ncount: ", any(is.na(all.monomac.int@meta.data$nCount_RNA))))
  print(paste0("mito: ", any(is.na(all.monomac.int@meta.data$percent.mito))))
  print(paste0("ribo: ", any(is.na(all.monomac.int@meta.data$percent.ribo))))
} else {
  message("already created integrated object, loading object...")
  all.monomac.int <- readRDS(paste0(dir.main, tempfile.name))
}

message("planning...")
plan("multicore", workers = 12, future.seed=T)
message("scaling...")
all.monomac.int <- ScaleData(all.monomac.int, features = all.monomac.int@assays$integrated@var.features,
                    vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ribo"))
message("pca'ing...")
all.monomac.int <- RunPCA(all.monomac.int, npcs = 20, verbose = FALSE)
message("reducing dimensions...")
all.monomac.int <- RunUMAP(all.monomac.int, reduction = "pca", dims = 1:20)
all.monomac.int <- RunTSNE(all.monomac.int, reduction = "pca", dims = 1:20)
message("clustering...")
all.monomac.int <- FindNeighbors(all.monomac.int, dims = 1:20)
all.monomac.int <- FindClusters(all.monomac.int, resolution = c(0.4, 0.6, 0.8))

message("saving...")
saveRDS(all.monomac.int, rdsOut)

message("done :)")

