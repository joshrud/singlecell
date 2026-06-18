# for running FindAllMarkers() in parallel
# Instructions for running: 
# 1. command is run on command line with 3 arguments: 
#    - arg1: the rds input object 
#    - arg2: output marker text file location 
#    - arg3: resolution to use (must have been calculated in loaded object) note: will prioritize "integrated_snn_res" over "RNA_snn_res"
# 2. example command line: Rscript parallelMarkersOnly.R /path/to/object/with/clusters/calculated.rds /path/to/marker/file/output.txt 0.2


library(Seurat)
library(future)
options(future.globals.maxSize= 10*1024^3) #10 GiB
args <- commandArgs(trailingOnly = T)

message("loading object...")  
if (args[1] == "" || args[2] == "" || args[3] == "") {
  stop("no argument specified...")
}
if (!file.exists(args[1])) {
  stop("file doesn't exist....")
}
if (!dir.exists(gsub("[^/]+$", "", args[2]))) {
  stop("can't write to dir: ", args[2])
}

so <- readRDS(args[1])

#check for the user-specified clustering resolution and set as Ident(so)
meta.cols <- colnames(so@meta.data)
if (suppressWarnings(!is.na(as.numeric(args[3])))) {
  print(paste0("selected column: ", args[3], " is a number, finding clust. res. "))
  if ("integrated" %in% names(so@assays)) {
    cluster.cols <- meta.cols[grep("integrated_snn_res", meta.cols)]
    cluster.cols <- gsub("integrated_snn_res\\.", "", cluster.cols)
    print(paste0("found cluster resolutions: ", paste0(cluster.cols, collapse=" ")))
    print(paste0("using resolution: ", args[3]))
    if (!(args[3] %in% cluster.cols)) {
      stop("not an option...")
    } else {
      Idents(so) <- so@meta.data[,paste0("integrated_snn_res.",args[3])]
    }
  } else if ("SCT" %in% names(so@assays)) {
    cluster.cols <- meta.cols[grep("SCT_snn_res", meta.cols)]
    cluster.cols <- gsub("SCT_snn_res\\.", "", cluster.cols)
    print(paste0("found cluster resolutions: ", paste0(cluster.cols, collapse=" ")))
    print(paste0("using resolution: ", args[3]))
    if (!(args[3] %in% cluster.cols)) {
      stop("not an option...")
    } else {
      Idents(so) <- so@meta.data[,paste0("SCT_snn_res.",args[3])]
    }
  } else {
    cluster.cols <- meta.cols[grep("RNA_snn_res", meta.cols)]
    cluster.cols <- gsub("RNA_snn_res\\.", "", cluster.cols)
    print(paste0("found cluster resolutions: ", paste0(cluster.cols, collapse=" ")))
    print(paste0("using resolution: ", args[3]))
    if (!(args[3] %in% cluster.cols)) {
      stop("not an option...")
    } else {
      Idents(so) <- so@meta.data[,paste0("RNA_snn_res.",args[3])]
    }
  }
} else {
  print(paste0("selected column: ", args[3], " is a character, finding matching cols "))
  cluster.cols <- meta.cols[grep(args[3], meta.cols)]
  if (length(cluster.cols) > 1) {
    print(paste0("input column matches more than one metadata column, choosing the first result: ", cluster.cols[1]))
    cluster.cols <- cluster.cols[1]
  }
  Idents(so) <- so@meta.data[,cluster.cols]
}

message("planning...")
plan("multiprocess", workers = 12)
if ("RNA" %in% names(so@assays)) {
  markers <- FindAllMarkers(so, only.pos = F, assay = "RNA", logfc.threshold=0.1)
} else {
  markers <- FindAllMarkers(so, only.pos = F, assay = "Spatial", logfc.threshold=0.1)
}
write.table(markers, args[2], quote = F, sep = "\t", col.names = NA)
message("done :)")

