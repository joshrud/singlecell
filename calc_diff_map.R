library(destiny)
library(Seurat)
library(ggplot2)

# args[1] is the input rds file 
# args[2] is the folder to output ggplots, but can be empty

args <- commandArgs(trailingOnly = T)

if (args[1] == "" || args[2] == "") {
  stop("arguments 1 and 2 are required...")
}

timeFmt <- function() {
  format(Sys.time(), "%y%m%d_%H%M%S")
}

dir.out <-  gsub("[^\\/]+$", "", x=args[1])
savefile.name = gsub(".*\\/", "", x=args[1])
tempfile.name <- gsub("\\.rds", "_withdiffmap.rds", savefile.name)

if (args[2] == "") {
  outdir <- paste0(dir.out, "/diffmap_outs/")
  message(paste0("no diffusion map directory supplied, creating directory: ", outdir))
  dir.create(outdir)
} else {
  outdir <- paste0(args[2], "/diffmap_outs/")
}

getColname <- function(x,y) {
  if (all(as.character(x) == as.character(y))) {
    return(TRUE)
  } else return(FALSE)
}

message(paste0(timeFmt(), " Reading input Seurat object..."))
a <- readRDS(args[1])
if (!any(sapply(a@meta.data, getColname, y=Idents(a)))) {
  message(unique(Idents(a)))
  stop("where are your idents coming from???")
}
a.metacol <- colnames(a@meta.data)[sapply(a@meta.data, getColname, y=Idents(a))]
message(paste0(timeFmt(), " Converting Seurat object to SingleCellExperiment..."))
message(paste0(timeFmt(), " Using ", a.metacol, " as meta data column..."))
a.sce <- as.SingleCellExperiment(a)
message(paste0(timeFmt(), " Calculating diffusion map..."))
dm <- DiffusionMap(t(a@assays$integrated@scale.data))
message(paste0(timeFmt(), " Saving diffusion map..."))
saveRDS(dm, paste0(dir.out, tempfile.name))

# create ggplots for all diffusion maps
dm.pd <- dm@eigenvectors
dm.pd <- as.data.frame(dm.pd, col.names = paste0("DC", 1:20))
dm.pd$labels <- a@meta.data[[a.metacol]] # this should be a user option somehow
all_coms <- combn(1:20, 2)
for (i in 1:ncol(all_coms)) {
  message(paste0("at ", timeFmt(), " DC.x = ", all_coms[1,i], " DC.y = ", all_coms[2,i]))
  p <- ggplot(dm.pd, aes_string(x=paste0("DC", all_coms[1,i]), y=paste0("DC", all_coms[2,i]),  #1 is first element in comb., 2 is 2nd
                                color="labels")) +
    geom_point(size = 0.5, stroke = 0, shape = 16) +
    guides(colour = guide_legend(override.aes = list(size=6)))
  theme_classic()
  pdf(paste0(outdir, "difmap_dc", all_coms[1,i], "_dc", all_coms[2,i], ".pdf"))
  print(p)
  dev.off()
}
message(paste0(timeFmt(), " Complete! "))

