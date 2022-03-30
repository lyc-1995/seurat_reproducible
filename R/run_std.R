
bin <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))

source(file.path(bin, "base.R"), chdir = TRUE)

args <- commandArgs(TRUE)
modify_seurat <- args[1]

analysis <- "std"
if ( !is.na( modify_seurat ) ) {
  source(modify_seurat, chdir = T)
  analysis <- paste0(analysis, "_modify")
}

input <- dirname(bin) %>% file.path("data")
outdir <- dirname(bin) %>% file.path("output")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

host_name <- Sys.info()["nodename"]
pid <- Sys.getpid()

# Let's shake it! ####################
obj <- readRDS(file.path(input, "rna_obj.Rds"))

obj <- NormalizeData(obj) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)

outfile <- paste0("obj_", analysis, "_", host_name, "_", pid, ".Rds")
saveRDS(obj, file = file.path(outdir, outfile))

