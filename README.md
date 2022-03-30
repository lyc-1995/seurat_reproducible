# seurat_reproducible
Explore the reproducibility of Seurat pipeline
## Instructions
### Dependency
docker and singularity

### Install
`git clone https://github.com/lyc-1995/seurat_reproducible.git`

### Directories and files
`data/rna_obj.Rds`: Seurat object to input. Created from the thymus datasets in E-MTAB-8581. Not included in the github repo since it will exceed the maximum allowed size (2.00 GiB).

`R/run_std.R`: R script to perform the standard pipeline of Seurat

`R/modify_seurat.R`: Modified Seurat functions, which are re-imported into the namespace of Seurat.

`run_all.sh`: Main .sh script to run original Seurat.

`run_all_modify.sh`: Main .sh script to run modified Seurat.

`output/`: R Objects output from Seurat pipelines. The file names are formatted as `obj_${analysis}_${modified_or_not}_${hostname}_${pid}.Rds`

`nb/`: Jupyter notebooks for analysis results.

## Test for reproducibility
Basically, I use docker and singularity to make sure that the pipelines are run in the same environment. I also used `samba` to mount the same working directory onto two different devices (told by `hostname`). The objects output from each device were collected and compared. 

## Log
2022-03-29
* Currently only the standard pipeline (without clustering) was tested.
