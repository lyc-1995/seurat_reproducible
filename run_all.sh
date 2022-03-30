#!/bin/bash -e

wd=$(dirname $(readlink -f $0))

IMAGE=lyc1995/bioinfo:1.0.1

singularity run -B ${wd}:${wd} docker://${IMAGE} Rscript ${wd}/R/run_std.R
