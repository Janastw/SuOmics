#!/bin/bash

docker run --rm -p 8787:8787 -e PASSWORD=yourpassword -v "$PWD":/home/rstudio/project seurat_qc rstudio