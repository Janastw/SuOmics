# Quick Setup and Run Locally


This early workflow uses docker on a local machine. Be sure to have docker installed and build the image prior to running the workflow.
```
docker build -t seurat_qc -f Dockerfile.signac 
```

Place your samples from Cell Ranger (directories) in `data/` and then run
```
nextflow run main.nf
```
