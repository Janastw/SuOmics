FROM rocker/rstudio:4.4.0

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libhdf5-dev \
    liblzma-dev \
    libbz2-dev \
    build-essential \
    libglpk-dev \
    libgsl0-dev \
    libudunits2-dev \
    zlib1g-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libjpeg-dev \
    libtiff5-dev \
    libpng-dev \
    tabix \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('remotes', 'patchwork', 'ggplot2', 'data.table', 'vctrs', 'devtools', 'matrixStats', 'Seurat', 'SeuratObject'), repos='https://cloud.r-project.org')"

RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

RUN R -e "BiocManager::install(c( \
  'Rsamtools', 'GenomicRanges', 'GenomeInfoDb', \
  'EnsDb.Mmusculus.v79', 'Biostrings', 'biovizBase' \
), ask=FALSE, update=TRUE)"

RUN R -e "install.packages('hdf5r')"

RUN R -e "BiocManager::install(c('chromVAR', 'JASPAR2020', 'motifmatchr', 'BSgenome.Mmusculis.ECSC.mm10', 'glmGamPoi', 'BSgenome.Mmusculus.UCSC.mm10'))"
RUN R -e "remotes::install_github('timoast/signac', dependencies = TRUE, upgrade = 'always', force = TRUE)"
RUN R -e "BiocManager::install(c('ensembldb', 'AnnotationHub', 'EnsDb.Hsapiens.v86'))"