## Goal: build docker image that activates existing conda environment in the form of a yaml file

Once built, I push the docker image to DockerHub where I can pull it into a Snakemake worflow and run the container via Singularity. 

### Example dataset/script for testing: `Plot_SeuratData.R`

```
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunUMAP(pbmc, dims = 1:10)

pdf("./Tutorial_PBMC_UMAP.pdf")
DimPlot(pbmc, reduction = "umap")
dev.off()

rm(.RData)
```

### Snakefile structure (similar to make):

```
container: # docker or singularity container to run

rule all: # final all-encompassing rule
    input:
        # Some downstream output

rule a:
    input:
        # source data/input data
    output:
        # output data/visualization
    script:
        # script used to process input into output

```

## Testing conda environment with the Seurat docker image Dockerfile (three variations of the Dockerfile):

1. only_seurat: this is the original Dockerfile built by the developer of Seurat 
2. seurat_plus_conda: this is the original Dockerfile with some code to install conda and activate my conda env
3. only_conda: this Dockerfile only has code to install conda and activate my conda env (excludes the original Dockerfile content)


## `only_seurat` Dockerfile

```
FROM rocker/r-ver:4.1.0

# Set global R options
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV RETICULATE_MINICONDA_ENABLED=FALSE

# Install Seurat's system dependencies
RUN apt-get update
RUN apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    libboost-all-dev \
    libxml2-dev \
    openjdk-8-jdk \
    python3-dev \
    python3-pip \
    wget \
    git \
    libfftw3-dev \
    libgsl-dev

RUN apt-get install -y llvm-10

# Install UMAP
RUN LLVM_CONFIG=/usr/lib/llvm-10/bin/llvm-config pip3 install llvmlite
RUN pip3 install numpy
RUN pip3 install umap-learn

# Install FIt-SNE
RUN git clone --branch v1.2.1 https://github.com/KlugerLab/FIt-SNE.git
RUN g++ -std=c++11 -O3 FIt-SNE/src/sptree.cpp FIt-SNE/src/tsne.cpp FIt-SNE/src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm

# Install bioconductor dependencies & suggests
RUN R --no-echo --no-restore --no-save -e "install.packages('BiocManager')"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install(c('multtest', 'S4Vectors', 'SummarizedExperiment', 'SingleCellExperiment', 'MAST', 'DESeq2', 'BiocGenerics', 'GenomicRanges', 'IRanges', 'rtracklayer', 'monocle', 'Biobase', 'limma'))"

# Install CRAN suggests
RUN R --no-echo --no-restore --no-save -e "install.packages(c('VGAM', 'R.utils', 'metap', 'Rfast2', 'ape', 'enrichR', 'mixtools'))"

# Install hdf5r
RUN R --no-echo --no-restore --no-save -e "install.packages('hdf5r')"

# Install Seurat
RUN R --no-echo --no-restore --no-save -e "install.packages('remotes')"
RUN R --no-echo --no-restore --no-save -e "install.packages('Seurat')"

# Install SeuratDisk
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('mojaveazure/seurat-disk')"

CMD [ "R" ]
```

## `only_seurat` Snakefile

```
container: "docker://regnerm/only_seurat"

rule all:
    input:
        "Tutorial_PBMC_UMAP.pdf"

rule a:
    input:
        "filtered_gene_bc_matrices/hg19/"
    output:
        "Tutorial_PBMC_UMAP.pdf"
    script:
        "./Plot_SeuratData.R"
```

## `only_seurat` output
```
Building DAG of jobs...
Pulling singularity image docker://regnerm/only_seurat.
Using shell: /usr/bin/bash
Provided cores: 8
Rules claiming more threads will be scaled down.
Job stats:
job      count    min threads    max threads
-----  -------  -------------  -------------
a            1              1              1
all          1              1              1
total        2              1              1

Select jobs to execute...

[Thu Jan  6 11:29:09 2022]
rule a:
    input: filtered_gene_bc_matrices/hg19
    output: Tutorial_PBMC_UMAP.pdf
    jobid: 1
    resources: tmpdir=/datastore/scratch/users/regnerm

Activating singularity image /home/regnerm/only_seurat/.snakemake/singularity/d70de55de5d67bbe976ee1b7dcd14382.simg

Attaching package: ‘dplyr’

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

Attaching SeuratObject
Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
^M  |                                                                            ^M  |                                                                      |   0%^M  |                                                                            ^M  |=====                                                                 |   7%^M  |                                                                            ^M  |==========                                                            |  14%^M  |                                                                            ^M  |===============                                                       |  21%^M  |                                                                            ^M  |====================                                                  |  29%^M  |                                                                            ^M  |=========================                                             |  36%^M  |                                                                            ^M  |==============================                                        |  43%^M  |                                                                            ^M  |===================================                                   |  50%^M  |                                                                            ^M  |========================================                              |  57%^M  |                                                                            ^M  |=============================================                         |  64%^M  |                                                                            ^M  |==================================================                    |  71%^M  |                                                                            ^M  |=======================================================               |  79%^M  |                                                                            ^M  |============================================================          |  86%^M  |                                                                            ^M  |=================================================================     |  93%^M  |                                                                            ^M  |======================================================================| 100%
PC_ 1
Positive:  CST3, TYROBP, LST1, AIF1, FTL, FTH1, LYZ, FCN1, S100A9, TYMP
           FCER1G, CFD, LGALS1, S100A8, CTSS, LGALS2, SERPINA1, IFITM3, SPI1, CFP
           PSAP, IFI30, SAT1, COTL1, S100A11, NPC2, GRN, LGALS3, GSTP1, PYCARD
Negative:  MALAT1, LTB, IL32, IL7R, CD2, B2M, ACAP1, CD27, STK17A, CTSW
           CD247, GIMAP5, AQP3, CCL5, SELL, TRAF3IP3, GZMA, MAL, CST7, ITM2A
           MYC, GIMAP7, HOPX, BEX2, LDLRAP1, GZMK, ETS1, ZAP70, TNFAIP8, RIC3
PC_ 2
Positive:  CD79A, MS4A1, TCL1A, HLA-DQA1, HLA-DQB1, HLA-DRA, LINC00926, CD79B, HLA-DRB1, CD74
           HLA-DMA, HLA-DPB1, HLA-DQA2, CD37, HLA-DRB5, HLA-DMB, HLA-DPA1, FCRLA, HVCN1, LTB
           BLNK, P2RX5, IGLL5, IRF8, SWAP70, ARHGAP24, FCGR2B, SMIM14, PPP1R14A, C16orf74
Negative:  NKG7, PRF1, CST7, GZMB, GZMA, FGFBP2, CTSW, GNLY, B2M, SPON2
           CCL4, GZMH, FCGR3A, CCL5, CD247, XCL2, CLIC3, AKR1C3, SRGN, HOPX
           TTC38, APMAP, CTSC, S100A4, IGFBP7, ANXA1, ID2, IL32, XCL1, RHOC
PC_ 3
Positive:  HLA-DQA1, CD79A, CD79B, HLA-DQB1, HLA-DPB1, HLA-DPA1, CD74, MS4A1, HLA-DRB1, HLA-DRA
           HLA-DRB5, HLA-DQA2, TCL1A, LINC00926, HLA-DMB, HLA-DMA, CD37, HVCN1, FCRLA, IRF8
           PLAC8, BLNK, MALAT1, SMIM14, PLD4, P2RX5, IGLL5, LAT2, SWAP70, FCGR2B
Negative:  PPBP, PF4, SDPR, SPARC, GNG11, NRGN, GP9, RGS18, TUBB1, CLU
           HIST1H2AC, AP001189.4, ITGA2B, CD9, TMEM40, PTCRA, CA2, ACRBP, MMD, TREML1
           NGFRAP1, F13A1, SEPT5, RUFY1, TSC22D1, MPP1, CMTM5, RP11-367G6.3, MYL9, GP1BA
PC_ 4
Positive:  HLA-DQA1, CD79B, CD79A, MS4A1, HLA-DQB1, CD74, HIST1H2AC, HLA-DPB1, PF4, SDPR
           TCL1A, HLA-DRB1, HLA-DPA1, HLA-DQA2, PPBP, HLA-DRA, LINC00926, GNG11, SPARC, HLA-DRB5
           GP9, AP001189.4, CA2, PTCRA, CD9, NRGN, RGS18, CLU, TUBB1, GZMB
Negative:  VIM, IL7R, S100A6, IL32, S100A8, S100A4, GIMAP7, S100A10, S100A9, MAL
           AQP3, CD2, CD14, FYB, LGALS2, GIMAP4, ANXA1, CD27, FCN1, RBP7
           LYZ, S100A11, GIMAP5, MS4A6A, S100A12, FOLR3, TRABD2A, AIF1, IL8, IFI6
PC_ 5
Positive:  GZMB, NKG7, S100A8, FGFBP2, GNLY, CCL4, CST7, PRF1, GZMA, SPON2
           GZMH, S100A9, LGALS2, CCL3, CTSW, XCL2, CD14, CLIC3, S100A12, RBP7
           CCL5, MS4A6A, GSTP1, FOLR3, IGFBP7, TYROBP, TTC38, AKR1C3, XCL1, HOPX
Negative:  LTB, IL7R, CKB, VIM, MS4A7, AQP3, CYTIP, RP11-290F20.3, SIGLEC10, HMOX1
           LILRB2, PTGES3, MAL, CD27, HN1, CD2, GDI2, CORO1B, ANXA5, TUBA1B
           FAM110A, ATP1A1, TRADD, PPA1, CCDC109B, ABRACL, CTD-2006K23.1, WARS, VMO1, FYB
Computing nearest neighbor graph
Computing SNN
Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

Number of nodes: 2638
Number of edges: 95927

Running Louvain algorithm...
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Maximum modularity in 10 random starts: 0.8728
Number of communities: 9
Elapsed time: 0 seconds
Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
This message will be shown once per session
11:29:33 UMAP embedding parameters a = 0.9922 b = 1.112
11:29:33 Read 2638 rows and found 10 numeric columns
11:29:33 Using Annoy for neighbor search, n_neighbors = 30
11:29:33 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
11:29:33 Writing NN index file to temp file /tmp/Rtmpoxy0sY/file17227e6db72e
11:29:33 Searching Annoy index using 1 thread, search_k = 3000
11:29:35 Annoy recall = 100%
11:29:35 Commencing smooth kNN distance calibration using 1 thread
11:29:36 Initializing from normalized Laplacian + noise
11:29:36 Commencing optimization for 500 epochs, with 105140 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
11:29:41 Optimization finished
null device
          1
Warning message:
In rm(.RData) : object '.RData' not found
[Thu Jan  6 11:29:42 2022]
Finished job 1.
1 of 2 steps (50%) done
Select jobs to execute...

[Thu Jan  6 11:29:42 2022]
localrule all:
    input: Tutorial_PBMC_UMAP.pdf
    jobid: 0
    resources: tmpdir=/datastore/scratch/users/regnerm

[Thu Jan  6 11:29:42 2022]
Finished job 0.
2 of 2 steps (100%) done
Complete log: /home/regnerm/only_seurat/.snakemake/log/2022-01-06T112846.674032.snakemake.log

```

## `seurat_plus_conda` Dockerfile

```
FROM rocker/r-ver:4.1.0

# Set global R options
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV RETICULATE_MINICONDA_ENABLED=FALSE

# Install Seurat's system dependencies
RUN apt-get update
RUN apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    libboost-all-dev \
    libxml2-dev \
    openjdk-8-jdk \
    python3-dev \
    python3-pip \
    wget \
    git \
    libfftw3-dev \
    libgsl-dev

RUN apt-get install -y llvm-10

# Install UMAP
RUN LLVM_CONFIG=/usr/lib/llvm-10/bin/llvm-config pip3 install llvmlite
RUN pip3 install numpy
RUN pip3 install umap-learn

# Install FIt-SNE
RUN git clone --branch v1.2.1 https://github.com/KlugerLab/FIt-SNE.git
RUN g++ -std=c++11 -O3 FIt-SNE/src/sptree.cpp FIt-SNE/src/tsne.cpp FIt-SNE/src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm

# Install bioconductor dependencies & suggests
RUN R --no-echo --no-restore --no-save -e "install.packages('BiocManager')"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install(c('multtest', 'S4Vectors', 'SummarizedExperiment', 'SingleCellExperiment', 'MAST', 'DESeq2', 'BiocGenerics', 'GenomicRanges', 'IRanges', 'rtracklayer', 'monocle', 'Biobase', 'limma'))"

# Install CRAN suggests
RUN R --no-echo --no-restore --no-save -e "install.packages(c('VGAM', 'R.utils', 'metap', 'Rfast2', 'ape', 'enrichR', 'mixtools'))"

# Install hdf5r
RUN R --no-echo --no-restore --no-save -e "install.packages('hdf5r')"

# Install Seurat
RUN R --no-echo --no-restore --no-save -e "install.packages('remotes')"
RUN R --no-echo --no-restore --no-save -e "install.packages('Seurat')"

# Install SeuratDisk
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('mojaveazure/seurat-disk')"

#CMD [ "R" ]

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
     /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# Create the environment:
COPY scBreast_2021_v2.yml .
RUN conda env create -f scBreast_2021_v2.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "scBreast_2021_v2", "/bin/bash", "-c"]

# The code to run when container is started:
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "scBreast_2021_v2"]
```


## `seurat_plus_conda` Snakefile

```
container: "docker://regnerm/seurat_plus_conda"

rule all:
    input:
        "Tutorial_PBMC_UMAP.pdf"

rule a:
    input:
        "filtered_gene_bc_matrices/hg19/"
    output:
        "Tutorial_PBMC_UMAP.pdf"
    script:
        "./Plot_SeuratData.R"
```


## `seurat_plus_conda` output

```
Building DAG of jobs...
Pulling singularity image docker://regnerm/seurat_plus_conda.
Using shell: /usr/bin/bash
Provided cores: 8
Rules claiming more threads will be scaled down.
Job stats:
job      count    min threads    max threads
-----  -------  -------------  -------------
a            1              1              1
all          1              1              1
total        2              1              1

Select jobs to execute...

[Thu Jan  6 12:13:45 2022]
rule a:
    input: filtered_gene_bc_matrices/hg19
    output: Tutorial_PBMC_UMAP.pdf
    jobid: 1
    resources: tmpdir=/datastore/scratch/users/regnerm

Activating singularity image /home/regnerm/seurat_plus_conda/.snakemake/singularity/48a9e3532533b6057ebb485122652b93.simg

Attaching package: ‘dplyr’

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

Attaching SeuratObject
Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
^M  |                                                                            ^M  |                                                                      |   0%^M  |                                                                            ^M  |=====                                                                 |   7%^M  |                                                                            ^M  |==========                                                            |  14%^M  |                                                                            ^M  |===============                                                       |  21%^M  |                                                                            ^M  |====================                                                  |  29%^M  |                                                                            ^M  |=========================                                             |  36%^M  |                                                                            ^M  |==============================                                        |  43%^M  |                                                                            ^M  |===================================                                   |  50%^M  |                                                                            ^M  |========================================                              |  57%^M  |                                                                            ^M  |=============================================                         |  64%^M  |                                                                            ^M  |==================================================                    |  71%^M  |                                                                            ^M  |=======================================================               |  79%^M  |                                                                            ^M  |============================================================          |  86%^M  |                                                                            ^M  |=================================================================     |  93%^M  |                                                                            ^M  |======================================================================| 100%
PC_ 1
Positive:  CST3, TYROBP, LST1, AIF1, FTL, FTH1, LYZ, FCN1, S100A9, TYMP
           FCER1G, CFD, LGALS1, S100A8, CTSS, LGALS2, SERPINA1, IFITM3, SPI1, CFP
           PSAP, IFI30, SAT1, COTL1, S100A11, NPC2, GRN, LGALS3, GSTP1, PYCARD
Negative:  MALAT1, LTB, IL32, IL7R, CD2, B2M, ACAP1, CD27, STK17A, CTSW
           CD247, GIMAP5, AQP3, CCL5, SELL, TRAF3IP3, GZMA, MAL, CST7, ITM2A
           MYC, GIMAP7, HOPX, BEX2, LDLRAP1, GZMK, ETS1, ZAP70, TNFAIP8, RIC3
PC_ 2
Positive:  CD79A, MS4A1, TCL1A, HLA-DQA1, HLA-DQB1, HLA-DRA, LINC00926, CD79B, HLA-DRB1, CD74
           HLA-DMA, HLA-DPB1, HLA-DQA2, CD37, HLA-DRB5, HLA-DMB, HLA-DPA1, FCRLA, HVCN1, LTB
           BLNK, P2RX5, IGLL5, IRF8, SWAP70, ARHGAP24, FCGR2B, SMIM14, PPP1R14A, C16orf74
Negative:  NKG7, PRF1, CST7, GZMB, GZMA, FGFBP2, CTSW, GNLY, B2M, SPON2
           CCL4, GZMH, FCGR3A, CCL5, CD247, XCL2, CLIC3, AKR1C3, SRGN, HOPX
           TTC38, APMAP, CTSC, S100A4, IGFBP7, ANXA1, ID2, IL32, XCL1, RHOC
PC_ 3
Positive:  HLA-DQA1, CD79A, CD79B, HLA-DQB1, HLA-DPB1, HLA-DPA1, CD74, MS4A1, HLA-DRB1, HLA-DRA
           HLA-DRB5, HLA-DQA2, TCL1A, LINC00926, HLA-DMB, HLA-DMA, CD37, HVCN1, FCRLA, IRF8
           PLAC8, BLNK, MALAT1, SMIM14, PLD4, P2RX5, IGLL5, LAT2, SWAP70, FCGR2B
Negative:  PPBP, PF4, SDPR, SPARC, GNG11, NRGN, GP9, RGS18, TUBB1, CLU
           HIST1H2AC, AP001189.4, ITGA2B, CD9, TMEM40, PTCRA, CA2, ACRBP, MMD, TREML1
           NGFRAP1, F13A1, SEPT5, RUFY1, TSC22D1, MPP1, CMTM5, RP11-367G6.3, MYL9, GP1BA
PC_ 4
Positive:  HLA-DQA1, CD79B, CD79A, MS4A1, HLA-DQB1, CD74, HIST1H2AC, HLA-DPB1, PF4, SDPR
           TCL1A, HLA-DRB1, HLA-DPA1, HLA-DQA2, PPBP, HLA-DRA, LINC00926, GNG11, SPARC, HLA-DRB5
           GP9, AP001189.4, CA2, PTCRA, CD9, NRGN, RGS18, CLU, TUBB1, GZMB
Negative:  VIM, IL7R, S100A6, IL32, S100A8, S100A4, GIMAP7, S100A10, S100A9, MAL
           AQP3, CD2, CD14, FYB, LGALS2, GIMAP4, ANXA1, CD27, FCN1, RBP7
           LYZ, S100A11, GIMAP5, MS4A6A, S100A12, FOLR3, TRABD2A, AIF1, IL8, IFI6
PC_ 5
Positive:  GZMB, NKG7, S100A8, FGFBP2, GNLY, CCL4, CST7, PRF1, GZMA, SPON2
           GZMH, S100A9, LGALS2, CCL3, CTSW, XCL2, CD14, CLIC3, S100A12, RBP7
           CCL5, MS4A6A, GSTP1, FOLR3, IGFBP7, TYROBP, TTC38, AKR1C3, XCL1, HOPX
Negative:  LTB, IL7R, CKB, VIM, MS4A7, AQP3, CYTIP, RP11-290F20.3, SIGLEC10, HMOX1
           LILRB2, PTGES3, MAL, CD27, HN1, CD2, GDI2, CORO1B, ANXA5, TUBA1B
           FAM110A, ATP1A1, TRADD, PPA1, CCDC109B, ABRACL, CTD-2006K23.1, WARS, VMO1, FYB
Computing nearest neighbor graph
Computing SNN
Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

Number of nodes: 2638
Number of edges: 95927

Running Louvain algorithm...
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Maximum modularity in 10 random starts: 0.8728
Number of communities: 9
Elapsed time: 0 seconds
Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
This message will be shown once per session
12:14:09 UMAP embedding parameters a = 0.9922 b = 1.112
12:14:09 Read 2638 rows and found 10 numeric columns
12:14:09 Using Annoy for neighbor search, n_neighbors = 30
12:14:09 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
12:14:10 Writing NN index file to temp file /tmp/Rtmppwylf2/file285629674755
12:14:10 Searching Annoy index using 1 thread, search_k = 3000
12:14:11 Annoy recall = 100%
12:14:12 Commencing smooth kNN distance calibration using 1 thread
12:14:12 Initializing from normalized Laplacian + noise
12:14:13 Commencing optimization for 500 epochs, with 105140 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
12:14:18 Optimization finished
null device
          1
Warning message:
In rm(.RData) : object '.RData' not found
[Thu Jan  6 12:14:19 2022]
Finished job 1.
1 of 2 steps (50%) done
Select jobs to execute...

[Thu Jan  6 12:14:19 2022]
localrule all:
    input: Tutorial_PBMC_UMAP.pdf
    jobid: 0
    resources: tmpdir=/datastore/scratch/users/regnerm

[Thu Jan  6 12:14:19 2022]
Finished job 0.
2 of 2 steps (100%) done
Complete log: /home/regnerm/seurat_plus_conda/.snakemake/log/2022-01-06T113641.904575.snakemake.log

```


## `only_conda` Dockerfile
```
FROM rocker/r-ver:4.1.0


# Install python and wget
RUN apt-get update
RUN apt-get install -y \
    python3-dev \
    python3-pip \
    wget

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
     /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# Create the environment:
COPY scBreast_2021_v2.yml .
RUN conda env create -f scBreast_2021_v2.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "scBreast_2021_v2", "/bin/bash", "-c"]

# The code to run when container is started:
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "scBreast_2021_v2"]

```
## `only_conda` Snakefile

```
container: "docker://regnerm/only_conda"

rule all:
    input:
        "Tutorial_PBMC_UMAP.pdf"

rule a:
    input:
        "filtered_gene_bc_matrices/hg19/"
    output:
        "Tutorial_PBMC_UMAP.pdf"
    script:
        "./Plot_SeuratData.R"
```


## `only_conda` output

```
Building DAG of jobs...
Pulling singularity image docker://regnerm/only_conda.
Using shell: /usr/bin/bash
Provided cores: 8
Rules claiming more threads will be scaled down.
Job stats:
job      count    min threads    max threads
-----  -------  -------------  -------------
a            1              1              1
all          1              1              1
total        2              1              1

Select jobs to execute...

[Thu Jan  6 13:21:08 2022]
rule a:
    input: filtered_gene_bc_matrices/hg19
    output: Tutorial_PBMC_UMAP.pdf
    jobid: 1
    resources: tmpdir=/datastore/scratch/users/regnerm

Activating singularity image /home/regnerm/only_conda/.snakemake/singularity/6ccfa8449f3aee8c61bab3a1c829e98d.simg
Error in library(dplyr) : there is no package called ‘dplyr’
Execution halted
[Thu Jan  6 13:21:10 2022]
Error in rule a:
    jobid: 1
    output: Tutorial_PBMC_UMAP.pdf

RuleException:
CalledProcessError in line 13 of /home/regnerm/only_conda/Snakefile:
Command ' singularity  exec --home /home/regnerm/only_conda  /home/regnerm/only_conda/.snakemake/singularity/6ccfa8449f3aee8c61bab3a1c829e98d.simg bash -c 'set -euo pipefail;  Rscript --vanilla /home/regnerm/only_conda/.snakemake/scripts/tmpb5q5n5_a.Plot_SeuratData.R'' returned non-zero exit status 1.
  File "/home/regnerm/only_conda/Snakefile", line 13, in __rule_a
  File "/home/regnerm/anaconda3/envs/snakemake/lib/python3.7/concurrent/futures/thread.py", line 57, in run
Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
Complete log: /home/regnerm/only_conda/.snakemake/log/2022-01-06T130026.696479.snakemake.log
```


## Conclusions:

only_seurat runs to completion. seurat_plus_conda runs to completion, but only_conda fails to find the R package dplyr.

This means that the conda env is not being correctly loaded/built in the Docker container. 
