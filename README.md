# Spatial Transcriptomic Tools (STtools)

STtools is a software package that is designed to process spatial
transciriptomics (ST) data from various platforms including
[Seq-Scope](https://www.cell.com/cell/fulltext/S0092-8674(21)00627-9?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421006279%3Fshowall%3Dtrue), 
[SlideSeq](https://www.cell.com/cell/fulltext/S0092-8674(21)00627-9?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421006279%3Fshowall%3Dtrue), 
and [VISIUM](https://www.nature.com/articles/s42003-020-01247-y). 
The STTools pipeline includes preprocessing of raw sequence reads, 
alignment, collapsing barcodes into grids, clustering cell types, and
high-resolution analysis with sliding window strategy.
STTools leverages many existing software tools for single-cell and
spatial transcriptomic analysis, such as 
[STARSolo](https://github.com/alexdobin/STAR),
[Seurat](https://satijalab.org/seurat/articles/spatial_vignette.html),
[BayesSpace](https://www.nature.com/articles/s41587-021-00935-2), and
[Seqtk](https://github.com/lh3/seqtk).

## Getting Started

We recommend running STTools in a linux operating system (e.g. Ubuntu
18.04). See [Installtion](#installation) for required software tools
to run STTools.

```sh
## clone the repository
git clone https://github.com/jyxi7676/STtools.git
cd STtools
## install required python packages
python -m pip install -r requirements.txt
## download example data and decompress
gdown https://drive.google.com/uc?id=1e0u57Yu_fVKFvs-UA7WYfj-vgm8Nd2y4
unzip STtools_example_data.zip 
## create output directory and set environment variables
mkdir out
export STHOME=$(pwd)
export STDATA=$STHOME/STtools_example_data ## directory containing data
export STOUT=$STHOME/out             ## output directory
export SEQTKPATH=/path/to/seqtk/bin  ## path that contains seqtk binary
export STARPATH=/path/to/STAR/bin    ## path that contains STAR binary
export GENOMEINDEX=/path/to/STAR/index ## path that contains STAR index
## UNCOMMENT if you need to build STAR index yourself for the example data,
## mkdir -p $STHOME/STtools_example_data/geneIndex/STARIndex
## $STARPATH/STAR --runThreadN 6 --runMode genomeGenerate --genomeDir $STHOME/STtools_example_data/geneIndex/STARIndex \
##     --genomeFastaFiles $STHOME/STtools_example_data/geneIndex/mm10.fasta \
##     --sjdbGTFfile $STHOME/STtools_example_data/geneIndex/mm10.gtf --sjdbOverhang 99
## export GENOMEINDEX=$STDATA/geneIndex/STARIndex/
## 
## Run STTools - step 1 to 7
python3 $STHOME/sttools.py --run-all --STtools $STHOME \
  --first-fq $STDATA/step1_extractCoordinates/liver-MiSeq-tile2106-sub-R1.fastq.gz \
  --second-fq1 $STDATA/step3_align/liver_tile2106_sub_R1.fastq.gz \
  --second-fq2 $STDATA/step3_align/liver_tile2106_sub_R2.fastq.gz \
  --outdir $STOUT --genome $GENOMEINDEX --star-path $STARPATH --seqtk-path $SEQTKPATH \
  --seqscope1st 'HiSeq' --clustering False --lane-tiles 1_2106 \
  --binsize 300 --window 150 -l 20 -o 'Sample' -c 2
```

STtools package have flexible options for the user to run **all
steps**, **specificn steps**, or **consecutive steps**. 
Several examples from various scenarios are given below for illustratrion. 
* [Running all steps (1-7)](./doc/readme1.md)
* [Running consecutive steps](./doc/readme2.md)
* [Running specific step](./doc/readme3.md)

## Overview of STtools

This image below illustrates the overall workflow for STtools. 

<p align="center">
    <img src="STtools_workflow.png" width="1550" height="700" />
</p>

There are 7 steps in total. 
Each step takes input from either the raw data or outputs of the
previous steps. Please see a brief explanation on each step:

* **Step 1** takes `fastq.gz` files as input and output spatial coordinates `.txt` files and whitelist used for `STARsolo` alignemnt in the current working directory.
* **Step 2** takes barcode info, and spatial coordinates file to generate a barcode/HDMI density plot which can be compared with HE images for an estimation of tissue boundary
* **Step 3** takes valid barcodes `whitelist.txt`, 2nd-seq `fastq.gz`
  files, and the STAR indices of reference genome as input to run
  `STARsolo` alignment; this step outputs digital expression matrix (DGE).
* **Step 4** takes DGE from **Step 3** and output `Seurat` object with collapsed DGE of simple square grids.
* **Step 5** takes DGE from **Step 3** and output Seurat object with collapsed DGE of square grids from sliding window strategy
* **Step 6** takes in RDS file from **Step 4** and **Step 5** as input and performs dimension reduction, clustering and conducts refernece mapping with simple square grids as query
* **Step 7** takes DGE (Velocyto-format) from **Step 3** and generate subcellular plots showing pattern of spliced/unspliced reads

## Installation
Linux operatin system is necessary to run STtools package. You also need to install the following software tools and librares/modules before using this package.
* STAR>=2.7.5c (Click for instructions to install [STAR](https://github.com/alexdobin/STAR))
* seqtk (Click for instructions to install [seqtk](https://github.com/lh3/seqtk))
* R>=4.0.0 (STtools will install packages automatically if not installed. Please refer to the  list of [packages](./doc/Rpackages))
* Python>= 3.0 (STtools will install modules automatically if not installed, refer to the list of [modules](./doc/PythonModules))
* perl(Click for instructions for installing [perl](https://learn.perl.org/installing/unix_linux.html) )
* pigz(Click for instructions for installing [pigz](https://zlib.net/pigz/))


To install **STtools**, please run:
```
git clone https://github.com/jyxi7676/STtools.git
```


## Example Data
* SeqScope exmaple data for each step can be found at [example data 1](https://drive.google.com/file/d/1e0u57Yu_fVKFvs-UA7WYfj-vgm8Nd2y4/view?usp=sharing), please download the zip files. For each step, the example input data is stored in the corresponding subdirectories. 
* VISIUM digital expresstion data and spatial coordinates are available at [example data 2](https://drive.google.com/drive/folders/130ENNRBEi7kCOXDnGZlHUnuf4CD3_JEI?usp=sharing)
* SlideSeq digital expression data and spatial coordinates are avaialbel at [example data 3](https://drive.google.com/drive/folders/1IktkJgDLnYS0fcW65xgHC04S-Mr8ciwf?usp=sharing)

## Input Data Format
Please refer to [data formats](./doc/fileformats.md) for an illustration of required input data format for each step.

## External links
Here are some useful external links:
* To generate gene index for STARsolo alignment: https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html
* Multimodal reference mapping: https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
* Incoporate transgenes to alignment: Please modify the gtf and fasta files according to https://github.com/igordot/genomics/blob/master/workflows/ref-genome-gfp.md before generating  genome  index in STAR.
