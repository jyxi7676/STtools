
# User Manual
## Running specific step in STtools
The user need to use the option --run-step if interested in running one step only per command .   We will give a detailed illustration about how to run STtools per step with several examples. Please modify your input files according to  [the link](./doc/fileformats.md) 
## Step 1
Step 1 aims to extracts spatial coordinates and whitelist info from the sequenced raw FASTQ.gz file. It assumes the user has the SeqScope data format (see xxxx for stucture, cite seqscope paper), where spatial  realated information such as HDMI/Barcode, lane, tile, X and Y coordinates can be retrieved from 1st-Seq and transcriptomic information can be retrieved from 2nd-Seq. 
### *Input*
  The user need to provide some of the following files in order for STtools to run Step 1. 
  *   --first--fq: Path to 1st-Seq FASTQ.gz file. STtools takes in fastq.gz files with SeqScope sequence design structure. If the barcode/UMI/randomer location is different, please see this link xxxxx for an example to make the inputs compatible to STtools package. **Required**. 
  *   --second-fq1 : Path to 2nd-Seq Read 1 FASTQ.gz file. If the barcode/UMI/randomer location is different, please see this link for an example to make it compatible to STtools package.**Required**. 
  *   --hdmilength: An integer of the length of HDMI/Barcode.(modify to take any number below 32?). By default hdmilength=20.
  *   --STtools: Path to the STtools package. If not given, using current working directory (add this to pkg)
  *   --outdir: Path to output files. If not given, using current working directory
 ### *Code*
 ```sh
 ## $STHOME indicates the path to the directory of STtools repository
 export STHOME=/path/to/STtools 
 ## $STDATA indicates the directory containing the example input files
 export STDATA=/path/to/data
 ## $STOUT indicates the directory containing the example output files.
 export STOUT=/path/to/outdir
 
 python3 $STHOME/sttools.py --run-steps 1 --first-fq $STDATA/liver-MiSeq-tile2106-sub-R1.fastq.gz --second-fq1 $STDATA/liver-HiSeq-tile2106-sub-R1.fastq.gz --STtools $STHOME  -l 20 --outdir $STOUT
 ```
 ### *Output*
 This step output two useful files in the current working directory, and are taken as input the next steps.
 * spatialcoordinates.txt 
 * whitelist.txt
 * HDMI_SeqScope_2nd.txt
 * summary_step1.txt
 
## Step 2
Step2 visualize the barcode/HDMI density discovery plot, with which the user are able to compare with HE images to estimate the tissue boundary. The alignment is done manually and automatic alignment is under development. 
### *Input*
  * --STtools: Path to the STtools package. If not given, using current working directory (add this to pkg)
  * --hdmi2ndSeq: : txt file with the HDMIs from the --second-fq1.
  * --spatial: txt file with HDMI and spatial information from --first-fq: HDMI, lane, tile, X and Y.
  * --maxScale: max scale value for color bar, if not specified, using default value in matplotlib.
  * --outdir: Path to output files. If not given, using current working directory

### *Code*
 ```sh
  ## $STHOME indicates the path to the directory of STtools repository
  export STHOME=/path/to/STtools
  ## $STDATA indicates the directory containing the example input files
  export STDATA=/path/to/data
  ## $STOUT indicates the directory containing the example output files.
  export STOUT=/path/to/outdir
  python3 $STHOME/sttools.py --run-steps 2 --STtools $STHOME --spatial $STDATA/spatialcoordinates.txt --hdmi2ndSeq $STDATA/HDMI_SeqScope_2nd.txt --outdir $STOUT
  ```
### *Output*
tile_lane*.png

## Step 3
Step 3 aligns the data the reference genome using STARsolo software and output digital expression matrix under Gene,GeneFull, and Velocyto options (see link xxxxx).
Before running step3, it is the required the user to generate the genome reference following this link  xxxxx.
### *Input*
 * --second-fq1: Path to 2nd-Seq FASTQ.gz file of read 1. Required.
 * --second-fq2: Path to 2nd-Seq FASTQ.gz file of read 2. Required.
 * --whitelist: Txt file of the whitelist, if not given, searching whitelist.txt in the current folder.
 * --genome or -g: STAR geneome index. Required
 * --STtools: Path to the STtools package. If not given, using currently working directory (add this to pkg)
 * -o: Output prefix of alignment, if not given, set to 'Sample'
 * --sesqtk-path: Path to seqtk executable. Required
 * --star-path: Path to STAR executable. Required
 * --outdir: Path to output files. If not given, using current working directory

 
### *Code*
``` sh
 ## $STHOME indicates the path to the directory of STtools repository
 export STHOME=/path/to/STtools
 ## $STDATA indicates the directory containing the example input files
 export STDATA=/path/to/data
 ## $STOUT indicates the directory containing the example output files.
 export STOUT=/path/to/outdir
 ## $GENEINDEX indicates the path to genome index for STARsolo alignment
 export GENOMEINDEX=/path/to/genomeIndex
 ## $SEQTKPATH indicates the path to the seqtk executive
 export SEQTKPATH=/path/to/seqtk/executive
 ## $STARPATH indicates the path to the STAR executive
 export STARPATH=/path/to/star/executive
 
python3 $STHOME/sttools.py --run-steps 3 --STtools $STHOME --whitelist $STDATA/whitelist.txt --second-fq1 $STDATA/liver_tile2106_sub_R1.fastq.gz --second-fq2 $STDATA/liver_tile2106_sub_R2.fastq.gz --outdir $STOUT --genome $GENOMEINDEX --star-path $STARPATH --seqtk-path $SEQTKPATH

```
### *Output*
This step outputs folders with STARsolo summary statistics, bam file, DGE, etc.
* SampleSolo.out/
* summary_step3.txt
  
## Step 4
Step 4 bins the DGE into simple square gridded data and collapses the reads counts within each grid. A new DGE is created and spatial information is updated with the center of each bin, and is saved as an RDS file for output in the current working directory.
### *input*
* --STtools: Path to STtools package. If not given, the current working directory is used.(add this to pkg).
* --tiles: Tiles that the user is insterested in. Multiple tile numbers need to be separated by comma. For example: --tiles 2106,2107,2108. Required.
* --binsize: The size of the square grids side. By default, it is set to 300 units.
* --DGEdir: Path to the digital expression matrix. If not given, use the path to DGE from previous steps.
* --spatial: Path to the txt file of spatial coordinates. If not given, use the path to spatialCoordinates.txt from previous steps.
* --nrow: number of rows when generating the super tile. If not give, we assume there is only one tile, and nrow=1.
* --ncol: number of cols when generating the super tile. If not give, we assume there is only one tile, and ncol=1.
* --outdir: Path to output files. If not given, using current working directory
### *code*
```sh
## $STHOME indicates the path to the directory of STtools repository
export STHOME=/path/to/STtools
## $STDATA indicates the directory containing the example input files
export STDATA=/path/to/data
## $STOUT indicates the directory containing the example output files
export STOUT=/path/to/outdir
## $STOUT indicates the path to digital expression matrix(DGE)
export STDGE=/path/to/DGE/
python3  $STHOME/sttools.py --run-steps 4 --STtools $STHOME --spatial $STDATA/spatialcoordinates.txt   --outdir $STOUT --tiles 2106 --binsize 300  -l 20 --DGEdir $STDGE

```
### *output*
* SimpleSquareGrids.RDS
* summary_step4.txt

## Step 5
Step 5 bins the DGE into simple square gridded data using sliding window strategy and outputs RDS file with collapsed barcodes and spatial information. All the results are store in the current working directory.
### *input*
* --STtools: Path to STtools package. If not given, the current working directory is used.(add this to pkg).
* --tiles: Tiles that the user is insterested in. Multiple tile numbers need to be separated by comma. For example: --tiles 2106,2107,2108. Required.
* --binsize: The size of the square grids side. By default, it is set to 300 units.
* --window: The side of sliding window. By default it is set to by 150 units.
* --DGEdir: Path to the digital expression matrix. If not given, use the path to DGE from previous steps.
* --spatial: Path to the txt file of spatial coordinates. If not given, use the path to spatialCoordinates.txt from previous steps.
* --nrow: number of rows when generating the super tile. If not give, we assume there is only one tile, and nrow=1.
* --ncol: number of cols when generating the super tile. If not give, we assume there is only one tile, and ncol=1.
* --outdir: Path to output files. If not given, using current working directory
### *code*
```sh
## $STHOME indicates the path to the directory of STtools repository
export STHOME=/path/to/STtools
## $STDATA indicates the directory containing the example input files
export STDATA=/path/to/data
## $STOUT indicates the path to digital expression matrix(DGE)
export STOUT=/path/to/outdir
## $STOUT indicates the path to digital expression matrix(DGE)
export STDGE=/path/to/DGE/
python3 $STHOME/sttools.py --run-steps 5 --STtools $STHOME --spatial $STDATA/spatialcoordinates.txt   --outdir $STOUT --tiles 2106 --binsize 300  -l 20 --window 150 --DGEdir $STDGE

```
### *output*
* SlidingSqureGrids.RDS



## Step 6
This step conducts clustering pipeline based on Seurat tutorial. And mapping the celltype to sliding grids using simple grids as query. For more details about Seurat clustering and maping, please refer to https://satijalab.org/seurat/articles/spatial_vignette.html  and https://satijalab.org/seurat/articles/integration_mapping.html.  Step 6 outputs two RDS files in the current working directory.
### *Input*
* --STtools: Path to STtools package. If not given, the current working directory is used.(add this to pkg).
* --simpleGridsPath: Path to SimpleSquareGrids.RDS from step 4.
* --slidingGridsPath: Path to SlidingSquareGrids.RDS from step 5,
* --geneCount1: Cutoff of nFeatures for SimpleSquareGrids.RDS .
* --geneCount2: Cutoff of nFeatures for SlidingSquareGrids.RDS.
* --nFeaturePlotOnly: If TRUE, output violin plot of nFeatures only. If FALSE, generate mapped RDS as well. 
### *Code*
```sh
## $STHOME indicates the path to the directory of STtools repository
export STHOME=/path/to/STtools
## $STDATA indicates the directory containing the example input files
export STDATA=/path/to/spatial/data
## $STSIMPLE indicates the path to SimpleSquareGrids.RDS
export STSIMPLE=/path/to/SimpleSquareGrids
## $STSIMPLE indicates the path to SlidingSquareGrids.RDS
export STSLIDING=/path/to/SlidingSquareGrids
## $STOUT indicates the path to digital expression matrix(DGE)
export STOUT=/path/to/outdir
python3 $STHOME/sttools.py --run-steps 6 --STtools $STHOME   --outdir $STOUT --simpleGridsPath $STSIMPLE/SimpleSquareGrids.RDS --slidingGridsPath $STSLIDING/SlidingSquareGrids.RDS


```
### *Output*
*  slidingGrid_mapping.RDS: Seurat object with cell type/clustering mapping
*  simpleGrid_clus.RDS: Seurat object with Seurat clusterings

## Step7
This step generate the spliced and unspliced plots when the genes are divided into three subsets. The plots are stored in the current working directory. (add spliced and unspliced without subsetting genes)
### *Input*
* --STtools: Path to STtools package. If not given, the current working directory is used.(add this to pkg).
* --spatial: Path to the txt file of spatial coordinates. If not given, use the path to spatialCoordinates.txt from previous steps.
* --subDGEdir:  Path to the digital expression of spliced and unspliced matrix generated by STARsolo Velocyto .
* --tiles: Tiles that the user is insterested in. Multiple tile numbers need to be separated by comma. For example: --tiles 2106,2107,2108. Required.
* --alpha: Transparency. If not given, alpha=0.01.
* --py: path to python executives

### *Code*
```sh
## $STHOME indicates the path to the directory of STtools repository
export STHOME=/path/to/STtools
## $STDATA indicates the directory containing the example input files
export STDATA=/path/to/spatial/data
## $STDATA indicates the path to unspliced digital expression matrix
export STSUBDGE=/path/to/unspliced/DGE
## $STOUT indicates the path to digital expression matrix(DGE)
export STOUT=/path/to/outdir
python3 $STHOME/sttools.py --run-steps 7 --STtools $STHOME   --outdir $STOUT --spatial $STDATA/spatialcoordinates.txt --subDGEdir $STSUBDGE -alpha 0.02 --tiles 2106,2107


```
### *Output*
* splice_subset*.png
* unsplice_subset*.png


  
  
 
