
# User Manual
## Running specific step in STtools
The user needs to use the option --run-steps if interested in running steps separately. (We will give a detailed illustration about how to run STtools per step with several examples,link). Please modify your input files format according to  [the link](./doc/fileformats.md) 

## Step 1
Step 1 aims to extracts spatial coordinates, whitelist and HDMIs from the sequenced raw FASTQ.gz file. It assumes the user has the SeqScope data format (see [Seq-Scope](https://www.cell.com/cell/fulltext/S0092-8674(21)00627-9?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421006279%3Fshowall%3Dtrue) paper), where spatial information such as HDMI/Barcode, lane, tile, X and Y coordinates can be retrieved from 1st-Seq and transcriptomic information can be retrieved from 2nd-Seq. If the data are from VISIUM or SlideSeq, we need to add this to our tool?

### *Input*
  The user needs to provide some of the following files to run Step 1. 
  *   --first--fq: Path to 1st-Seq FASTQ.gz file. STtools takes in fastq.gz files with SeqScope sequence design structure. If the barcode/UMI/randomer location is different, please see this link xxxxx (we need to add this ???) for an example to make the inputs compatible to STtools package. **Required**. 
  *   --second-fq1 : Path to 2nd-Seq Read 1 FASTQ.gz file. If the barcode/UMI/randomer location is different, please see this link for an example to make it compatible to STtools package.**Required**. 
  *   --hdmilength or -l: An integer of the length of HDMI/Barcode.If the length of HDMIs will be set to 30 if it is bigger than 30. By default hdmilength=20.
  *   --STtools: Path to the STtools package. If not given, using current working directory. **Required**
  *   --outdir: Path to output files. If not given, using current working directory
 ### *Code*
 ```sh
 ## $STHOME indicates the path to the directory of STtools repository
 export STHOME=/path/to/STtools 
 ## $STDATA indicates the directory containing the 1st-Seq FASTQ.gz and 2nd-Seq Read 1 FASTQ.gz file
 export STDATA=/path/to/data
 ## $STOUT indicates the directory containing the output files.
 export STOUT=/path/to/outdir
 
 python3 $STHOME/sttools.py --run-steps 1 --first-fq $STDATA/liver-MiSeq-tile2106-sub-R1.fastq.gz --second-fq1 $STDATA/liver-HiSeq-tile2106-sub-R1.fastq.gz --STtools $STHOME  -l 20 --outdir $STOUT
 ```
 ### *Output*
 This step outputs the following files in the current working directory, and are taken as input for next steps.
 * spatialcoordinates.txt 
 * whitelist.txt
 * HDMI_SeqScope_2nd.txt
 
## Step 2
Step2 visualizes the barcode/HDMI density discovery plot, with which the user is able to compare with histology images to estimate the tissue boundary. The alignment is done manually and automatic alignment is under development. 
### *Input*
  * --STtools: Path to the STtools package. If not given, using current working directory.
  * --hdmi2ndSeq: txt file with the HDMIs from the --second-fq1. For exmaple, HDMI_SeqScope_2nd.txt from step 1. **Required**
  * --spatial: txt file with HDMI and spatial information from --first-fq: HDMI, lane, tile, X and Y. For example, spatialcoordinates.txt from step 1 **Required**
  * --maxScale: max scale value for color bar, if not specified, using default value in matplotlib.
  * --outdir: Path to output files. If not given, using current working directory

### *Code*
 ```sh
  ## $STHOME indicates the path to the directory of STtools repository
  export STHOME=/path/to/STtools
  ## $STDATA indicates the directory containing the example input files: --spatial and --hdmi2ndSeq
  export STDATA=/path/to/data
  ## $STOUT indicates the directory containing the example output files.
  export STOUT=/path/to/outdir
  
  python3 $STHOME/sttools.py --run-steps 2 --STtools $STHOME --spatial $STDATA/spatialcoordinates.txt --hdmi2ndSeq $STDATA/HDMI_SeqScope_2nd.txt --outdir $STOUT
  ```
### *Output*
lane_tile*.png

## Step 3
 Step 3 maps the reads to the reference genome using the standard STAR  read alignment algorithm and outputs digital expression matrix under Gene,GeneFull, and Velocyto options (change here? add option?).
 To run Step 3, the user needs to genete gene index as an input for STARsolo alignment.
 
### *Input*
 * --second-fq1: Path to 2nd-Seq FASTQ.gz file of read 1. Required.
 * --second-fq2: Path to 2nd-Seq FASTQ.gz file of read 2. Required.
 * --hdmilength or -l: An integer of the length of HDMI/Barcode.If the length of HDMIs will be set to 30 if it is bigger than 30. By default hdmilength=20.
 * --whitelist: Txt file of the whitelist, if not given, searching whitelist.txt in the current folder.
 * --genome or -g: STAR geneome index. **Required**
 * --STtools: Path to the STtools package. If not given, using currently working directory.
 * -o: Output prefix of alignment, if not given, set to 'Sample'
 * --sesqtk-path: Path to seqtk executable. **Required**
 * --star-path: Path to STAR executable. **Required**
 * --outdir: Path to output files. If not given, using current working directory

 
### *Code*
``` sh
 ## $STHOME indicates the path to the directory of STtools repository
 export STHOME=/path/to/STtools
 ## $STDATA indicates the directory containing 2nd-Seq R1 and R2 files in FASTQ.GZ format
 export STDATA=/path/to/data
 ## $STWHITELIST indicates the directory containing the whiteslist.txt. If not given,searching in the current directory
 export STWHITELIST=/path/to/whitelist
 ## $STOUT indicates the directory containing the example output files.
 export STOUT=/path/to/outdir
 ## $GENOMEINDEX indicates the path to genome index for STARsolo alignment
 export GENOMEINDEX=/path/to/genomeIndex
 ## $SEQTKPATH indicates the path to the seqtk executive
 export SEQTKPATH=/path/to/seqtk/executive
 ## $STARPATH indicates the path to the STAR executive
 export STARPATH=/path/to/star/executive
 
python3 $STHOME/sttools.py --run-steps 3 --STtools $STHOME --whitelist $STWHITELIST/whitelist.txt --second-fq1 $STDATA/liver_tile2106_sub_R1.fastq.gz --second-fq2 $STDATA/liver_tile2106_sub_R2.fastq.gz --outdir $STOUT --genome $GENOMEINDEX --star-path $STARPATH --seqtk-path $SEQTKPATH

```
### *Output*
This step outputs folders with STARsolo summary statistics, bam file, DGE, etc. The key outputs are as follows:
* SampleSolo.out/
* SampleAligned.sortedByCoord.out.bam
  
## Step 4
Step 4 conducts different functions depending on the input --datasource.

If datasource is 'SeqScope', then Step 4 bins the DGE into simple square gridded data and collapses the reads counts within each grid. A new DGE is generated and spatial information is updated with the center of each bin. This step outputs a RDS file with collapsed barcodes and spatial information. 
If datasource is 'VISIUM', then Step 4 will implement either Seurat pipeline or BayesSpace pipeline depending on the --algo option.
If datasource is 'SlideSeq', then Step 4 will implement either simple square gridding (by default) or Seurat pipeline depending on --algo option.

### *input*
* --STtools: Path to STtools package. If not given, the current working directory is used.
* --binsize: The size of the square grids side. By default, it is set to 300 pixel units, which is 10um in SeqScope data.
* --DGEdir: Path to the digital expression matrix. **Required**
* --spatial: Path to the txt file of spatial coordinates. **Required**
* --lane-tiles: Lane and tiles that the users are interested,for exaample: 1_2106 are lane 1 and tile 2106. By default, all lanes and tiles will be used. Otherwise, the user should give a value for the parameter. Multiple lane and tiles should be separated by comma. For instance: --lane-tiles 1_2106, 2_2106, 1_2107. 
* --seqscope1st: String indicating the layout template, either 'MiSeq' or 'HiSeq' or 'Custom'. If the seqscope1st is 'MiSeq' or 'HiSeq', the supertiles will be generated according to the template layout(link) among --lane-tiels. If seqscope1st is 'Custom' then the user should give the --layout with user specified layout if more than 1 tiles are used. 
* --layout: User can have a customized arrangment of lanes and tiles, and --layout is the path to the file. Please click [layout](./doc/Layout_template.png) for the input format. 
* --outdir: Path to output files. If not given, using current working directory.
### *code*
```sh
## $STHOME indicates the path to the directory of STtools repository
export STHOME=/path/to/STtools
## $STDATA indicates the directory containing the example input files:spatialcoordinates.txt
export STDATA=/path/to/data
## $STOUT indicates the directory containing the example output files
export STOUT=/path/to/outdir
## $STDGE indicates the path to digital expression matrix(DGE)
export STDGE=/path/to/DGE/

python3  $STHOME/sttools.py --run-steps 4 --STtools $STHOME --spatial $STDATA/spatialcoordinates.txt   --outdir $STOUT --binsize 300 --DGEdir $STDGE --layout HiSeq 
```
### *output*
* SimpleSquareGrids.RDS
* collapsedBarcodes.csv			
* collapsedGenes.csv			 
* collapsedMatrix.mtx

## Step 5
Step 5 bins the DGE into simple square gridded data using sliding window strategy and outputs a RDS file with collapsed barcodes and spatial information. 
### *input*
* --STtools: Path to STtools package. If not given, the current working directory is used.(add this to pkg).
* --lane-tiles: Tiles that the user is insterested in. Multiple tile numbers need to be separated by comma. For example: --lane-tiles 1_2106,1_2107,1_2108. **Required**.
* --binsize: The size of the square grids side. By default, it is set to 300 units.
* --window: The side of sliding window. By default it is set to by 150 units.
* --DGEdir: Path to the digital expression matrix. If not given, use the path to DGE from previous steps.
* --spatial: Path to the txt file of spatial coordinates. If not given, use the path to spatialCoordinates.txt from previous steps.
* --nrow: number of rows when generating the super tile. If not give, we assume there is only one tile, and nrow=1.
* --ncol: number of cols when generating the super tile. If not give, we assume there is only one tile, and ncol=1.
* --outdir: Path to output files. If not given, using current working directory
*  --layout: User can have a customized arrangment of lanes and tiles, and --layout is the path to the file. Please click [layout](./Layout_template.png) for the input format.
* --order: How the lane and tiles are stacked when --layout is not given. If order is 'top' then the --lane-tiles will be ordered in increasing order and stacked in 2 rows. If order is 'bottom', they will be ordered in a decreasing order.  
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
python3 $STHOME/sttools.py --run-steps 5 --STtools $STHOME --spatial $STDATA/spatialcoordinates.txt   --outdir $STOUT --lane-tiles 1_2106 --binsize 300  -l 20 --window 150 --DGEdir $STDGE

```
### *output*
* SlidingSqureGrids.RDS



## Step 6
Step 6 conducts clustering pipeline based on Seurat tutorial. And mapping the celltype to sliding grids using simple grids as query. For more details about Seurat clustering and maping, please refer to https://satijalab.org/seurat/articles/spatial_vignette.html  and https://satijalab.org/seurat/articles/integration_mapping.html.  This step requires the user to annotate the clustering results from simpleSquareGrids.RDS for a better result. Need to make it compatible to take in single cell dataset as query dataset.
### *Input*
* --STtools: Path to STtools package.  **Required**.
* --simpleGridsPath: Path to SimpleSquareGrids.RDS, **Required**.
* --slidingGridsPath: Path to SlidingSquareGrids.RDS,**Required**.
* --geneCount1: Cutoff of nFeatures for SimpleSquareGrids.RDS. If not given, geneCount1=0.
* --geneCount2: Cutoff of nFeatures for SlidingSquareGrids.RDS. If not givn geneCount1=0.
* --nFeaturePlotOnly: If TRUE, output violin plot of nFeatures only. If FALSE, generate mapped RDS as well. 
### *Code*
```sh
## $STHOME indicates the path to the directory of STtools repository
export STHOME=/path/to/STtools
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
This step generates the spliced and unspliced plots of the genes that are  divided into three subsets. To generate the plot as in Seq-Scope papers, do we need to modify the plots?
### *Input*
* --STtools: Path to STtools package.  **Required**.
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
python3 $STHOME/sttools.py --run-steps 7 --STtools $STHOME   --outdir $STOUT --spatial $STDATA/spatialcoordinates.txt --subDGEdir $STSUBDGE --alpha 0.02 --tiles 2106,2107


```
### *Output*
* splice_subset*.png
* unsplice_subset*.png


  
  
 
