
# Input File formarts 
## Step 1
* --first-fq: FASTQ.GZ file. The first N bases are HDMIs. Coordinates information (lane, tile, X and Y) of HDMIs  is stored in the header of FASTQ.GZ file. Please click the [FASTQ file format](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm) for more details. An example of our sample data is shown below:
<p align="center">
    <img src="./firstseq.png" width="800" height="250" />
</p>
 

* --second-fq1: FASTQ.GZ file. The first N bases are HDMIs. For example, in the example data 1, we have N=20.
<p align="center">
    <img src="./secondseq_fq1.png" width="800" height="250" />
</p>

--second-fq2: FASTQ.GZ file. Read 2 are cDNA sequences and the first 9 bases are randomers. 
<p align="center">
    <img src="./secondseq_fq2.png" width="800" height="250" />
</p> 
 
## Step 2
* --HDMI2ndSeq: txt file with only one column, representing the HDMIs from --first-fq. See following:
<p align="center">
    <img src="./HDMI2ndSeq.png" width="200" height="400" />
</p>
* --spatial: txt file with five columns: HDMI, lane, tile, X and Y (no header). See the following:
<p align="center">
    <img src="./spatialcoor.png" width="400" height="400" />
</p>

## Step 3
* --second-fq1: FASTQ.GZ file. See Step 1 --second-fq1 input format.
* --second-fq2: FASTQ.GZ file. The first N bases are randomers or UMIs, followed by cDNA sequences.

## Step 4
* --spatial: txt file with five columns: HDMI, lane, tile, X and Y. Same as step 2
* --DGEdir: This path should contain 3 files: barcodes.tsv, features.tsv and matrix.mtx. They are typical format for single cell RNA-seq. The format can be seen at the [link](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices)
* --layout: Path to a .csv file with the arrangement information of super tiles. 
<p align="center">
    <img src="./layout.png" width="400" height="200" />
</p>

## Step 5
* --spatial:txt file with five columns: HDMI, lane, tile, X and Y. Same as step 2 
* --DGEdir: This path should contain 3 files: barcodes.tsv, features.tsv and matrix.mtx, the same as Step 4.
## Step 6
* --simpleGridsPath: Path to SimpleSquareGrids.RDS. The RDS file will contain the "image" slots with the spatial coordinates and meta data along with the collapsed count data from Step 4.
* --slidingGridsPath: Path to slidingSquareGrids.RDS. The RDS file will contain the "image" slots with the spatial coordinates and meta data along with the collapsed count data from Step 5.

## Step 7
* --spatial: txt file with five columns: HDMI, lane, tile, X and Y.Same as step 2
* --subDGEdir: This path should contain 3 files: barcodes.tsv, features.tsv and matrix.mtx. The first two files are typical file format from scRNA-seq. As for the matrix.mtx, it contain  spliced, unspliced, and ambiguous reads (col3=spliced, col4=unspliced, col5=ambiguous) as the following. 
<p align="center">
    <img src="./unspliced.png" width="800" height="400" />
</p>


