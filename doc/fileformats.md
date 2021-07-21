
# Input File formarts 
## Step 1
* --first-fq: FASTQ.GZ file. The first N bases are HDMIs. Coordinates information of HDMIs  is stored in the header of FASTQ.GZ file.
<p align="center">
    <img src="./firstseq.png" width="800" height="400" />
</p>
 

* --second-fq1: FASTQ.GZ file. The first N bases are HDMIs, which are reverse complementary of the HDMIs from --first-fq.
<p align="center">
    <img src="./secondseq_fq1.png" width="800" height="400" />
</p>
 
 
## Step 2
* --HDMI2ndSeq: txt file with only one column, representing the HDMIs from --first-fq. See following:
<p align="center">
    <img src="./HDMI2ndSeq.png" width="800" height="400" />
</p>
* --spatial: txt file with five columns: HDMI, lane, tile, X and Y (no header). See the following:
<p align="center">
    <img src="./spatialcoor.png" width="400" height="400" />
</p>

## Step 3
* --second-fq1: FASTQ.GZ file. See Step 1 --second-fq1 input format.
* --second-fq1: FASTQ.GZ file. The first N bases are randomers or UMIs, followed by cDNA sequences.
<p align="center">
    <img src="./secondseq_fq2.png" width="800" height="400" />
</p>

## Step 4
* --spatial: txt file with five columns: HDMI, lane, tile, X and Y. Same as step 2
* --DGEdir: This path should contain 3 files: barcodes.tsv, features.tsv and matrix.mtx. 
## Step 5
* --spatial:txt file with five columns: HDMI, lane, tile, X and Y. Same as step 2 
* --DGEdir: This path should contain 3 files: barcodes.tsv, features.tsv and matrix.mtx, the same as Step 4.
## Step 6
* --

## Step 7
* --spatial: txt file with five columns: HDMI, lane, tile, X and Y.Same as step 2
* --subDGEdir: This path should contain 3 files: barcodes.tsv, features.tsv and matrix.mtx. The last file matrix.mtx should contain unspliced and spliced reads as following:

