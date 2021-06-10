
# Input File formarts 
## Step 1
* --fq1: FASTQ.GZ file. Sequence structure is as follows: 
* --seq1: FASTQ.GZ file. seq1 has the spatial information and HDMI information as the following:
 
## Step 2
* --HDMI2ndSeq: txt file with only one column, representing the HDMIs from --seq1. See following:
* --spatial: txt file with five columns: HDMI, lane, tile, X and Y. See the following:

## Step 3
* --fq1: FASTQ.GZ file. Sequence structure is as follows: 
* --fq2: FASTQ.GZ file. Sequence structure is as follows: 
## Step 4
* --spatial: txt file with five columns: HDMI, lane, tile, X and Y. Same as step 2
* --DGEdir: This path should contain 3 files: barcodes.tsv, features.tsv and matrix.mtx. See the following:
## Step 5
* --spatial:txt file with five columns: HDMI, lane, tile, X and Y. Same as step 2 
* --DGEdir: This path should contain 3 files: barcodes.tsv, features.tsv and matrix.mtx, the same as Step 4.
## Step 6
* --

## Step 7
* --spatial: txt file with five columns: HDMI, lane, tile, X and Y.Same as step 2
* --subDGEdir: This path should contain 3 files: barcodes.tsv, features.tsv and matrix.mtx. The last file matrix.mtx should contain unspliced and spliced reads as following:

