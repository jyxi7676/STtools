
# User Manual
## Running specific step in STtools
The user need to use the option --run-step if interested in running one step only per command . The required input/output data format can be found at https://github.com/jyxi7676/STtools/blob/main/scipts/fileformats.md.  We will give a detailed illustration about how to run STtools per step with several examples. 
## Step 1
Step 1 aims to extracts spatial coordinates and whitelist info from the sequenced raw FASTQ.gz file. It assumes the user has the SeqScope data format (see xxxx for stucture, cite seqscope paper), where spatial  realated information such as HDMI/Barcode, lane, tile, X and Y coordinates can be retrieved from 1st-Seq and transcriptomic information can be retrieved from 2nd-Seq. 

### Input
  The user need to provide some of the following files in order for STtools to run all the steps sequentially. 
  *   --seq1: Path to 1st-Seq FASTQ.gz file. STtools takes in fastq.gz files with specific structure(link here). If the barcode/UMI/randomer location is different, please see this link for an example to make it compatible to STtools package. 
  *   --fq1 : Path to 2nd-Seq Read 1 FASTQ.gz file. If the barcode/UMI/randomer location is different, please see this link for an example to make it compatible to STtools package.
  *   --fq2: Path to 2nd-Seq Read 2 FASTQ.gz file. If the barcode/UMI/randomer location is different, please see this link for an example to make it compatible to STtools package.
## Step 2
## Step 3
## Step 4
## Step 5
## Step 6
  
  
 
