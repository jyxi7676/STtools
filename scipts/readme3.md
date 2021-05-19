
# User Manual
## Running specific step in STtools
The user need to use the option --run-step if interested in running one step only per command . The required input/output data format for each step can be found at https://github.com/jyxi7676/STtools/blob/main/scipts/fileformats.md.  We will give a detailed illustration about how to run STtools per step with several examples. 
## Step 1
Step 1 aims to extracts spatial coordinates and whitelist info from the sequenced raw FASTQ.gz file. It assumes the user has the SeqScope data format (see xxxx for stucture, cite seqscope paper), where spatial  realated information such as HDMI/Barcode, lane, tile, X and Y coordinates can be retrieved from 1st-Seq and transcriptomic information can be retrieved from 2nd-Seq. 
### *Input*
  The user need to provide some of the following files in order for STtools to run Step 1. 
  *   --seq1: Path to 1st-Seq FASTQ.gz file. STtools takes in fastq.gz files with SeqScope sequence design structure. If the barcode/UMI/randomer location is different, please see this link xxxxx for an example to make the inputs compatible to STtools package. **Required**. 
  *   --fq1 : Path to 2nd-Seq Read 1 FASTQ.gz file. If the barcode/UMI/randomer location is different, please see this link for an example to make it compatible to STtools package.**Required**. 
  *   --hdmilength: An integer of the length of HDMI/Barcode.(modify to take any number below 32?). By default hdmilength=20.
  *   --STtools: Path to the STtools package. If not given, using currently working directory (add this to pkg)
 ### *Code*
 ```
 python3 /net/fantasia/home/jyxi/scrna/leejun/ngst/STtools//sttools_v4.py --run-steps 1 --seq1 '/net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/RDuo3/RDuo3_S1_L001_R1_001.fastq.gz' --fq1 '/net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129-15/19129FL-15-01-02_S87_L008_R1_001.fastq.gz'  --STtools '/net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/'  

 ```
 ### *Output*
 This step output two useful files in the current working directory, and are taken as input the next steps.
 * spatialcoordinates.txt: 
 * whitelist.txt:
 
## Step 2
## Step 3
## Step 4
## Step 5
## Step 6
  
  
 
