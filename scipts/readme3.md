
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
Step2 visualize the barcode/HDMI density discovery plot, with which the user are able to compare with HE images to estimate the tissue boundary. The alignment is done manually and automatic alignment is under development. 
### *Input*
  * -fq1
  * -fq2
### *Code*
 ```
 python3 /net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/sttools_v4.py --run-steps 3 --STtools '/net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/'  --py 'python3' --spatial /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129-15/spatialcoordinates.txt --hdmi2ndSeq /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129-15/HDMI_SeqScope_2nd.txt

 
 python3 /net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/sttools_v4.py --run-steps 2 --fq1 '/net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129-15/19129FL-15-01-02_S87_L008_R1_001.fastq.gz' --fq2 '/net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129-15/19129FL-15-01-02_S87_L008_R2_001.fastq.gz'  -g '/net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/geneIndex/' --STtools '/net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/' --star-path '/net/fantasia/home/jyxi/STAR-2.7.5c/source/' --seqtk-path '/net/fantasia/home/jyxi/seqtk/' --py 'python3'  --py 'python3' -o 'ColonWTA'

 ```
### *Output*

## Step 3
Step 3 aligns the data the reference genome using STARsolo software and output digital expression matrix under Gene,GeneFull, and Velocyto options (see link xxxxx).
Before running step3, it is the required the user to generate the genome reference following this link  xxxxx.
### *Input*
 * --fq1: Path to 2nd-Seq FASTQ.gz file of read 1. Required.
 * --fq2: Path to 2nd-Seq FASTQ.gz file of read 2. Required.
 * -g: genome reference. Required
 * --STtools: Path to the STtools package. If not given, using currently working directory (add this to pkg)
 * -o: Output prefix of alignment, if not given, set to 'Sample'
 * --sesqtk-path: Path to seqtk executable. Required
 * --star-path: Path to STAR executable. Required
 * --whitelist: Txt file of the whitelist, if not given, searching whitelist.txt in the current folder.

 
### *Code*
``` 
python3 /net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/sttools_v4.py --run-steps 2 --fq1 '/net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129-15/19129FL-15-01-02_S87_L008_R1_001.fastq.gz' --fq2 '/net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129-15/19129FL-15-01-02_S87_L008_R2_001.fastq.gz'  -g '/net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/geneIndex/' --STtools '/net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/' --star-path '/net/fantasia/home/jyxi/STAR-2.7.5c/source/' --seqtk-path '/net/fantasia/home/jyxi/seqtk/' --py 'python3'  --py 'python3' -o 'ColonWTA'
```
### *Output*
This step outputs folders with STARsolo summary statistics, bam file, DGE, etc.
* SampleSolo.out/
  
## Step 4
Step 4 bins the DGE into simple square gridded data and collapses the reads counts within each grid. A new DGE is created and spatial information is updated with the center of each bin, and is saved as an RDS file for output.
### *input*
### *code*
### *output*
* SimpleSquareGrids.RDS

## Step 5
## Step 6
  
  
 
