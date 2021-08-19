
# User Manual
## Running all steps in STtools
Before running all steps, it is suggested that the user take a look at the  [running specific step](./doc/readme3.md).
The user need to use the option --run-all to run all steps and to provide the following input files. 
### Input
  The user need to provide some of the following files in order for STtools to run all the steps sequentially. 
  *   --STtools: Path to the STtools package. If not given, using current working directory.
  *   --outdir: Path to output files. If not given, using current working directory
  *   --first--fq: Path to 1st-Seq FASTQ.gz file. **Required**
  *   --second-fq1 : Path to 2nd-Seq Read 1 FASTQ.gz file. **Required**
  *   --second-fq2: Path to 2nd-Seq Read 2 FASTQ.gz file. If the barcode/UMI/randomer location is different, please see this link for an example to make it compatible to STtools package.** Required **
  *   --star-path: Path to STAR executable. **Required**
  *   --seqtk-path: Path to seqtk executable. **Required**
  *   --genome or -g: Path to genome index needed for alignement in Step 2. **Required**
  *   --hdmilength or -l: An integer of the barcode/HDMI length. By default is 20. 
  *   -whitelist:txt file of the whitelist
  *   --outprefix or -o: Output prefix for STARsolo alignemnt, if not given, outprefix will be set to 'Sample'
  *   --tiles: The tiles that the user are interested in for Step 3,4,5,6. If more than one tiles are given, separate them by comma.
  *   --binsize: Side length of square grid. By default, sidesize=300
  *   --window: Sise of sliding window. By default, window=150
  *   --cores: Number of cores for parralell computing of step 5. By default is 5.
  *   --maxScale: max color bar value for HDMI discovery plot in step 2, if not given, will use the max number of HDMI reads.
  *   --alpha: Transparency for plotting in step 6. If not given, alpha=0.01
### Output: 
  There are intermediate outputs for each step. Please refer to https://github.com/jyxi7676/STtools/blob/main/scipts/readme3.md for details of the outputs in each step. Among all the inputs, the most usefull files are:
  * Digital expression matrix (from step 3)
  * RDS file of Seurat object(from step 4 and 5)
  
### Example Code:
   ```
  ## $STHOME indicates the path to the directory of STtools repository
  export STHOME=/path/to/STtools
  ## $STDATA indicates the directory containing the example input files:spatialcoordinates.txt
  export STDATA=/path/to/data
  ## $STOUT indicates the directory containing the example output files
  export STOUT=/path/to/outdir
  ## $GENEINDEX indicates the path to genome index for STARsolo alignment
  export GENOMEINDEX=/path/to/genomeIndex
  ## $SEQTKPATH indicates the path to the seqtk executive
  export SEQTKPATH=/path/to/seqtk/executive
  ## $STARPATH indicates the path to the STAR executive
  export STARPATH=/path/to/star/executive
  
  python3 $STHOME/sttools.py --run-all --STtools $STHOME  --first-fq $STDATA/liver-MiSeq-tile2106-sub-R1.fastq.gz --second-fq1 $STDATA/liver_tile2106_sub_R1.fastq.gz --second-fq2 $STDATA/liver_tile2106_sub_R2.fastq.gz --outdir $STOUT --genome $GENOMEINDEX --star-path $STARPATH --seqtk-path $SEQTKPATH --tiles 2106 --binsize 300 --window 150 -l 20 -o 'Sample' -c 2
  ```
