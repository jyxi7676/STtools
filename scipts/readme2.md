
# User Manual
## Running user defined consecutive steps in STtools
### Input
  The user need to provide the following files in order for STtools to run all the steps sequentially. 
  *   --seq1: Path to 1st-Seq FASTQ.gz file. STtools takes in fastq.gz files with specific structure(link here). If the barcode/UMI/randomer location is different, please see this link for an example to make it compatible to STtools package. 
  *   --fq1 : Path to 2nd-Seq Read 1 FASTQ.gz file. If the barcode/UMI/randomer location is different, please see this link for an example to make it compatible to STtools package.
  *   --fq2: Path to 2nd-Seq Read 2 FASTQ.gz file. If the barcode/UMI/randomer location is different, please see this link for an example to make it compatible to STtools package.
  *   --star-path: Path to STAR executable
  *   --seqtk-path: Path to seqtk executable
  *   --genome: Path to gene index folder needed for alignement in Step 2. To generate the genome index, one example can be found at https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html
  *   --hdmilength: An integer of the barcode/HDMI length. By default is 20. This version only takes 20 or 32(need to modify this)
  *   --outprefix: Output prefix for STARsolo alignemnt, if not given, outprefix will be set to 'Sample'
  *   --STtools: Path to the STtools package (modify it with a default path)
  *   --py: Binary path to python executable
  *   --tiles: The tiles that the user are interested in for Step 3,4,5,6. If more than one tiles are given, separate them by comma.
  *   --sidesize: Side length of square grid. By default, sidesize=300
  *   --window: Sise of sliding window. By default, window=150
  *   --cores: Number of cores for parralell computing of step 5. By default is 5.
  *   --maxScale: max color bar value for HDMI discovery plot in step 2, if not given, will use the max number of HDMI reads.
  *   --alpha: Transparency for plotting in step 6. If not given, alpha=0.01
### Output: 
  There are intermediate outputs for each step. Please refer to https://github.com/jyxi7676/STtools/blob/main/scipts/readme3.md for details of the outputs in each step. Among all the inputs, the most usefull files are:
  * Digital expression matrix (from step 3)
  * Seurat object xx.RDS(from step 4 and 5)
### Example Data:
  
### Example Code:
  * ```
  python3 sttools_v3.py --run-all 3 --seq1 '/net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/extractCoord/input/liver-MiSeq-tile2106-sub-R1.fastq.gz' --fq1 '/net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/align/input/liver_tile2106_sub_R1.fastq.gz' --fq2 '/net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/align/input/liver_tile2106_sub_R2.fastq.gz' -g '/net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/geneIndex/' --STtools '/net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/' --star-path '/net/fantasia/home/jyxi/STAR-2.7.5c/source/' --seqtk-path '/net/fantasia/home/jyxi/seqtk/' --py 'python3' --tiles 2106 --sidesize 300 --window 150 -l 20 -o 'Sample' -c 2
  ```
