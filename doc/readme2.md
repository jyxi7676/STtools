
# User Manual
## Running user defined consecutive steps in STtools
STtools can run user specified consecutive steps with the option: --run-steps
For example: --run-steps A1,A2,A3 or --run-steps C1,C2. The steps must be separated by comma and must be consecutive. The command will raise an error if the steps are not continuous  such as --run-steps A3,C2,C3.

### Input&Output
 The inputs and outputs for the command depend on the optional steps. Please refer to [running specific step](./doc/readme3.md)
 for more details of each step. Here, we will give you some examples for better illustration.
 

### Example 1:
If the user wants to output the DGE from raw fastq.gz files, then the user can specify --run-steps 1,2,3 as follows:
```
## $STHOME indicates the path to the directory of STtools repository
export STHOME=/path/to/STtools
## $STDATA indicates the directory containing the 1st-Seq FASTQ.gz and 2nd-Seq Read 1 FASTQ.gz file
export STDATA=/path/to/data
## $STOUT indicates the directory containing the output files.
export STOUT=/path/to/outdir
## $GENOMEINDEX indicates the path to genome index for STARsolo alignment
export GENENOME=/path/to/geneIndex
## $SEQTKPATH indicates the path to the seqtk executive
export SEQTKPATH=/path/to/seqtk/executive
## $STARPATH indicates the path to the STAR executive
export STARPATH=/path/to/star/executive

python3 $STHOME/sttools.py --run-steps A1,A2,A3 --first-fq $STDATA/liver-MiSeq-tile2106-sub-R1.fastq.gz --second-fq1 $STDATA/liver_tile2106_sub_R1.fastq.gz --second-fq2 $STDATA/liver_tile2106_sub_R2.fastq.gz  --STtools $STHOME  -l 20  --genome $GENENOME --star-path $STARPATH --seqtk-path $SEQTKPATH  --outdir $STOUT
```
### Example 2
If the user wants to start from DGE and spatial coordinates and generate simple square grids, sliding square grids and reference mapping. Then the user can specify --run-steps C1,C2,C3 as follows:


```
## $STHOME indicates the path to the directory of STtools repository
export STHOME=/path/to/STtools
## $STDATA indicates the directory containing the spatialcoordinats.txt
export STDATA=/path/to/spatial/data
## $STDGE indicates the path to the digital expression matrix 
export STDGE=/path/to/DGE
## $STOUT indicates the directory containing the output files.
export STOUT=/path/to/outdir

python3 $STHOME/sttools.py --run-steps C1,C2,C3  --STtools $STHOME   --spatial $STDATA/spatialcoordinates.txt --tiles 2106 -l 20 --window 150 --binsize 300  -l 20 --DGEdir $STDGE --outdir $STOUT
```



 
