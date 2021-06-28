
# User Manual
## Running user defined consecutive steps in STtools
STtools can run user specified consecutive steps with the option: --run-steps
For example: --run-steps 1,2,3 or --run-steps 4,5. The steps must be separated by comma and must be consecutive. The command will raise an error if the steps are not continuous integers such as --run-steps 3,5,6.

### Input&Output
  The inputs and outputs for the command depend on the optional steps. Please refer to https://github.com/jyxi7676/STtools/blob/main/scipts/readme3.md for more details about the details for each step. Here, we will give you some examples for better illustration.
 

### Example 1:
If the user wants to output the DGE from raw fastq.gz files, then the user can specify --run-steps 1,2,3 as follows:
```
export STHOME=/path/to/STtools
export STDATA=/path/to/data
export STOUT=/path/to/outdir
export GENEINDEX=/path/to/geneIndex
export SEQTKPATH=/path/to/seqtk/executive
export STARPATH=/path/to/star/executive
python3 $STHOME/sttools_v6.py --run-steps 1,2,3 --first-fq $STDATA/liver-MiSeq-tile2106-sub-R1.fastq.gz --second-fq1 $STDATA/liver_tile2106_sub_R1.fastq.gz --second-fq2 $STDATA/liver_tile2106_sub_R2.fastq.gz  --STtools $STHOME  -l 20  --genome $GENEINDEX --star-path $STARPATH --seqtk-path $SEQTKPATH  --whitelist /net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/temp/whitelist.txt --outdir $STOUT
```
### Example 2
If the user wants to start from DGE and spatial coordinates and generate simple square grids, sliding square grids and reference mapping. Then the user can specify --run-steps 4,5,6 as follows:

```
export STHOME=/path/to/STtools
export STDATA=/path/to/spatial/data
export STDGE=/path/to/DGE
export STOUT=/path/to/outdir
python3 $STHOME/sttools_v6.py --run-steps 4,5,6  --STtools $STHOME   --spatial $STDATA/spatialcoordinates.txt --tiles 2106 -l 20 --window 150 --sidesize 300  -l 20 --DGEdir $STDGE --outdir $STOUT
```



 
