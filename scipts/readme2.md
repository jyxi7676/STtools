
# User Manual
## Running user defined consecutive steps in STtools
STtools can run user specified consecutive steps with the option: --run-steps
For example: --run-steps 1,2,3 or --run-steps 4,5. The steps must be separated by comma and must be consecutive. The command will raise an error if the steps are not continuous integers such as --run-steps 3,5,6.

### Input&Output
  The inputs and outputs for the command depend on the optional steps. Please refer to https://github.com/jyxi7676/STtools/blob/main/scipts/readme3.md for more details about the details for each step. Here, we will give you some examples for better illustration.
 

### Example 1:

```
 python3 sttools_v4.py --run-steps 2,3,4,5 --fq1 '/net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/align/input/liver_tile2106_sub_R1.fastq.gz' --fq2 '/net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/align/input/liver_tile2106_sub_R2.fastq.gz' -g '/net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/geneIndex/' --STtools '/net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/' --star-path '/net/fantasia/home/jyxi/STAR-2.7.5c/source/' --seqtk-path '/net/fantasia/home/jyxi/seqtk/' --py 'python3' --tiles 2106 --sidesize 300 --window 150 -l 20 -o 'Sample' -c 2 --whitelist '/net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/whitelist.txt' --spatial '/net/fantasia/home/jyxi/scrna/leejun/ngst/STtools/spatialcoordinates.txt'
```
### Example 2

### Example 3
 
