#make all-steps running sequentially!!!!!!

import argparse
import os, sys
import subprocess as sp

parser = argparse.ArgumentParser(description="STtools Main Launcher")
parser.add_argument("--run-all", default=False, action='store_true', help="Run all steps sequentially")
parser.add_argument("--run-steps", type=str, help="Comma-separated string containing the steps to run")
parser.add_argument("--seq1",type=str, help="Path to 1nd-Seq FASTQ.gz file")
parser.add_argument("-l","--hdmilength",type=int,default=20, help="An integer indicating the length of the HDMIs") #if not specify -l erro None occurs
parser.add_argument("-s", "--star-path", type=str, default="STAR", help="Binary folder path to STAR executable") #if not given use current directory
parser.add_argument("-q", "--seqtk-path", type=str, default="seqtk", help="Binary folder path to seqtk executable")
parser.add_argument("-g", "--genome", type=str, help="Path to the reference files")
parser.add_argument("--fq1", type=str, help="Path to 2nd-Seq Read 1 FASTQ.gz file")
parser.add_argument("--fq2", type=str, help="Path to 2nd-Seq Read 2 FASTQ.gz file")
parser.add_argument("-o","--outprefix" ,type=str, help="Prefix for STARsolo output",default='Sample')
parser.add_argument("--STtools",type=str,help='Path to the scripts of STtools scripts',default='STtools')
parser.add_argument("--py",type=str,help='Binary path to python executable')
parser.add_argument("--tiles",type=str,help='tiles')
parser.add_argument("--sidesize",type=int,help='tiles',default=300)
parser.add_argument("--window",type=int,help='tiles',default=150)

##probabily add python binaries?

#parser.add_argument("--cmd", type=str, help="Command to run in step 2")
#parser.add_argument("--out", required=True, type=str, help="Output prefix")

args = parser.parse_args()
steps = []

if ( args.run_all ):
    steps = [1,2,3,4,5]
    print("Running the following steps: ", steps)

    ## check whether parameter is given
    if ( args.seq1 is None ):
        raise ValueError("Cannot find --seq1 argument. Need to specify for running all")
    if ( args.fq1 is None ):
        raise ValueError("Cannot find --fq1 argument. Need to specify for running all")
    if ( args.fq2 is None ):
        raise ValueError("Cannot find --fq2 argument. Need to specify for running all")
    if ( args.genome is None ):
        raise ValueError("Cannot find --genome argument. Need to specify for running all")

    #check file/directory path exist:
    if(os.path.isfile(args.seq1)==False):
        raise ValueError("File --seq1 does not exist")
    if(os.path.isfile(args.fq1)==False):
        raise ValueError("File --fq1 does not exist")
    if(os.path.isfile(args.fq2)==False):
        raise ValueError("File --fq2 does not exist")
    if(os.path.isdir(args.star_path)==False):
        raise ValueError("Directory --star-path does not exist")
    if(os.path.isdir(args.seqtk_path)==False):
        raise ValueError("Directory --seqtk-path does not exist")
    if(os.path.isdir(args.genome)==False):
        raise ValueError("Directory --genome does not exist")
    if(os.path.isdir(args.STtools)==False):
        raise ValueError("Directory --STtools does not exist")



    ###############STEP1##################
   # modify: if files already exist, start from the next steps.....
    # print('Start running Step 1: extract coordinates')
    # cmd1 = "sbatch extractCoord/extractCoord_v2.sh -m {args.seq1} -h {args.fq1} -l {args.hdmilength}".format(args=args)
    # ret = os.system(cmd1)
    # if ( ret != 0 ):
    #     raise ValueError(f"ERROR in running {cmd1}, returning exit code {ret}") 

    ###############STEP2a##################
    # print('Start running Step 2a: STARsolo alignment')
    # whitelist=os.getcwd()+"/whitelist.txt"
    # spatial=os.getcwd()+"/spatialcoordinates.txt"
    # print(spatial)
    # if(os.path.isfile(whitelist)==False):
    #     raise ValueError("Whitelist from STEP 1 output does not exist")
    # if(os.path.isfile(spatial)==False):
    #     raise ValueError("Spatial coordinates from STEP 1 output does not exist")

    # cmd2a = "sbatch align/align_v2.sh -a {args.fq1} -b {args.fq2} -w {args.whitelist} -o {args.outprefix} -t {args.star_path} -q {args.seqtk_path} -g {args.genome}".format(args=args) 
    # ret = os.system(cmd2a)
    # if ( ret != 0 ):
    #     raise ValueError(f"ERROR in running {cmd2a}, returning exit code {ret}")

    ################STEP2b##################

    # args.pos1stSeq=os.getcwd()+"/spatialcoordinates.txt"
    # args.hdmi2ndSeq=os.getcwd()+"/HDMI_SeqScope_2nd.txt"
    # args.maxScale=2
    args.outpath=os.getcwd()
    # cmd2b="{args.py} estimateTissueBoundary/estimateTissueBoundary.py {args.pos1stSeq} {args.hdmi2ndSeq} {args.maxScale} {args.outpath} ".format(args=args)
    # ret = os.system(cmd2b)
    # if ( ret != 0 ):
    #     raise ValueError(f"ERROR in running {cmd2b}, returning exit code {ret}")
    
    ################STEP3a##################
    # print('Start Simple Gridding')
    # args.seqscope1st="MiSeq"
    # args.DGEdir=args.outpath+"Solo.out/GeneFull/raw/"
    # args.nrow=1
    # args.ncol=1
    # args.collapsePath=args.STtools+"/getSimpleGrid/colalpse.cpp"



    # cmd3a="Rscript getSimpleGrid/simpleGrid_v2.R {args.seqscope1st} {args.DGEdir} {args.spatial} {args.tiles} {args.nrow} {args.ncol} {args.sidesize} {args.outpath}".format(args=args)
    # ret = os.system(cmd3a)
    # if ( ret != 0 ):
    #     raise ValueError(f"ERROR in running {cmd3a}, returning exit code {ret}")
    
    ################STEP3b##################
    print("Start Sliding Window Gridding")
    print('Start Simple Gridding')

    args.DGEdir=args.outpath+"Solo.out/GeneFull/raw/"
    args.nrow=1
    args.ncol=1
    args.collapsePath=args.STtools+"/getSimpleGrid/colalpse.cpp"
    args.slidingPath=args.STtools+"/getSlidingGrid/slidingWindow.cpp"

    #cmd3b=cat tiles_input.txt | xargs -I {} -P 2 bash -c "Rscript getSlidingGrids_v2.R $seqscope1st $DGEdir $spatial $nrow $ncol $sidesize $outpath $window $collapsePath $slidingPath {} "


    args.input=os.getcwd+'/tiles_input.txt'
    print(args.input)
    cmd3b="cat {args.input} | xargs -I {} -P 2 bash -c "Rscript getSlidingGrid/getSlidingGrid_V2.R {args.seqscope1st} {args.DGEdir} {args.spatial} {args.nrow} {args.ncol} {args.sidesize} {args.outpath} {args.collapsePath} {args.slidingPath} {}"".format(args=args)
    ret = os.system(cmd3b)
    if ( ret != 0 ):
        raise ValueError(f"ERROR in running {cmd3b}, returning exit code {ret}")
    
    


    #print('Start running Step 1')


    #if ( args.star_path is None ):
     #   raise ValueError("Cannot find --star-path argument. Need to specify for running all")






# if ( args.run_steps is not None ):
#     steps = [int(x) for x in args.run_steps.split(',')]
#     print("Running the following steps: ", steps)


# if ( 1 in steps ): ## run step 1
#     print("Running step 1")
#     ## check whether parameter is given
#     if ( args.fq1 is None ):
#         raise ValueError("Cannot find --fq1 argument. Need to specify for running step 1")
#     with sp.Popen(f"zcat {args.fq1}", shell=True, encoding='utf-8', stdout=sp.PIPE).stdout as fh:
#         with sp.Popen(f"gzip -c > {args.out}.seq.gz", shell=True, encoding='utf-8', stdin=sp.PIPE).stdin as wh:
#             lineno = 0
#             for line in fh:
#                 if ( lineno % 4 == 1 ):
#                     wh.write(line)
#                 lineno += 1

# if ( 2 in steps ): ## run step2
#     print("Running step 2")
#     if ( args.cmd is None ):
#         raise ValueError("Cannot find --cmd argument. Need to specify for running step 2")    
#     ret = os.system(args.cmd)
#     if ( ret != 0 ):
#         raise ValueError(f"ERROR in running {args.cmd}, returning exit code {ret}") 

    
