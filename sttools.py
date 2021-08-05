def str2bool(value):
    if value.lower() in {'false', 'f', '0', 'no', 'n'}:
        return False
    elif value.lower() in {'true', 't', '1', 'yes', 'y'}:
        return True
    raise ValueError(f'{value} is not a valid boolean value')

def import_or_install(package):
    try:
        #print('try')
        __import__(package)
    except ImportError:
        #print('install')
        pip.main(['install', package])

#modules
import pip
modules = ["argparse","os",'sys','subprocess','math']
#map(import_or_install,modules)
md=[import_or_install(i) for i in modules]
import argparse
import os, sys
import subprocess as sp
import math



parser = argparse.ArgumentParser(description="STtools Main Launcher")
parser.add_argument("--run-all", default=False, action='store_true', help="Run all steps sequentially")
parser.add_argument("--run-steps", type=str, help="Comma-separated string containing the steps to run")
parser.add_argument("--first-fq",type=str, help="Path to 1nd-Seq FASTQ.gz file")
parser.add_argument("-l","--hdmilength",type=int,default=20, help="An integer indicating the length of the HDMIs") #if not specify -l erro None occurs
parser.add_argument("-s", "--star-path", type=str, default="STAR", help="Binary folder path to STAR executable") #if not given use current directory
parser.add_argument("-q", "--seqtk-path", type=str, default="seqtk", help="Binary folder path to seqtk executable")
parser.add_argument("-g", "--genome", type=str, help="Path to the reference files")
parser.add_argument("--second-fq1", type=str, help="Path to 2nd-Seq Read 1 FASTQ.gz file")
parser.add_argument("--second-fq2", type=str, help="Path to 2nd-Seq Read 2 FASTQ.gz file")
parser.add_argument("-o","--outprefix" ,type=str, help="Prefix for STARsolo output",default='Sample')
parser.add_argument("--STtools",type=str,help='Path to the scripts of STtools scripts',default='STtools')
parser.add_argument("--py",type=str,default=sys.executable,help='Binary path to python executable')
parser.add_argument("--lane-tiles",type=str,help='lane and tiles')
parser.add_argument("--binsize",type=int,help='size of the side of  square grids',default=300)
parser.add_argument("--window",type=int,help='sliding size of sliding window',default=150)
parser.add_argument("-c", "--cores",type=int,help='number of cores',default=5)
parser.add_argument("-m", "--maxScale",type=float,help='max number for color bar of HDMI density plot')
parser.add_argument("--DGEdir",type=str,help='Directory of digital expression matrix for griding')
parser.add_argument("--subDGEdir",type=str,help='Directory of digiral expression matrix containing spliced and unspliced reads for subcellular analysis')
#parser.add_argument("--subDGE",type=str,help='Directory of digiral expression matrix containing spliced and unspliced reads for subcellular analysis')
parser.add_argument("--hdmi2ndSeq",type=str,help='Path to .txt file with HDMIs from SeqScope 2nd-Seq')
parser.add_argument("-a","--alpha",type=float,help='transparency when ploting subCelluar images')
parser.add_argument("--spatial",type=str,help='Path to .txt file of spatial coordinates. Please see readme for file format')
parser.add_argument("--whitelist",type=str,help='Path to .txt file of whitelist for alignment. Please see xxx for file format')
parser.add_argument("--nFeaturePlotOnly",type=str,help='Plot feature plot only during clustering and mapping',default='FALSE')
parser.add_argument("--geneCount1",type=float,help='cutoff of nFeature_Spatial of Seurat spatial object of simple square grids for clustering',default=0)
parser.add_argument("--geneCount2",type=float,help='cutoff of nFeature_Spatial of Seurat spatial object of sliding square grids for clustering',default=0)
parser.add_argument("--simpleGridsPath",type=str,help='path to rds file of simple squrae grids')
parser.add_argument("--slidingGridsPath",type=str,help='path to rds file of sliding squrae grids')
parser.add_argument("--outdir",type=str,help='path to the output directory, if not specified, output in the current working directory')
#parser.add_argument("--nrow",type=int,help='number of rows when stacking the tiles')
parser.add_argument("--VISIUMpath",type=str,help='path to VISIUM data folder with ontaining the spatial/ and filtered_feature_bc_matrix/ subdirectories')
parser.add_argument("--seed",type=float, help="seed for running bayespace on visium data")
parser.add_argument("--nPCs",type=float, help='number of PCs used for running bayespace on visium data')
parser.add_argument("--nHVGs",type=float, help='number of highly variable genes for running bayespace on visium data')
parser.add_argument('--logNormalize',type=str2bool,default=True,help='boolean value indicating wheather to log normalize data or not before running bayespace on visium data')
parser.add_argument('--enhancedRes',type=str2bool,default=False, help='boolean valud indicating to use enhance clusterin or not for running bayespace on visium data')
parser.add_argument('--nCluster',type=float,help='number of expected clusters for visium data')
parser.add_argument('--nrep',type=float, help='number of repetion for running bayespace, it is suggested that at least 10000 is needed for a full data')
parser.add_argument('--datasource',type=str,help='string of the data source, for example: SeqScope, VISIUM, etc.')
parser.add_argument('--layout',type=str,help='path of layout files of super tiles')
parser.add_argument('--order',type=str,help='either bottom(bottom tiles at bottom) or top(bottom tiles at top)')

#Functions
def step1():

    print('Start running Step 1: extract coordinates')
    #check files
    if(os.path.isfile(args.first_fq)==False):
        raise ValueError("File --first-fq does not exist")
    if(os.path.isfile(args.second_fq1)==False):
        raise ValueError("File --second-fq1 does not exist")
    if(os.path.isdir(args.STtools)==False):
        raise ValueError("Directory --STtools does not exist")
    if ( args.hdmilength is None ):
        args.hdmilength=20
    if(args.outdir is None):
        args.outdir=os.getcwd()
    if(os.path.isdir(args.outdir)==False):
        raise ValueError("Directory --outdir does not exist")
    
    args.extract=args.STtools+"/extractCoord/extractCoord_v2.sh"
    cmd1 = "bash {args.extract} -m {args.first_fq} -h {args.second_fq1} -l {args.hdmilength} -o {args.outdir}".format(args=args)
    ret = os.system(cmd1)
    if ( ret != 0 ):
        raise ValueError(f"ERROR in running {cmd1}, returning exit code {ret}")


    no_miseq=sp.getoutput("wc -l ./HDMIs-MiSeq-temp.txt")
    #print(no_miseq)
    no_miseq_uniq=sp.getoutput("wc -l ./whitelist.txt")
    f=open("./summary_step1.txt","w")
    f.write('MiSeq reads:'+str(no_miseq)+" \n")
    f.write('Unique MiSeq:' +str(no_miseq_uniq) +' \n')
    f.close()    



def step2():
    print("Tissue boundary estimation")
    if(args.outdir is None):
        args.outdir=os.getcwd()
    if(os.path.isdir(args.outdir)==False):
        raise ValueError("Directory --outdir does not exist")

    if(args.hdmi2ndSeq is None):
        args.hdmi2ndSeq=args.outdir+"/HDMI_SeqScope_2nd.txt"
    if(args.spatial is None):
        args.spatial=args.outdir+"/spatialcoordinates.txt"
    if(args.maxScale is None):
        args.maxScale='0'
    if(os.path.isfile(args.spatial)==False):
        raise ValueError("Spatial coordinates from STEP 1 output does not exist")
    if(os.path.isfile(args.hdmi2ndSeq)==False):
        raise ValueError("HDMIs from 2nd-Seq from STEP 2 output does not exist")
    #args.maxScale=2
    #if(args.outdir is None):
     #   args.outdir=os.getcwd()
    print(args.maxScale)
    args.estTB=args.STtools+"/estimateTissueBoundary/estimateTissueBoundary.py"
    cmd3="{args.py} {args.estTB} {args.spatial} {args.hdmi2ndSeq} {args.maxScale} {args.outdir} ".format(args=args)
    ret = os.system(cmd3)
    if ( ret != 0 ):
        raise ValueError(f"ERROR in running {cmd3}, returning exit code {ret}")






    
def step3():

    print('Start running Step 2: STARsolo alignment')
    if ( args.second_fq1 is None ):
        raise ValueError("Cannot find --second-fq1 argument. Need to specify for running all")
    if ( args.second_fq2 is None ):
        raise ValueError("Cannot find --second-fq2 argument. Need to specify for running all")
    if ( args.genome is None ):
        raise ValueError("Cannot find --genome argument. Need to specify for running all")

#    if(args.whitelist is None):
 #      args.whitelist=os.getcwd()+"/whitelist.txt"
  #  if(os.path.isfile(args.whitelist)==False):
   #      raise ValueError("Whitelist from STEP 1 output does not exist")
    #if(os.path.isfile(spatial)==False):
     #    raise ValueError("Spatial coordinates from STEP 1 output does not exist")

    if ( args.hdmilength is None ):
        args.hdmilength=20
    if(os.path.isdir(args.star_path)==False):
        raise ValueError("Directory --star-path does not exist")
    if(os.path.isdir(args.seqtk_path)==False):
        raise ValueError("Directory --seqtk-path does not exist")
    if(os.path.isdir(args.STtools)==False):
        raise ValueError("Directory --STtools does not exist")
    if(args.outdir is None):
        args.outdir=os.getcwd()
    if(os.path.isdir(args.outdir)==False):
        raise ValueError("Directory --outdir does not exist")

    if(args.whitelist is None):
       args.whitelist=args.outdir+"/whitelist.txt"
    if(os.path.isfile(args.whitelist)==False):
         raise ValueError("Whitelist from STEP 1 output does not exist")
    #if(os.path.isfile(spatial)==False):
     #    raise ValueError("Spatial coordinates from STEP 1 output does not exi\
#st")


    args.align=args.STtools+"/align/align_v2.sh"
    print(args.outdir)
    cmd2 = "bash {args.align} -a {args.second_fq1} -b {args.second_fq2} -l {args.hdmilength} -w {args.whitelist} -o {args.outprefix} -t {args.star_path} -q {args.seqtk_path} -g {args.genome} -d {args.outdir}".format(args=args) 
    ret = os.system(cmd2)
    if ( ret != 0 ):
        raise ValueError(f"ERROR in running {cmd2}, returning exit code {ret}")


    #with open(args.outprefix+'Solo.out/GeneFull/Summary.csv','r') as firstfile, open('summary_step3.txt','a') as secondfile:
        #for line in firstfile:
               
             # append content to second file
            # secondfile.write(line)

   # no_hiseq=sp.getoutput("wc -l ./HDMI_SeqScope_2nd.txt")
    #print(no_miseq)
  #  args.out=os.getcwd()+"/"+args.outprefix+"Solo.out/GeneFull/Summary.csv"
 #   print(args.out)
#    os.system("scp {args.out} summary_step3.txt")
    #print(out)
    #print(args.out)
    
    #starsolo=sp.getoutput('less {args.out}')
    #f=open("summary_step3.txt","w")
    #f.write('HiSeq reads:'+str(no_hiseq)+" \n")
    #f.write('StarSolo summmaries:' +' \n')
    #f.write((starsolo))
    #f.close()
    #copyfile(args.out, 'summary_step3.txt')

    ##no_hiseq=sp.getoutput("wc -l HDMI_SeqScope_2nd.txt")
   # print(no_hiseq)
   # out=os.getcwd()+"/"+args.outprefix+"Solo/GeneFull/Summary.csv"
   # print(out)
   # args.out=out
   # print(args.out)
   # starsolo=sp.getoutput('less {args.out}')
   # f=open("./summary_step3.txt","w")
   # f.write('HiSeq reads:'+str(no_hiseq)+" \n")
   # f.write('StarSolo summmaries:' +' \n')
   # f.write((starsolo))
    #f.close()




def step4():
    print('Start Simple Gridding')
    #if(args.nrow is None):
     #   args.nrow=2
      #  tiles_vec=args.tiles.split(',')
      #  args.ncol=math.ceil(len(tiles_vec)/args.nrow)


    
    args.nrow=1
    args.ncol=1
    args.collapsePath=args.STtools+"/getSimpleGrid/collapse.cpp"
#    args.outpath=os.getcwd()
    if(args.binsize is None):
        args.binsize=300
    if(args.outdir is None):
        args.outdir=os.getcwd()
    if(os.path.isdir(args.outdir)==False):
        raise ValueError("Directory --outdir does not exist")

    if(args.DGEdir is None):
        args.DGEdir=args.outdir+'/'+args.outprefix+"Solo.out/GeneFull/raw/"
    if(os.path.isdir(args.DGEdir)==False):
        raise ValueError("Directory --DGEdir does not exist")
    if(('seqscope1st' in vars(args))==False):
        args.seqscope1st='MiSeq'
    if(os.path.isfile(args.collapsePath)==False):
        raise ValueError("File --collapsePath does not exist")
    if(args.spatial is None):
        args.spatial=args.outdir+"/spatialcoordinates.txt"
    if(args.order is None):
        args.order='top'
    if(args.layout is None):
        args.layout='FALSE'

    
    # if(args.outdir is None):
    #    args.outdir=os.getcwd()
    #if(os.path.isdir(args.outdir)==False):
     #   raise ValueError("Directory --outdir does not exist")
    #if(args.lanes is None):
        
    args.simple=args.STtools+"/getSimpleGrid/simpleGrid_v3.R"
    #cmd4="Rscript {args.simple} {args.seqscope1st} {args.DGEdir} {args.spatial} {args.tiles} {args.nrow} {args.ncol} {args.binsize} {args.outdir} {args.collapsePath}".format(args=args)
    print('here')
    cmd4="Rscript {args.simple} {args.seqscope1st} {args.DGEdir} {args.spatial} {args.lane_tiles} {args.nrow} {args.ncol} {args.binsize} {args.outdir} {args.collapsePath} {args.layout} {args.order}".format(args=args)

    ret = os.system(cmd4)
    if ( ret != 0 ):
        raise ValueError(f"ERROR in running {cmd4}, returning exit code {ret}")
    #args.bc=args.outdir+'/collapsedBarcodes.csv'
    #args.gene=args.outdir+'/collapsedGenes.csv'
    #no_bc=sp.getoutput("wc -l collapsedBarcodes.csv")
    #no_gene=sp.getoutput("wc -l collapsedGenes.csv")

    #print(no_bc)
    #print(args.outdir)
   # fout=args.outdir+"/summary_step4.txt"
   # f=open(fout,"w")
   # f.write('Total number of collapsed grids:'+str(no_bc)+" \n")
   # f.write('Total number of genes:' +str(no_gene) +' \n')
   # f.close()

def step5():
    print("Start Sliding Window Gridding")
    args.nrow=1
    args.ncol=1
    args.collapsePath=args.STtools+"/getSimpleGrid/collapse.cpp"
    #args.slidingPath=args.STtools+"/getSlidingGrid/slidingWindow.cpp"
 #   args.outpath=os.getcwd()
    if(args.outdir is None):
        args.outdir=os.getcwd()
    if(os.path.isdir(args.outdir)==False):
        raise ValueError("Directory --outdir does not exist")
    if(args.DGEdir is None):
        args.DGEdir=args.outdir+'/'+args.outprefix+"Solo.out/GeneFull/raw/"
    if(os.path.isdir(args.DGEdir)==False):
        raise ValueError("Directory --DGEdir does not exist")
    if(os.path.isfile(args.collapsePath)==False):
        raise ValueError("File --collapsePath does not exist")
    #if(os.path.isfile(args.slidingPath)==False):
     #   raise ValueError("File --slidingPath does not exist")
    if(('seqscope1st' in vars(args))==False):
        args.seqscope1st='MiSeq'
    if(args.spatial is None):
        args.spatial=args.outdir+"/spatialcoordinates.txt"
    if(os.path.isfile(args.spatial)==False):
        raise ValueError("File --spatial does not exist")
    if(args.cores is None):
        args.cores=5
        
    if(args.window is None):
        args.window=150
    if(args.binsize is None):
        args.binsize=300
    if(args.lane_tiles is None):
        raise ValueError("Lane and tiles are required")
    if(args.order is None):
        args.order='top'
    if(args.layout is None):
        args.layout='FALSE'
    #tiles_vec =args.tiles.split(',')
   # inputF=os.getcwd()+'/inputTiles.txt'
   # with open(inputF, 'w') as f:
   #     for l in tiles_vec:
   #         f.write(l+'\n')
   # f.close()
   # args.input=inputF

    #args.input=os.getcwd()+'/input.txt'
    #print(args.input)
    args.sliding_P1=args.STtools+"/getSlidingGrid/slidingGrid_P1_V1.R"
    args.sliding_P2=args.STtools+"/getSlidingGrid/slidingGrid_P2_V1.R"
    args.sliding_P3=args.STtools+"/getSlidingGrid/slidingGrid_P3_V1.R"
    if(os.path.isfile(args.sliding_P1)==False):
        raise ValueError("File --sliding_P1 does not exist")
   
    if(os.path.isfile(args.sliding_P2)==False):
        raise ValueError("File --sliding_P2 does not exist")
    
    if(os.path.isfile(args.sliding_P3)==False):
        raise ValueError("File --sliding_P3 does not exist")
    
    cmd5_1="Rscript {args.sliding_P1} {args.seqscope1st} {args.DGEdir} {args.spatial} {args.outdir} {args.lane_tiles}".format(args=args)
    ret = os.system(cmd5_1)
    if ( ret != 0 ):
        raise ValueError(f"ERROR in running {cmd5_1}, returning exit code {ret}")

    args.input=args.outdir+'/groupgrids_tile.txt'
    cmd5_2='cat {args.input} | xargs -I [] -P {args.cores} bash -c "Rscript {args.sliding_P2} {args.collapsePath} {args.DGEdir} {args.outdir} {args.window} {args.binsize} []"'.format(args=args)
    ret = os.system(cmd5_2)
    if ( ret != 0 ):
        raise ValueError(f"ERROR in running {cmd5_2}, returning exit code {ret}")

    cmd5_3 = "Rscript {args.sliding_P3} {args.outdir} {args.ncol} {args.nrow} {args.layout} {args.order} {args.lane_tiles}".format(args=args)
    ret = os.system(cmd5_3)
    if ( ret != 0 ):
        raise ValueError(f"ERROR in running {cmd5_3}, returning exit code {ret}\
")



def step6():
    print('Start clustering and mapping')
    if(args.outdir is None):
        args.outdir=os.getcwd()
    if(os.path.isdir(args.outdir)==False):
        raise ValueError("Directory --outdir does not exist")


   # args.workingdir = os.getcwd()
    if(args.simpleGridsPath is None):
        args.simpleGridsPath = args.outdir+'/SimpleSquareGrids.RDS'
    
    if(args.slidingGridsPath is None):
        args.slidingGridsPath = args.outdir+'/SlidingSquareGrids.RDS'
    if(args.geneCount1 is None):
        args.geneCount1=0
    if(args.geneCount2 is None):
        args.geneCount2=0
    if(args.nFeaturePlotOnly is None):
        args.nFeaturePlotOnly = 'FALSE'
    print(args.slidingGridsPath)
    print(args.simpleGridsPath)
    if(os.path.isfile(args.slidingGridsPath)==False):
        raise ValueError("Path --slidingGridsPath does not exist")
    if(os.path.isfile(args.simpleGridsPath)==False):
        raise ValueError("Path --simpleGridsPath does not exist")
#    if(args.outdir is None):
 #       args.outdir=os.getcwd()
  #  if(os.path.isdir(args.outdir)==False):
   #     raise ValueError("Directory --outdir does not exist")

    #print(args.slidingGridsPath)
    args.clus=args.STtools+'/clusterAndMap/clusterAndMap.R'
    
    cmd = "Rscript {args.clus} {args.outdir} {args.simpleGridsPath} {args.slidingGridsPath} {args.geneCount1} {args.geneCount2} {args.nFeaturePlotOnly}".format(args=args)
    ret = os.system(cmd)
    if(ret!=0):
        raise ValueError (f"ERROR in running {cmd}, returning exit code {ret}")


    
       
def step7():
    print('Start subcellular analysis!')
    if(args.outdir is None):
        args.outdir=os.getcwd()
    if(os.path.isdir(args.outdir)==False):
        raise ValueError("Directory --outdir does not exist")

    if(args.subDGEdir is None):
        args.subDGEdir=args.outdir+"/"+args.outprefix+"Solo.out/Velocyto/raw/"
    print(args.subDGEdir)
    if(os.path.isdir(args.subDGEdir)==False):
        raise ValueError("Directory --subDGEdir does not exist")
    if(args.spatial is None):
        args.spatial=args.outdir+"/spatialcoordinates.txt"
    if(os.path.isfile(args.spatial)==False):
        raise ValueError("File --spatial does not exist")
    args.subana=args.STtools+"/subCellularAna/subCellularAna_v2.py"
    #args.workingdir=os.getcwd()
    if(('seqscope1st' in vars(args))==False):
        args.seqscope1st='MiSeq'  
    if(args.tiles is None):
        raise ValueError("Tiles are required")
    if(args.alpha is None):
        args.alpha=0.01
   # args.tiles_vec =args.tiles.split(',')
    #if(args.outdir is None):
     #   args.outdir=os.getcwd()
    #if(os.path.isdir(args.outdir)==False):
     #   raise ValueError("Directory --outdir does not exist")


    #print(args.py)
    cmd7="{args.py} {args.subana} {args.subDGEdir} {args.outdir} {args.spatial} {args.seqscope1st} {args.tiles} {args.alpha}".format(args=args)
    ret = os.system(cmd7)
    if ( ret != 0 ):
        raise ValueError(f"ERROR in running {cmd7}, returning exit code {ret}")



args = parser.parse_args()
steps = []
if ( args.run_steps is not None ):
    steps = [int(x) for x in args.run_steps.split(',')]
    n_steps=len(steps)
    print(n_steps)
    s_steps=sorted(steps)
  #if only one step
    if (set(steps).issubset(set(range(1,8))) ==False):
        print('Steps must be in the range from 1 to 7!')
    if (n_steps==1):
        print('Run step', steps[0])
        func=eval("step"+str(steps[0]))
        func()
    elif (s_steps == list(range(min(steps), max(steps)+1))):
        print('Run the following consecutive steps', s_steps)
        for i in s_steps:
            print(i)
            func=eval("step"+str(i))
            func()
    else:
        print('Please provide consecutive steps!')

#     print("Running the following steps: ", steps)
if ( args.run_all ):
   # steps = [1,2,3,4,5]
    print("Running the following steps: ", steps)

    ## check whether parameter is given
    if ( args.first-fq is None ):
        raise ValueError("Cannot find --first-fq argument. Need to specify for running all")
    if ( args.second-fq1 is None ):
        raise ValueError("Cannot find --second-fq1 argument. Need to specify for running all")
    if ( args.second-fq2 is None ):
        raise ValueError("Cannot find --second-fq2 argument. Need to specify for running all")
    if ( args.genome is None ):
        raise ValueError("Cannot find --genome argument. Need to specify for running all")

    #check file/directory path exist:
    if(os.path.isfile(args.first-fq)==False):
        raise ValueError("File --first-fq does not exist")
    if(os.path.isfile(args.second-fq1)==False):
        raise ValueError("File --second-fq1 does not exist")
    if(os.path.isfile(args.second-fq2)==False):
        raise ValueError("File --second-fq2 does not exist")
    if(os.path.isdir(args.star_path)==False):
        raise ValueError("Directory --star-path does not exist")
    if(os.path.isdir(args.seqtk_path)==False):
        raise ValueError("Directory --seqtk-path does not exist")
    if(os.path.isdir(args.genome)==False):
        raise ValueError("Directory --genome does not exist")
    if(os.path.isdir(args.STtools)==False):
        raise ValueError("Directory --STtools does not exist")

    for i in range(1,8):
            print('Running step:' + i)
            func=eval("step"+str(i))
            func()
