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
parser.add_argument("--geneCount1",type=float,help='cutoff of nFeature_Spatial of Seurat spatial object of simple square grids for clustering')
parser.add_argument("--geneCount2",type=float,help='cutoff of nFeature_Spatial of Seurat spatial object of sliding square grids for clustering')
parser.add_argument("--simpleGridsPath",type=str,help='path to rds file of simple squrae grids')
parser.add_argument("--slidingGridsPath",type=str,help='path to rds file of sliding squrae grids')
parser.add_argument("--outdir",type=str,help='path to the output directory, if not specified, output in the current working directory')
#parser.add_argument("--nrow",type=int,help='number of rows when stacking the tiles')
parser.add_argument("--VISIUMpath",type=str,help='path to VISIUM data folder with ontaining the spatial/ and filtered_feature_bc_matrix/ subdirectories')
parser.add_argument("--seed",type=float, help="seed for running bayespace on visium data")
parser.add_argument("--nPCs",type=float, help='number of PCs used for running bayespace on visium data')
parser.add_argument("--nHVGs",type=float, help='number of highly variable genes for running bayespace on visium data')
parser.add_argument('--logNormalize',type=str2bool,default=True,help='boolean value indicating wheather to log normalize data or not before running bayespace on visium data')
parser.add_argument('--enhancedRes',type=str2bool,default=False, help='boolean value indicating to use enhance clusterin or not for running bayespace on visium data')
parser.add_argument('--nCluster',type=float,help='number of expected clusters for visium data')
parser.add_argument('--nrep',type=float, help='number of repetion for running bayespace, it is suggested that at least 10000 is needed for a full data')
parser.add_argument('--datasource',type=str,help='string of the data source, for example: SeqScope, VISIUM, etc.')
parser.add_argument('--layout',type=str,help='path of layout files of super tiles')
parser.add_argument('--order',type=str,help='either bottom(bottom tiles at bottom) or top(bottom tiles at top)')
#parser.add_argument('--sliding',type=str,type=str2bool,default=True,help='boolean indicating wheather sliding window is run or not')
parser.add_argument('--algo',type=str,help=' A string indicating what clustering method is employed. Is is suggested that either Seurat or Bayespace for VISIUM data, slidingWindow for SlideSeq and SeqScope data')
parser.add_argument('--res',type=float,help='resolution for Seurat clustering')
parser.add_argument('--clustering',type=str2bool,help='wheather clustering for simplesquaregrids, default is False')
parser.add_argument('--seqscope1st',type=str,help='either MiSeq or HiSeq or Custom')
parser.add_argument('--annotatedSimpleGrids',type=str2bool,help='indecate weather cluster is annoted by user or not. If yes, will not run clustering in Step6, otehrwise, run clustering with simple steps')





parser.add_argument("--skip-errors-after", type=int, default=100, help="Skip printing error messages after displaying a certain number of times")
parser.add_argument("--sorted-coo", default=False, action='store_true', help="Indicates that coordinate file is already sorted")
#outdir:parser.add_argument("-o", "--out", type=str, required=True, help="Output directory")
parser.add_argument("--pigz", type=str, default="pigz", help="Path to pigz binary")
parser.add_argument("--sort", type=str, default="sort", help="Path to sort binary (with options if needed)")
parser.add_argument("--ncpus", type=int, default=4, help="Number of CPUs to use for parallelized jobs")
parser.add_argument("--visualizelayout", type=str, help="Layout file of tiles to draw [lane] [tile] [row] [col] format in each line")
parser.add_argument("--predir", type=str,  help="Input (STTools output) directory")
parser.add_argument("--scale", type=float, default=20.0, help="Scale each color to have the same mean intensity")
parser.add_argument("--inv-weight", default=False, action='store_true', help="Weight each gene inverse")
parser.add_argument("--min-tpm", default=1000, type=float, help="Minimum TPM value for inverse weighting (effective only with --inv-weight)")
parser.add_argument("--red", type=str, help="Comma-separate list of genes (colon with weights) for red color")
parser.add_argument("--green",type=str, help="Comma-separate list of genes (colon with weights) for green color")
parser.add_argument("--blue", type=str, help="Comma-separate list of genes (colon with weights) for blue color")
parser.add_argument("--pixel", type=int, default=80, help="Resolution of pixel (how to bin each pixel) - default: 80 (um2 per pixel)")
parser.add_argument("--outfilePrefix", type=str, help="Output file prefix")
parser.add_argument("--tmpdir", type=str, help="temporary directory")
parser.add_argument("--buffer-size", type=str, help="Buffer size for sorting DGE")
parser.add_argument("--nFeature-cutoff",type=float,help='cutoff of nFeature_Spatial of Seurat object for visualization V2')
parser.add_argument("--objpath",type=str,help='path of Seurat Object for visualization V2')
#Functions
def stepA1():

    #print('Start running Step 1: extract coordinates')
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
    
    args.extract=args.STtools+"/extractCoord/extractCoord.sh"
    cmd1 = "bash {args.extract} -m {args.first_fq} -h {args.second_fq1} -l {args.hdmilength} -o {args.outdir}".format(args=args)
    ret = os.system(cmd1)
    if ( ret != 0 ):
        raise ValueError(f"ERROR in running {cmd1}, returning exit code {ret}")


    #no_miseq=sp.getoutput("wc -l ./HDMIs-MiSeq-temp.txt")
    #print(no_miseq)
    #no_miseq_uniq=sp.getoutput("wc -l ./whitelist.txt")
    #f=open("./summary_step1.txt","w")
    #f.write('MiSeq reads:'+str(no_miseq)+" \n")
    #f.write('Unique MiSeq:' +str(no_miseq_uniq) +' \n')
    #f.close()    



def stepA2():
    #print("Tissue boundary estimation")
    if(args.outdir is None):
        args.outdir=os.getcwd()
    if(os.path.isdir(args.outdir)==False):
        raise ValueError("Directory --outdir does not exist")

    if(args.hdmi2ndSeq is None):
        args.hdmi2ndSeq=args.outdir+"/HDMI_SeqScope_2nd.txt.gz"
    if(args.spatial is None):
        args.spatial=args.outdir+"/spatialcoordinates.txt.gz"
    if(args.maxScale is None):
        args.maxScale='0'
    if(os.path.isfile(args.spatial)==False):
        raise ValueError("Spatial coordinates from STEP 1 output does not exist")
    if(os.path.isfile(args.hdmi2ndSeq)==False):
        raise ValueError("HDMIs from 2nd-Seq from STEP 2 output does not exist")
    #args.maxScale=2
    #if(args.outdir is None):
     #   args.outdir=os.getcwd()
    #print(args.maxScale)
    args.estTB=args.STtools+"/estimateTissueBoundary/estimateTissueBoundary.py"
    cmd2="{args.py} {args.estTB} {args.spatial} {args.hdmi2ndSeq} {args.maxScale} {args.outdir} ".format(args=args)
    ret = os.system(cmd2)
    if ( ret != 0 ):
        raise ValueError(f"ERROR in running {cmd2}, returning exit code {ret}")






    
def stepA3():

    #print('Start running Step 2: STARsolo alignment')
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


    args.align=args.STtools+"/align/align.sh"
    #print(args.outdir)

    cmd3 = "bash {args.align} -a {args.second_fq1} -b {args.second_fq2} -l {args.hdmilength} -w {args.whitelist} -o {args.outprefix} -t {args.star_path} -q {args.seqtk_path} -g {args.genome} -d {args.outdir}".format(args=args) 
    ret = os.system(cmd3)
    if ( ret != 0 ):
        raise ValueError(f"ERROR in running {cmd3}, returning exit code {ret}")
    args.orgDGE=args.STtools+"/align/merge-dge-hdmi.py"
    if(args.ncpus is None):
        args.ncpus=5
    if(args.DGEdir is None):
        args.DGEdir=args.outdir+'/'+args.outprefix+"Solo.out/GeneFull/raw/"
    args.predir=args.outdir+'/'+args.outprefix+"Solo.out/GeneFull/ordered"
    if not args.outdir:
        os.makedirs(args.outdir)

    print('temp')
    print(args.tmpdir)
    print(args.buffer_size)
    args.tmpdir='/new/tmp_dir'
    args.buffer_size='10G'
   # if (args.tmpdir is None) & (args.buffer_size is None):
    cmd3_5="{args.py} {args.orgDGE} -c {args.spatial}  -d {args.DGEdir} -o {args.predir} --ncpus {args.ncpus} --sort {args.sort} --tmpdir {args.tmpdir} --buffer-size {args.buffer_size}".format(args=args)
    #cmd3_5="{args.py} {args.orgDGE} -c {args.spatial}  -d {args.DGEdir} -o {args.predir} --ncpus {args.ncpus} --sort {args.sort}".format(args=args)
    print(cmd3_5)
    ret = os.system(cmd3_5)
    if ( ret != 0 ):
        raise ValueError(f"ERROR in running {cmd3_5}, returning exit code {ret}")

    #args.outdir=args.outdir+'/'+args.outprefix+"Solo.out/Velocyto/ordered"
    #if not args.outdir:
    #    os.makedirs(args.outdir)
    #cmd3_5="{args.py} {args.orgDGE} -c {args.spatial}  -d {args.DGEdir} -o {arg\
#s.outdir} --ncpus {args.ncpus} --sort {args.sort}".format(args=args)
 #   print(cmd3_5)
  #  ret = os.system(cmd3_5)
   # if ( ret != 0 ):
    #    raise ValueError(f"ERROR in running {cmd3_5}, returning exit code {ret}")

def stepC1():
    if(args.datasource is None):
       args.datasource = 'SeqScope'
       args.nMax=100
    if(args.datasource =='SlideSeq'):
       args.nMax=None 
    if(args.datasource =='VISIUM'):
        print('VIUSIUM DATA')
        if(args.outdir is None):
            args.outdir=os.getcwd()
        if(os.path.isdir(args.outdir)==False):
            raise ValueError("Directory --outdir does not exist")

        if(args.algo is None):
            args.algo='Bayespace'

        print('stop 2')
        if(args.algo=='Bayespace'):
            print('Starting Clustering using Bayespace')
            if(os.path.isdir(args.DGEdir)==False):
                raise ValueError("Directory --DGEdir does not exist")
            if(args.spatial is None):
                args.spatial=os.getcwd()+"/spatialcoordinates.txt.gz"
            if(args.seed is None):
                args.seed=1234
            if (args.nPCs is None):
                args.nPCs = 10
            if (args.nHVGs is None):
                args.nHVGs=2000
            if (args.logNormalize is None):
                args.logNormalize=True
            if (args.enhancedRes is None):
                args.enhancedRes = True
            if (args.nCluster is None):
                args.nCluster=4
            if (args.nrep is None):
                args.nrep=10000
            args.bayespace=args.STtools+"/bayespace/bayespace_v2.R"
            print(args.bayespace)
            cmd4="Rscript {args.bayespace} {args.DGEdir} {args.spatial} {args.outdir} {args.seed} {args.nPCs} {args.nHVGs} {args.logNormalize} {args.enhancedRes} {args.nCluster} {args.nrep} ".format(args=args)
            print (cmd4)
            ret = os.system(cmd4)
            if ( ret != 0 ):
                raise ValueError(f"ERROR in running {cmd4}, returning exit code {ret}")


        if(args.algo=='Seurat'):
            if(os.path.isdir(args.DGEdir)==False):
                raise ValueError("Directory --DGEdir does not exist")
            if(args.spatial is None):
                args.spatial=os.getcwd()+"/spatialcoordinates.txt.gz"

            if (args.nPCs is None):
                args.nPCs = 10
            if (args.res is None):
                args.res=0.5
            args.seurat=args.STtools+"/STSeurat/runSeurat.R"
            cmd4="Rscript {args.seurat} {args.DGEdir} {args.spatial} {args.outdir} {args.nPCs} {args.res}".format(args=args)
            print (cmd4)
            ret = os.system(cmd4)
            if ( ret != 0 ):
                raise ValueError(f"ERROR in running {cmd4}, returning exit code {ret}")

    if(args.datasource !='VISIUM'):    
        print('Start Simple Squre Gridding')
        args.nrow=1
        args.ncol=1
        args.collapsePath=args.STtools+"/getSimpleGrid/collapse.cpp"
        if (args.clustering is None):
            args.clustering = False  
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

        if(os.path.isfile(args.collapsePath)==False):
            raise ValueError("File --collapsePath does not exist")
        if(args.spatial is None):
            args.spatial=args.outdir+"/spatialcoordinates.txt.gz"
        if(args.layout is None):
            args.layout='MiSeq'
        if(args.lane_tiles is None):
            args.lane_tiles='All'

        args.simple=args.STtools+"/getSimpleGrid/simpleGrid.R"
        print(args.lane_tiles)
        print(args.outdir)
        cmd4="Rscript {args.simple}  {args.DGEdir} {args.spatial} {args.lane_tiles} {args.nrow} {args.ncol} {args.binsize} {args.outdir} {args.collapsePath} {args.layout} {args.nMax} {args.clustering}".format(args=args)
        ret = os.system(cmd4)
        if ( ret != 0 ):
            raise ValueError(f"ERROR in running {cmd4}, returning exit code {ret}")

def stepC2():
    if(args.datasource is None):
       args.datasource = 'SeqScope'
       args.nMax=100
    if(args.datasource =='VISIUM'):
       print('Not a step for VISIUM data')
       exit(1)
  #  if (args.datasource == 'SlideSeq' ):
   #     args.sliding=
    else:
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
        if(args.spatial is None):
            args.spatial=args.outdir+"/spatialcoordinates.txt.gz"
        if(os.path.isfile(args.spatial)==False):
            raise ValueError("File --spatial does not exist")
        if(args.cores is None):
            args.cores=5
            
        if(args.window is None):
            args.window=150
        if(args.binsize is None):
            args.binsize=300
       # if(args.lane_tiles is None):
       #     raise ValueError("Lane and tiles are required")
       # if(args.layout is None):
        #    args.layout='HiSeq'
        if(args.layout is None):
            args.layout='MiSeq'
        if(args.lane_tiles is None):
            args.lane_tiles='All'
        args.sliding_P1=args.STtools+"/getSlidingGrid/slidingGrid_P1.R"
        args.sliding_P2=args.STtools+"/getSlidingGrid/slidingGrid_P2.R"
        args.sliding_P3=args.STtools+"/getSlidingGrid/slidingGrid_P3.R"
        if(os.path.isfile(args.sliding_P1)==False):
            raise ValueError("File --sliding_P1 does not exist")
       
        if(os.path.isfile(args.sliding_P2)==False):
            raise ValueError("File --sliding_P2 does not exist")
        
        if(os.path.isfile(args.sliding_P3)==False):
            raise ValueError("File --sliding_P3 does not exist")
        
        cmd5_1="Rscript {args.sliding_P1} {args.layout} {args.DGEdir} {args.spatial} {args.outdir} {args.lane_tiles}".format(args=args)
        ret = os.system(cmd5_1)
        if ( ret != 0 ):
            raise ValueError(f"ERROR in running {cmd5_1}, returning exit code {ret}")

        args.input=args.outdir+'/groupgrids_tile.txt'
        cmd5_2='cat {args.input} | xargs -I [] -P {args.cores} bash -c "Rscript {args.sliding_P2} {args.collapsePath} {args.DGEdir} {args.outdir} {args.window} {args.binsize} []"'.format(args=args)
        ret = os.system(cmd5_2)
        if ( ret != 0 ):
            raise ValueError(f"ERROR in running {cmd5_2}, returning exit code {ret}")
        print('fnish cmd52')
        cmd5_3 = "Rscript {args.sliding_P3} {args.outdir} {args.ncol} {args.nrow} {args.layout} {args.order} {args.lane_tiles}".format(args=args)
        ret = os.system(cmd5_3)
        if ( ret != 0 ):
            raise ValueError(f"ERROR in running {cmd5_3}, returning exit code {ret}")



def stepC3():
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
        args.geneCount1=-1
    if(args.geneCount2 is None):
        args.geneCount2=-1
    if(args.annotatedSimpleGrids is None):
        args.annotatedSimpleGrids=False
    # if(args.nFeaturePlotOnly is None):
    #     args.nFeaturePlotOnly = 'FALSE'
    # print(args.slidingGridsPath)
    # print(args.simpleGridsPath)
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
    print('args.genecount1')
    print(args.geneCount1)
    cmd6 = "Rscript {args.clus} {args.outdir} {args.simpleGridsPath} {args.slidingGridsPath} {args.geneCount1} {args.geneCount2} {args.annotatedSimpleGrids}".format(args=args)
    print(cmd6)
    ret = os.system(cmd6)
    if(ret!=0):
        raise ValueError (f"ERROR in running {cmd6}, returning exit code {ret}")


    
       
def stepV1():
    #print(args.red)
    #print(args.blue)
    if (args.red is None):
       # args.red="Alb"
        raise ValueError("Please speficy  --red for visualization")
    if (args.green is None):
       # args.red="GAPDH"
      	raise ValueError("Please speficy  --green for visualization")
    if (args.blue is None):
       # args.blue="PGK1"
       	raise ValueError("Please speficy  --blue for visualization")
    #print(args.red)
    #print(args.blue)

    if(args.layout=='MiSeq'):
        args.visualizelayout=args.STtools+'visualization/miseq_layout.tsv'
    if(args.layout=='HiSeq'):
        args.visualizelayout=args.STtools+'visualization/hiseq_layout.tsv'
    if(args.lane_tiles is None):
          args.lane_tiles='All'
    print('Start Visualization!')
    if(args.outdir is None):
        args.outdir=os.getcwd()
    print('a')
    if(os.path.isdir(args.outdir)==False):
        raise ValueError("Directory --outdir does not exist")
    print('b')
    if(args.outfilePrefix is None):
        args.outfilePrefix=args.outdir+'/Sample_vis'
    if(args.DGEdir is None):
        args.DGEdir=args.outdir+'/'+args.outprefix+"Solo.out/GeneFull/raw/"
    
    if(args.predir is None):
        args.predir=args.outdir+'/'+args.outprefix+"Solo.out/GeneFull/ordered"
    if(args.ncpus is None):
        args.ncpus=5
    if(args.spatial is None):
        args.spatial=args.outdir+'/spatialcoordinates.txt.gz'
        print('here')
    if not os.path.exists(args.predir):
        print('inside')
        os.makedirs(args.predir)
        args.orgDGE=args.STtools+"/align/merge-dge-hdmi.py"
        #args.spatial=
        print('temp')
        print(args.tmpdir)
        print(args.buffer_size)
        cmd3_5="{args.py} {args.orgDGE} -c {args.spatial}  -d {args.DGEdir} -o {args.predir} --ncpus {args.ncpus} --sort {args.sort} --tmpdir {args.tmpdir} --buffer-size {args.buffer_size}".format(args=args)
        print(cmd3_5)
        ret = os.system(cmd3_5)
        if ( ret != 0 ):
             raise ValueError(f"ERROR in running {cmd3_5}, returning exit code {ret}")
        
    args.rgb=args.STtools+'visualization/rgb-gene-image.py'
    #args.visualizelayout=args.STtools+'/visualization/layout.tsv'
    print(args.predir)
    print(args.outfilePrefix)
    cmd7="{args.py} {args.rgb} -o {args.outfilePrefix} -d {args.predir}  -g {args.green} -b {args.blue} -r {args.red} --scale {args.scale}  --min-tpm {args.min_tpm} --res {args.pixel} --layout {args.visualizelayout}".format(args=args)
    ret = os.system(cmd7)
    if ( ret != 0 ):
        raise ValueError(f"ERROR in running {cmd7}, returning exit code {ret}")

def stepV2():
    print('Run V2')
    if args.objpath is None:
        args.objpath=args.outdir+'/simpleGridWithClustering.RDS'
    if(args.outdir is None):
        args.outdir=os.getcwd()
    if(args.nFeature_cutoff is None):
        args.nFeature_cutoff= -1
    args.customPlot=args.STtools+"/visualization/customPlot.R"
    cmd8="Rscript {args.customPlot}  {args.objpath} {args.nFeature_cutoff} {args.outdir}".format(args=args)
    print(cmd8)
    ret = os.system(cmd8)
    if ( ret != 0 ):
        raise ValueError(f"ERROR in running {cmd8}, returning exit code {ret}")
    print('ret')
    print(ret)

args = parser.parse_args()
steps = []
allsets=['A1','A2','A3','C1','C2','C3','V1','V2']
if ( args.run_steps is not None ):
    steps = [x for x in args.run_steps.split(',')]
    n_steps=len(steps)
    #print(n_steps)
    s_steps=sorted(steps)
    #print('sorted steps')
    #print(s_steps)
    ind_min=allsets.index(min(s_steps))
    ind_max=allsets.index(max(s_steps))
    matched_sets=[allsets[i] for i in range(ind_min,ind_max+1)]
  #if only one step
    if (set(steps).issubset(set(allsets)) ==False):
        print('Steps must among A1,A2,A3,C1.C2,C3,V1,V2!')
    if (n_steps==1):
        print('Run step', steps[0])
        func=eval("step"+str(steps[0]))
        func()
        #print(func)
    elif (s_steps ==matched_sets):
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
    print("Running all steps: ", steps)

    ## check whether parameter is given
    if ( args.first_fq is None ):
        raise ValueError("Cannot find --first-fq argument. Need to specify for running all")
    if ( args.second_fq1 is None ):
        raise ValueError("Cannot find --second-fq1 argument. Need to specify for running all")
    if ( args.second_fq2 is None ):
        raise ValueError("Cannot find --second-fq2 argument. Need to specify for running all")
    if ( args.genome is None ):
        raise ValueError("Cannot find --genome argument. Need to specify for running all")

    #check file/directory path exist:
    if(os.path.isfile(args.first_fq)==False):
        raise ValueError("File --first-fq does not exist")
    if(os.path.isfile(args.second_fq1)==False):
        raise ValueError("File --second-fq1 does not exist")
    if(os.path.isfile(args.second_fq2)==False):
        raise ValueError("File --second-fq2 does not exist")
    if(os.path.isdir(args.star_path)==False):
        raise ValueError("Directory --star-path does not exist")
    if(os.path.isdir(args.seqtk_path)==False):
        raise ValueError("Directory --seqtk-path does not exist")
    if(os.path.isdir(args.genome)==False):
        raise ValueError("Directory --genome does not exist")
    if(os.path.isdir(args.STtools)==False):
        raise ValueError("Directory --STtools does not exist")

    for i in allsets:
            print('Running step:' + str(i))
            func=eval("step"+str(i))
            func()



# if ( args.run_steps is not None ):
#     steps = [int(x) for x in args.run_steps.split(',')]
#     n_steps=len(steps)
#     #print(n_steps)
#     s_steps=sorted(steps)
#   #if only one step
#     if (set(steps).issubset(set(range(1,8))) ==False):
#         print('Steps must be in the range from 1 to 7!')
#     if (n_steps==1):
#         print('Run step', steps[0])
#         func=eval("step"+str(steps[0]))
#         func()
#     elif (s_steps == list(range(min(steps), max(steps)+1))):
#         print('Run the following consecutive steps', s_steps)
#         for i in s_steps:
#             print(i)
#             func=eval("step"+str(i))
#             func()
#     else:
#         print('Please provide consecutive steps!')

# #     print("Running the following steps: ", steps)
# if ( args.run_all ):
#    # steps = [1,2,3,4,5]
#     print("Running the following steps: ", steps)

#     ## check whether parameter is given
#     if ( args.first_fq is None ):
#         raise ValueError("Cannot find --first-fq argument. Need to specify for running all")
#     if ( args.second_fq1 is None ):
#         raise ValueError("Cannot find --second-fq1 argument. Need to specify for running all")
#     if ( args.second_fq2 is None ):
#         raise ValueError("Cannot find --second-fq2 argument. Need to specify for running all")
#     if ( args.genome is None ):
#         raise ValueError("Cannot find --genome argument. Need to specify for running all")

#     #check file/directory path exist:
#     if(os.path.isfile(args.first_fq)==False):
#         raise ValueError("File --first-fq does not exist")
#     if(os.path.isfile(args.second_fq1)==False):
#         raise ValueError("File --second-fq1 does not exist")
#     if(os.path.isfile(args.second_fq2)==False):
#         raise ValueError("File --second-fq2 does not exist")
#     if(os.path.isdir(args.star_path)==False):
#         raise ValueError("Directory --star-path does not exist")
#     if(os.path.isdir(args.seqtk_path)==False):
#         raise ValueError("Directory --seqtk-path does not exist")
#     if(os.path.isdir(args.genome)==False):
#         raise ValueError("Directory --genome does not exist")
#     if(os.path.isdir(args.STtools)==False):
#         raise ValueError("Directory --STtools does not exist")

#     for i in range(1,8):
#             print('Running step:' + str(i))
#             func=eval("step"+str(i))
#             func()
