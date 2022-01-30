import os, sys, gzip, math, random
import argparse
import subprocess as sp

parser = argparse.ArgumentParser(description="STtools Post-alignment pipeline")
parser.add_argument("-d", "--dge-dir", type=str, required=True, help="10x formatted directory")
parser.add_argument("-m", "--mtx", type=str, default="matrix.mtx", help="matrix file for DGE inside --dge-dir. Must be sorted by barcode index (default by STARsolo)")
parser.add_argument("-b", "--bcd", type=str, default="barcodes.tsv", help="barcode.tsv file for DGE inside --dge-dir. Must be sorted alphabetically (default by STARsolo)")
parser.add_argument("-f", "--ftr", type=str, default="features.tsv", help="features.tsv file for DGE inside --dge-dir")
parser.add_argument("-c", "--coo", type=str, required=True, help="spatialcoordinates.txt file for barcodes")
parser.add_argument("--skip-errors-after", type=int, default=100, help="Skip printing error messages after displaying a certain number of times")
parser.add_argument("--sorted-coo", default=False, action='store_true', help="Indicates that coordinate file is already sorted")
parser.add_argument("-o", "--out", type=str, required=True, help="Output directory")
parser.add_argument("--pigz", type=str, default="pigz", help="Path to pigz binary")
parser.add_argument("--sort", type=str, default="sort", help="Path to sort binary")
parser.add_argument("--tmpdir", type=str, help="temporary directory (-T) for sort")
parser.add_argument("--buffer-size", type=str, help="buffer size (-S) for sort")
parser.add_argument("--ncpus", type=int, default=4, help="Number of CPUs to use for parallelized jobs")

args = parser.parse_args()

## Sort the coordinate file
if not os.path.exists(args.out):
    os.makedirs(args.out)    

if not args.sorted_coo:
    ## need to sort coo files
    coofile = f"{args.out}/spatialcoordinates.sorted.txt.gz"
    sortcmd = args.sort + (f" -T{args.tmpdir}" if ( args.tmpdir is not None ) else "") + (f" -S{args.buffer_size}" if ( args.buffer_size is not None) else "")
    cmd = f"zcat -f {args.coo} | {sortcmd} --parallel {args.ncpus} | {args.pigz} -p {args.ncpus} > {coofile}"
    print(f"Running the following command:\n{cmd}\nNOTE: This make take a tens or minutes or even hours...\nIf this job fails or runs very slow, try using [--tmpdir /new/tmp_dir] and/or [--buffer-size 10G'] options instead")
    exitcode = os.system(cmd)
    if ( exitcode != 0 ):
        raise ValueError(f"The command failed with exit code {exitcode}")
else:
    coofile = args.coo

## Read by barcodes
lane2tile2iftr2cnts = {} ## lane -> tile -> iftr -> cnt
iftr2cnts = {}
lane2tile2info = {} ## lane -> tile -> [nbcds, nlines, mtxh, bcdh]

g_mtx_fh = open(f"{args.out}/.tmp.matrix.mtx",'wt',encoding='utf-8')
g_bcd_fh = gzip.open(f"{args.out}/barcodes.tsv.gz",'wt',encoding='utf-8')
g_ftr_fh = gzip.open(f"{args.out}/features.tsv.gz",'wt',encoding='utf-8')
len_cnts = 0

mtxs = args.mtx.split(',')
if ( len(mtxs) == 1 ):
    if mtxs[0].endswith("gz"):
        mtxcmd = f"zcat {args.dge_dir}/{mtxs[0]}"
    else:
        mtxcmd = f"cat {args.dge_dir}/{mtxs[0]}"
else:
    mtxcmd = "paste "
    for m in mtxs:
        if m.endswith("gz"):
            mtxcmd += f" <(zcat {args.dge_dir}/{m})"
        else:
            mtxcmd += f" {args.dge_dir}/{m}"
    if "<(" in mtxcmd :
        mtxcmd = "bash -c '" + mtxcmd + "'"

with gzip.open(coofile,'rt',encoding='utf-8') if ( coofile.endswith(".gz") ) else open(coofile,'rt') as fch:
    with gzip.open(f"{args.dge_dir}/{args.bcd}",'rt',encoding='utf-8') if ( args.bcd.endswith(".gz") ) else open(f"{args.dge_dir}/{args.bcd}",'rt') as fbh:
        #with gzip.open(f"{args.dge_dir}/{args.mtx}",'rt',encoding='utf-8') if ( args.mtx.endswith(".gz") ) else open(f"{args.dge_dir}/{args.mtx}",'rt') as fmh:
        proc = sp.Popen(mtxcmd, shell=True, encoding='utf-8', stdout=sp.PIPE)
        with proc.stdout as fmh:
            ## read the header for mtx file
            for line in fmh:
                if ( line[0] == '%' ):
                    pass
                else:
                    #(nftr, nbcd, nlines) = [int(x) for x in line.rstrip().split()]
                    nums = [int(x) for x in line.rstrip().split()]
                    (nftr, nbcd, nlines) = nums[0:3]
                    for i in range(0,len(nums),3): ## make sure that numbers are all consistent between mtx files
                        assert nftr == nums[i]
                        assert nbcd == nums[i+1]
                        assert nlines == nums[i+2]
                        
                    break
            ## read actual contents
            cur_imtx = 0      ## index of barcode in matrix.tsv
            cur_ibcd = 0      ## index of barcode in barcodes.tsv -- must match to cur_imtx
            cur_sbcd = ""     ## string of barcodes in barcodes.tsv
            cur_ctoks = [""]  ## spatial coodinates
            g_ibcd = 0
            g_iline = 0            
            n_miss_bcd = 0
            for i in range(nlines):
                mtoks = [int(x) for x in fmh.readline().rstrip().split()]
                if len(mtxs) == 1:
                    (iftr, ibcd, cnts) = (mtoks[0], mtoks[1], mtoks[2:])
                else:
                    (iftr, ibcd, cnts) = (mtoks[0], mtoks[1], mtoks[2:3])
                    for j in range(1,len(mtxs)):
                        assert iftr == mtoks[3*j+0]
                        assert ibcd == mtoks[3*j+1]
                        cnts.append(mtoks[3*j+2])
                if i == 0:
                    len_cnts = len(cnts)
                    cur_sums = [0] * len_cnts
                assert len_cnts == len(cnts)
                
                if ( random.randint(0,1000000) == 0 ):
                #if ( random.randint(0,1) == 0 ):                    
                #    print(f"Reading {args.mtx} at line {i}: {iftr} {cur_imtx} {cur_ibcd} {cur_sbcd} {cur_ctoks[0]} {ibcd} {cnts}",file=sys.stderr)
                    print(f"Reading {args.mtx} at line {i}: {iftr} {ibcd} {cnts} at barcode {cur_sbcd}, with {n_miss_bcd} missing spatial coordinates",file=sys.stderr)                    
                if ( cur_imtx != ibcd ): ## register new barcode
                    if ( ( cur_imtx > 0 ) and ( cur_ctoks[0] == cur_sbcd ) ):
                        #print(f"{i}\t{cur_imtx}\t{cur_ibcd}\t{cur_sbcd}\n", file=sys.stderr)
                        str_sums = ",".join([str(x) for x in cur_sums])
                        cur_info[3].write(f"{cur_sbcd}\t{cur_info[0]}\t{cur_ibcd}\t{sum(cur_sums)}\t{lane}\t{tile}\t{x}\t{y}\t{str_sums}\n")
                        g_bcd_fh.write(f"{cur_sbcd}\t{g_ibcd}\t{cur_ibcd}\t{sum(cur_sums)}\t{lane}\t{tile}\t{x}\t{y}\t{str_sums}\n")
                        cur_sums = [0] * len_cnts                        
                
                    ## sanity check
                    if ( cur_imtx > ibcd ):
                        raise ValueError(f"The {args.mtx} file is not sorted by barcode index (second column)")
                    ## get the barcode from bcd
                    while( cur_ibcd < ibcd ):
                        cur_ibcd += 1
                        sbcd = fbh.readline().rstrip()
                        if ( cur_sbcd > sbcd ):
                            raise ValueError(f"The {args.bcd} file is not alphanumerically sorted")
                        cur_sbcd = sbcd
                    assert ibcd == cur_ibcd
                    ## get the barcode from coo files
                    while( cur_ctoks[0] < cur_sbcd ):
                        ctoks = fch.readline().rstrip().split()
                        if ( len(ctoks) > 0 ) and ( cur_ctoks[0] > ctoks[0] ):
                            raise ValueError(f"The {coofile} file is not alphanumerically sorted")
                        cur_ctoks = ctoks if ( len(ctoks) > 0 ) else ["ZZZ_EOF_BARCODE"]
                    cur_imtx = ibcd
                    if ( cur_ctoks[0] == cur_sbcd ): ## found barcode with spatial coordinates
                        (lane, tile, x, y) = cur_ctoks[1:]
                        if ( lane not in lane2tile2info ):
                            lane2tile2info[lane] = {}
                        if ( tile not in lane2tile2info[lane] ):
                            ## create the output directory
                            outdir = f"{args.out}/{lane}/{tile}"
                            if not os.path.exists(outdir):
                                os.makedirs(outdir)
                            mtx_fh = open(f"{outdir}/.tmp.matrix.mtx",'wt',encoding='utf-8')
                            bcd_fh = gzip.open(f"{outdir}/barcodes.tsv.gz",'wt',encoding='utf-8')
                            ftr_fh = gzip.open(f"{outdir}/features.tsv.gz",'wt',encoding='utf-8')
                            cur_info = lane2tile2info[lane][tile] = [0, 0, mtx_fh, bcd_fh, ftr_fh] ## ibcd, iline, ...
                        else:
                            cur_info = lane2tile2info[lane][tile]
                        ## write the barcodes
                        cur_info[0] += 1
                        g_ibcd += 1
                    else: ## could not find barcode with spatial coordinates
                        n_miss_bcd += 1
                        if ( n_miss_bcd < args.skip_errors_after ):
                            print(f"Could not find a spatial coodinates for barcode {cur_sbcd}. {n_miss_bcd} missing total so far")
                        elif ( n_miss_bcd == args.skip_errors_after ):
                            print(f"Skipping further warning messages (>{n_miss_bcd}) on missing spatial coordinates..") 
                        continue
                assert cur_imtx == cur_ibcd
                if ( cur_ctoks[0] != cur_sbcd ):
                    continue
                    #raise ValueError(f"{i} {cur_ctoks[0]} {cur_sbcd} {n_miss_bcd}")

                assert cur_ctoks[0] == cur_sbcd

                ## write matrix.mtx
                (lane, tile, x, y) = cur_ctoks[1:]
                str_cnts = " ".join([str(x) for x in cnts])
                cur_info[2].write(f"{iftr} {cur_info[0]} {str_cnts}\n")
                g_mtx_fh.write(f"{iftr} {g_ibcd} {str_cnts}\n")
                cur_info[1] += 1
                g_iline += 1
                if ( iftr not in iftr2cnts ):
                    iftr2cnts[iftr] = [0] * len_cnts
                for i in range(len_cnts):
                    iftr2cnts[iftr][i] += cnts[i]
                    cur_sums[i] += cnts[i]
                if ( lane not in lane2tile2iftr2cnts ):
                    lane2tile2iftr2cnts[lane] = {}
                if ( tile not in lane2tile2iftr2cnts[lane] ):
                    lane2tile2iftr2cnts[lane][tile] = {}
                d = lane2tile2iftr2cnts[lane][tile]
                if ( iftr not in d ):
                    d[iftr] = [0] * len_cnts
                for i in range(len_cnts):
                    d[iftr][i] += cnts[i]
            str_sums = ",".join([str(x) for x in cur_sums])                    
            cur_info[3].write(f"{cur_sbcd}\t{cur_info[0]}\t{cur_ibcd}\t{sum(cur_sums)}\t{lane}\t{tile}\t{x}\t{y}\t{str_sums}\n")
            ## KNOWN BUG : g_ibcd and cur_info[0] is 1 smaller than expected sometimes at the last line
            g_bcd_fh.write(f"{cur_sbcd}\t{g_ibcd}\t{cur_ibcd}\t{sum(cur_sums)}\t{lane}\t{tile}\t{x}\t{y}\t{str_sums}\n")
        proc.wait()

print(f"Writing features.. to {args.dge_dir}/{args.ftr}\n", file=sys.stderr)
## write feature.tsv file
with gzip.open(f"{args.dge_dir}/{args.ftr}",'rt',encoding='utf-8') if ( args.ftr.endswith(".gz") ) else open(f"{args.dge_dir}/{args.ftr}",'rt') as fth:
    iftr = 0
    for line in fth:
        toks = line.rstrip().split('\t')
        iftr += 1
        #print(iftr,file=sys.stderr)
        cnts = iftr2cnts[iftr] if ( iftr in iftr2cnts ) else [0] * len_cnts
        str_cnts = ",".join([str(x) for x in cnts])
        g_ftr_fh.write("\t".join(toks + [str(iftr), str(sum(cnts)), str_cnts])+"\n")
        for lane in lane2tile2iftr2cnts:
            for tile in lane2tile2iftr2cnts[lane]:
                d = lane2tile2iftr2cnts[lane][tile]
                cnts = d[iftr] if ( iftr in d ) else [0] * len_cnts
                str_cnts = ",".join([str(x) for x in cnts])                
                lane2tile2info[lane][tile][4].write("\t".join(toks + [str(iftr), str(sum(cnts)), str_cnts])+ "\n")


for lane in lane2tile2info:
    for tile in lane2tile2info[lane]:
        print(f"Re-writing matrix files for lane {lane} tile {tile}",file=sys.stderr)
        (nbcds, nlines, mtx_fh, bcd_fh, ftr_fh) = lane2tile2info[lane][tile]
        mtx_fh.close()
        bcd_fh.close()
        ftr_fh.close()
        proc = sp.Popen(f"cat - {args.out}/{lane}/{tile}/.tmp.matrix.mtx | {args.pigz} -p {args.ncpus} > {args.out}/{lane}/{tile}/matrix.mtx.gz", shell=True, encoding='utf-8', stdin=sp.PIPE)
        with proc.stdin as wh:
            wh.write(f"%%MatrixMarket matrix coordinate integer general\n%\n{iftr} {nbcds} {nlines}\n")
        proc.wait()
        os.unlink(f"{args.out}/{lane}/{tile}/.tmp.matrix.mtx")
g_mtx_fh.close()
g_bcd_fh.close()
g_ftr_fh.close()

print(f"Writing global matrix files for all lanes and tiles",file=sys.stderr)
#with sp.Popen(f"cat - {args.out}/.tmp.matrix.mtx | {args.pigz} -p {args.ncpus} > {args.out}/matrix.mtx.gz; rm {args.out}/.tmp.matrix.mtx", shell=True, encoding='utf-8', stdin=sp.PIPE).stdin as wh:
proc = sp.Popen(f"cat - {args.out}/.tmp.matrix.mtx | {args.pigz} -p {args.ncpus} > {args.out}/matrix.mtx.gz", shell=True, encoding='utf-8', stdin=sp.PIPE)
with proc.stdin as wh:
    wh.write(f"%%MatrixMarket matrix coordinate integer general\n%\n{iftr} {g_ibcd} {g_iline}\n")
proc.wait()
os.unlink(f"{args.out}/.tmp.matrix.mtx")
