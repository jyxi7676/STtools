
def import_or_install(package):
    try:                                                                                                                                                                                 
        __import__(package)
    except ImportError:
        #print('install')                                                                                                                                                                                 
        pip.main(['install', package])

#modules                                                                                                                                                                                                
import pip
modules = ["sys","os","gzip","math","ranodm","numpy","pandas","argparse","subprocess","PIL"]
#map(import_or_install,modules)                                                                                                                                                                           
md=[import_or_install(i) for i in modules]

import os, sys, gzip, math, random
import pandas as pd
import numpy as np
import argparse
import subprocess as sp
from PIL import Image

parser = argparse.ArgumentParser(description="Visualize tiles for SeqScope data")
parser.add_argument("-d", "--in-dir", type=str, required=True, help="Input (STTools output) directory")
parser.add_argument("-t", "--tile", type=str, help="A single tile to draw (in [lane_tile] format)")
parser.add_argument("-l", "--layout", type=str, help="Layout file of tiles to draw [lane] [tile] [row] [col] format in each line")
parser.add_argument("-s", "--scale", type=float, default=20.0, help="Scale each color to have the same mean intensity")
parser.add_argument("-i", "--inv-weight", default=False, action='store_true', help="Weight each gene inverse")
parser.add_argument("-m", "--min-tpm", default=1000, type=float, help="Minimum TPM value for inverse weighting (effective only with --inv-weight)")
parser.add_argument("-r", "--red", type=str, required=True, help="Comma-separate list of genes (colon with weights) for red color")
parser.add_argument("-g", "--green", type=str, required=True, help="Comma-separate list of genes (colon with weights) for green color")
parser.add_argument("-b", "--blue", type=str, required=True, help="Comma-separate list of genes (colon with weights) for blue color")
parser.add_argument("--res", type=int, default=80, help="Resolution of pixel (how to bin each pixel) - default: 80 (um2 per pixel)")
parser.add_argument("-o", "--out", type=str, required=True, help="Output file prefix")

args = parser.parse_args()

def update_xyrange(xyr, x, y):
    #print(f"{x} {y} {xyr}",file=sys.stderr)
    if xyr[0][0] > x:
        xyr[0][0] = x
    if xyr[0][1] < x:
        xyr[0][1] = x
    if xyr[1][0] > y:
        xyr[1][0] = y
    if xyr[1][1] < y:
        xyr[1][1] = y

def parse_gene_weights(s):
    ## expected arguments [GENE_1],[GENE_2],...,[GENE_N], or [GENE_1]:[W1],[GENE2]:[W2],... default weight is 1
    ## Use "_sum" to sum all genes, to count all pixels
    toks = s.split(',')
    d = {}
    for tok in toks:
        gtoks = tok.split(':')        
        if ( len(gtoks) == 1 ):
            d[gtoks[0]] = [1.0, 0]
        elif ( len(gtoks) == 2 ):
            d[gtoks[0]] = [float(gtoks[1]), 0]
        elif ( len(gtoks) == 3 ):
            d[gtoks[0]] = [float(gtoks[1]), int(gtoks[2])]
        else:
            raise ValueError(f"Cannot parse {tok} in {s}")
    return d

## parse --red, --blue, --green arguments
rgbGWs = [parse_gene_weights(args.red), parse_gene_weights(args.green), parse_gene_weights(args.blue)]

#print(rgbGWs, file=sys.stderr)

if ( args.layout is not None ):
    if ( args.tile is not None ):
        raise ValueError("Cannot use -t/--tile and -l/--layout options together")
    df = pd.read_csv(args.layout,sep="\t")
    nrows = max(df['row'])
    ncols = max(df['col'])
    lanes = []
    tiles = []
    for i in range(nrows):
        lanes.append( [None] * ncols )
        tiles.append( [None] * ncols )
    for index, row in df.iterrows():
        i = int(row['row'])-1
        j = int(row['col'])-1
        lanes[i][j] = str(row['lane'])
        tiles[i][j] = str(row['tile'])
    ## process later
elif ( args.tile is not None ):
    (lane, tile) = args.tile.split('_')
    ## matrix 1x1 matrix (list of list)
    lanes = [ [lane] ]
    tiles = [ [tile] ]
    nrows = 1
    ncols = 1
else:
    raise ValueError(f"ERROR: Either -l/--layout or -t/--tile option is required (but not both)")

## read each tile and construct bins of sums and npixs
xyrange = [[sys.maxsize, 0], [sys.maxsize, 0]]
bin2cnts = {}

## calculate per gene weights
geneTPMs = None
if args.inv_weight:
    geneTPMs = []
    for r in range(nrows):
        for c in range(ncols):
            with gzip.open(f"{args.in_dir}/{lanes[r][c]}/{tiles[r][c]}/features.tsv.gz",'rt',encoding='utf-8') as fh:
                idx = 0
                for line in fh:
                    toks = line.rstrip().split('\t')
                    cnts = [int(x) for x in toks[-1].split(',')]
                    if r == 0 and c == 0:
                        geneTPMs.append(cnts)
                    else:
                        for i in range(len(cnts)):
                            geneTPMs[idx][i] += cnts[i]
                    idx += 1
    geneSum = 0
    for i in range(len(geneTPMs)):
        for j in range(len(geneTPMs[i])):
            geneSum += geneTPMs[i][j]
    
    for i in range(len(geneTPMs)):
        for j in range(len(geneTPMs[i])):
            geneTPMs[i][j] = geneTPMs[i][j] / geneSum * 1e6
            if geneTPMs[j][j] < args.min_tpm:
                geneTPMs[i][j] = args.min_tpm

for r in range(nrows):
    for c in range(ncols):
        print(f"Processing ({r},{c}) : {lanes[r][c]}_{tiles[r][c]}",file=sys.stderr)
        ## read features.tsv
        colDict = {} ## iftr -> RGB weights
        with gzip.open(f"{args.in_dir}/{lanes[r][c]}/{tiles[r][c]}/features.tsv.gz",'rt',encoding='utf-8') as fh:
            for line in fh:
                toks = line.rstrip().split('\t')
                gidx = int(toks[3])-1
                rgbw = [0, 0, 0]
                rgbi = [0, 0, 0]
                for i in range(3):
                    g = "_all" if ( "_all" in rgbGWs[i] ) else toks[1]
                    if ( g in rgbGWs[i] ):
                        (weight, idx) = rgbGWs[i][g]
                        #print(f"{idx} {i}",file=sys.stderr)
                        rgbw[i] = weight / (1 if ( geneTPMs is None ) else geneTPMs[gidx][idx])
                        rgbi[i] = idx
                if ( rgbw[0] + rgbw[1] + rgbw[2] > 0 ):
                    colDict[toks[3]] = [rgbw, rgbi] ## make indices

        #print(colDict)
                        
        with gzip.open(f"{args.in_dir}/{lanes[r][c]}/{tiles[r][c]}/barcodes.tsv.gz",'rt',encoding='utf-8') as bh:
            with gzip.open(f"{args.in_dir}/{lanes[r][c]}/{tiles[r][c]}/matrix.mtx.gz",'rt',encoding='utf-8') as mh:
                for line in mh:  ## skip headers
                    if ( line[0] != '%' ):
                        mtoks = line.split()
                        nlines = int(mtoks[2])
                        break
                btoks = bh.readline().split('\t')
                update_xyrange(xyrange, int(btoks[6]), int(btoks[7]))
                for i in range(nlines):
                    mtoks = mh.readline().split()
                    while ( int(btoks[1]) < int(mtoks[1]) ):
                        btoks = bh.readline().split('\t')
                        update_xyrange(xyrange, int(btoks[6]), int(btoks[7]))                        
                    if ( btoks[1] != mtoks[1] ):
                        raise ValueError(f"Cannot find matched barcode between barcodes and matrix. {btoks[1]}!={mtoks[1]}")
                    if ( mtoks[0] in colDict ): ## gene has to be counted
                        (rgbW, rgbI) = colDict[mtoks[0]]
                        (lane, tile, x, y) = (btoks[4],btoks[5],int(btoks[6]),int(btoks[7]))
                        cnts = [int(z) for z in btoks[8].split(',')]
                        xbin = x // args.res
                        ybin = y // args.res
                        key = ":".join([str(r),str(c),str(xbin),str(ybin)])                    
                        if key not in bin2cnts:
                            bin2cnts[key] = [0, 0, 0]
                        for i in range(3):
                            bin2cnts[key][i] += (rgbW[i] * cnts[rgbI[i]])
                            #if ( cnts[1] > 0 ):
                            #    print(f"{key} {i} {bin2cnts[key][i]} {rgbWI[0][i]} {rgbWI[1][i]} {cnts}",file=sys.stderr)

## calculate the range of bins
## print(xyrange, file=sys.stderr)
(xbin_min, xbin_max, ybin_min, ybin_max) = (xyrange[0][0]//args.res, xyrange[0][1]//args.res, xyrange[1][0]//args.res, xyrange[1][1]//args.res)

## construct a sparse matrix
tile_h  = xbin_max - xbin_min + 1
total_h = tile_h * nrows
tile_w  = ybin_max - ybin_min + 1
total_w = tile_w * ncols

print(f"Calculate mean intensity of each color",file=sys.stderr)
sums = [0,0,0]
for (key, value) in bin2cnts.items():
#   print(f"{key} {value}",file=sys.stderr)
    sums[0] += value[0]
    sums[1] += value[1]
    sums[2] += value[2]
sums[0] = sums[0] / total_h / total_w + 1e-10
sums[1] = sums[1] / total_h / total_w + 1e-10
sums[2] = sums[2] / total_h / total_w + 1e-10
scale_factors = (args.scale / sums[0], args.scale / sums[1], args.scale / sums[2])

print(f"Scale_factors = {scale_factors}",file=sys.stderr)

print(f"Constructing an image of {total_h} x {total_w}",file=sys.stderr)             

data = np.zeros( (total_h, total_w, 3), dtype=np.uint8 )
for (key, value) in bin2cnts.items():
    (row, col, xbin, ybin) = [int(x) for x in key.split(':')]
    x_ext = row * tile_h + (tile_h - xbin + xbin_min - 1)
    y_ext = col * tile_w + (ybin - ybin_min)
    #print(f"{key} {tile_h} {x_ext} {y_ext}",file=sys.stderr)
    red = value[0] * scale_factors[0]
    grn = value[1] * scale_factors[1]
    blu = value[2] * scale_factors[2]
    data[x_ext,y_ext,0] = int(red) if ( red < 256 ) else 255
    data[x_ext,y_ext,1] = int(grn) if ( grn < 256 ) else 255
    data[x_ext,y_ext,2] = int(blu) if ( blu < 256 ) else 255
    #print(f"{x_ext}\t{y_ext}\t{red}\t{grn}",file=sys.stderr)
for r in range(1,nrows):
    data[tile_h * r, :, :] = 255
for c in range(1,ncols):
    data[:, tile_w * c, :] = 255

print(f"Converting the matrix to an image",file=sys.stderr)    
img = Image.fromarray(data)

print(f"Saving the image {args.out}.png",file=sys.stderr)    
img.save(f"{args.out}.png")
