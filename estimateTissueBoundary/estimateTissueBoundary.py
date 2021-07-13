import numpy as np
import pandas as pd
import pandas as pd
import os
import mpl_scatter_density # adds projection='scatter_density'
from matplotlib.colors import LinearSegmentedColormap
import pylab as plt
import sys

def readFiles(pos1stSeq,hdmi2ndSeq):
    """A simple function to read SeqScope files"""

    if os.path.isfile(pos1stSeq)==False:
        raise NameError('Not valide pos1stSeq file')
    if os.path.isfile(hdmi2ndSeq)==False:
        raise NameError('Not valide hdmi2ndSeq file')

    miseq_pos =  pd.read_csv(pos1stSeq,delim_whitespace=True, header=None)
    miseq_pos.columns = ['HDMI','lane','tile','x','y']
    hiseq = pd.read_csv(hdmi2ndSeq,delim_whitespace=True,header=None)
    hiseq.columns = ['HDMI']
    merge_df = pd.merge(hiseq,miseq_pos, on='HDMI',how="inner")
    return merge_df

def using_mpl_scatter_density(fig, x, y,maxScale):
    white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
        (0, '#ffffff'),
        (1e-20, '#440053'),
        (0.5, '#2a788e'),
        (1, '#fde624')], N=256)
    if maxScale is None:
        ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
        density = ax.scatter_density(x, y, cmap=white_viridis)
    else:
        maxScale=float(maxScale)
        ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
        density = ax.scatter_density(x, y, cmap=white_viridis,vmax =maxScale)
    fig.colorbar(density, label='Number of points per pixel')


def estimateTissueBoundary(pos1stSeq,hdmi2ndSeq,maxScale,outpath):
    """Generate plot to visualize tissue boundary"""

    if os.path.isdir(outpath)==True:
        os.chdir(outpath)
    else:
        raise NameError('Not valide output path')

    if os.path.isfile(pos1stSeq)==False:
        raise NameError('Not valide pos1stSeq file')
    if os.path.isfile(hdmi2ndSeq)==False:
        raise NameError('Not valide hdmi2ndSeq file')


    merge_df=readFiles(pos1stSeq,hdmi2ndSeq)
    merge_df["lane_tile"]=merge_df["lane"].astype(str)+"_"+merge_df["tile"].astype(str)
    print(merge_df.head())
    tiles_uniq = merge_df['lane_tile'].unique()
    
    tiles_uniq = list(tiles_uniq)
    print('maxscale')
    print(maxScale)
    if maxScale == '0':
        print('TRUE 0')
        maxScale=None
          
    for i in tiles_uniq:
        x = merge_df[merge_df.lane_tile.eq(i)]
        fig = plt.figure()
        using_mpl_scatter_density(fig, x['y'], x['x'],maxScale)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.savefig('lane_tile'+str(i)+'.png',dpi=1000)
        plt.show()

        # if maxScale == '0':
        #    print('TRUE 0')
         #   maxScale=None
          #  print('use none')
           # using_mpl_scatter_density(fig, x['y'], x['x'],maxScale)
           # print('ok')
#            plt.gca().set_aspect('equal', adjustable='box')
        #plt.savefig('tile'+str(i)+'.png',dpi=1000)
 #           plt.savefig('lane_tile'+str(i)+'.png',dpi=1000)
  #          plt.show()
   #     else:
    #        print('FALSE 0')
     #       using_mpl_scatter_density(fig, x['y'], x['x'],float(maxScale))
            
      #      plt.gca().set_aspect('equal', adjustable='box')
        #plt.savefig('tile'+str(i)+'.png',dpi=1000)
       #     plt.savefig('lane_tile'+str(i)+'.png',dpi=1000)
        #    plt.show()


#takes user-input arguments
pos1stSeq=sys.argv[1] 
hdmi2ndSeq=sys.argv[2] 
maxScale=sys.argv[3] 
outpath=sys.argv[4] 

estimateTissueBoundary(pos1stSeq,hdmi2ndSeq,maxScale,outpath)
