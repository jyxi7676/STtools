#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
DGEdir=args[1]
spatial=args[2]
tiles=args[3]
nrow=as.numeric(args[4])
ncol=as.numeric(args[5])
sidesize=as.numeric(args[6])
outpath=args[7]
collapsePath=args[8]
layout=args[9]
nMax=args[10]
clustering=args[11]
tiles=(unlist(strsplit(tiles,',')))
print(tiles)









####################################################################################################3
####################################################################################################3
#' This function collapse tiles into small grids and create Seurat objects
#' @param tile_df Dataframe including coordinates
#' @param i integer of tile number
#' @param binx binning size for x coordinates
#' @param biny binning size for y coordinates
#' @param m_tile sparse count matrix
#' @export
collapseTiles=function(tile_df,i,binx,biny,m_tile)
{

  #print('start collaping tile')
  #print(i)
  tile_df_i=tile_df[tile_df$tile_miseq==i,]
  if(dim(tile_df_i)[1]==0)
  {
    return(NULL)
  }
  miny = min(tile_df_i$y_miseq)
  maxy = max(tile_df_i$y_miseq)
  minx = min(tile_df_i$x_miseq)
  maxx = max(tile_df_i$x_miseq)
  xlim= c(min(tile_df_i$x_miseq),max(tile_df_i$x_miseq))
  ylim = c(min(tile_df_i$y_miseq),max(tile_df_i$y_miseq))
  xlim = c(min(tile_df_sub_wind$x_miseq),max(tile_df_sub_wind$x_miseq))
  ylim = c(min(tile_df_sub_wind$y_miseq),max(tile_df_sub_wind$y_miseq))
  print('xlim')
  print(xlim)
  print(ylim)
  print(binx)
                                        #print(xlim[0])                                                                                                                                                     
  gapx=xlim[2]-xlim[1]
  gapy=ylim[2]-ylim[1]
  if (gapx<300)
   {
       xlim[1]=xlim[1]-gapx/2
       xlim[2]=xlim[2]+gapx/2

   }
   if (gapy<300)
   {
       ylim[1]=ylim[1]-gapy/2
       ylim[2]=ylim[2]+gapy/2

   }

 grd = make.grid(tile_df_i$x_miseq,tile_df_i$y_miseq,tile_df_i$UMI, binx,biny, xlim, ylim)
 grd=t(grd)


 fn1=paste0('Temp_CollapsedHDMIsIndLength','.csv')
 fn2=paste0('Temp_CollapsedHDMIsInd','.txt')

  #fn1=paste0('Temp_CollapsedHDMIsIndLength',tile,'_',j,t,'.csv')
  #fn2=paste0('Temp_CollapsedHDMIsInd',tile,'_',j,t,'.txt')
  
  if (file.exists(fn1)) {
    #Delete file if it exists
    file.remove(fn1)
  }

  if (file.exists(fn2)) {
    #Delete file if it exists
    file.remove(fn2)
  }
  #grd_re = make.grid(tile_df_sub_wind$x_miseq,tile_df_sub_wind$y_miseq,tile_df_sub_wind$tileHDMIind, binx, biny, xlim2, ylim2,function(x) {write(length(x), file=fn1,append = T)})
  #grd_re = make.grid(tile_df_sub_wind$x_miseq,tile_df_sub_wind$y_miseq,tile_df_sub_wind$tileHDMIind, binx, biny, xlim2, ylim2,function(x) {cat(x,file=fn2,append=TRUE,sep='\n')})
  test.env <- new.env()
  test.env$len =c()
  test.env$no = c()
  grd_re = make.grid(tile_df_i$x_miseq,tile_df_i$y_miseq,tile_df_i$tileHDMIind, binx, biny, xlim, ylim,function(x) {test.env$len=c(test.env$len,length(x))})
  grd_re = make.grid(tile_df_i$x_miseq,tile_df_i$y_miseq,tile_df_i$tileHDMIind, binx, biny, xlim, ylim,function(x) {test.env$no=c(test.env$no,x)})  
 
  collapseLen = test.env$len
  collapseInd = test.env$no
  collapseLen = cbind(collapseLen,cumsum(collapseLen))
  colnames(collapseLen) =c("len","end")
  collapseLen=as.data.frame(collapseLen)
  interv = c(0,collapseLen$end)

  #sourceCpp(collapsePath)
  print('Start Simple Square Gridding!')
  tic();collapseM = collapse(m_tile,collapseInd,interv);toc()
  rownames(collapseM) = rownames(m_tile)
  colnames(collapseM) = paste0("Collapse_",i,'_',1:(length(interv)-1))
  sparse.gbm <- Matrix(collapseM , sparse = T )
  #writeMM(obj = sparse.gbm, file=paste0('tile',i,"collapsedMatrix.mtx"))
  #write.csv(rownames(collapseM),paste0('tile',i,'collapsedGenes.csv'))
  #write.csv(colnames(collapseM),paste0('tile',i,'collapsedBarcodes.csv'))

  grd=t(grd)
  pos = which(!is.na(grd), TRUE)
  pos_coor = t(sapply(1:(dim(pos)[1]),function(x) {c(as.numeric(rownames(grd)[as.numeric(pos[x,1])]),as.numeric(colnames(grd)[as.numeric(pos[x,2])]) )})) #here the grd is transposed
  colnames(pos_coor) =c("X","Y")
  coord.df = data.frame("Y"=pos_coor[,2], "X"=pos_coor[,1],"tile"=i, stringsAsFactors=FALSE)
  #write.csv(coord.df,paste0('tile',i,'coor_df_stratgy1.csv'))


  obj1 = CreateSeuratObject(counts=collapseM,assay='Spatial')
  #obj1$status = "Original"
  obj1@meta.data$tile = coord.df$tile
  obj1@meta.data$X = coord.df$X
  obj1@meta.data$Y = coord.df$Y
  obj1@meta.data$tile = i
  return(obj1)

}








mergeSeuratObj=function(seurat_object_list)
{
  seurat_object_list = unlist(seurat_object_list)

  for (i in names(seurat_object_list))
  {
    seurat_object_list[[i]] = RenameCells(seurat_object_list[[i]],
                                          add.cell.id = i)
  }
  merged_combined = suppressWarnings(expr=reduce(seurat_object_list,
                                                 merge,
                                                 do.normalize = FALSE))
  return(merged_combined)
}










####################################################################################################3
#' This function grids the data with user-defined binning size and collapsed HDMIs within each grid
#' @param seqscope1st Data flatform,
#' @param DGEdir folder that stores barcodes.tsv, features.tsv and matrix.mtx
#' @param spatial txt file stores spatial informaiton with four columns: 'HDMI','tile_miseq','x_miseq','y_miseq'
#' @param tiles a vector of tiles that the user is interested in collapsing
#' @param nrow an integer of how many rows to organize the tiles
#' @param ncol an integer of how many cols to organize the tiles
#' @param sidesize side size of the square grid 300 represents 10um
#' @param outpath path to store the output, you need to make sure the path exists before running the function
#' @import Seurat
#' @import Matrix
#' @import ggplot2
#' @export

getSimpleGrid = function(DGEdir,spatial,tiles,nrow,ncol,sidesize,outpath,collapsePath,layout,nMax,clustering)
{
  if (missing(layout))
  {
    layout='MiSeq'
  }

  if(missing(sidesize))
  {
    sidesize=300
  }

  if (!dir.exists(DGEdir)){
    stop("DGEdir does not exist")
  }

  if (!dir.exists(outpath)){
    stop("outpath does not exist")
  }

  if (!file.exists(spatial))
  {
    stop("spatial data does not exist")

  }


  #install required packages
  packages = c("R.utils","data.table","tidyverse","Matrix", "tictoc", "ggplot2", "ggsci","Seurat","mapplots","rlist","cowplot","dplyr","Rcpp")
  ## add more packages to load if needed
  ## Now load or install&load all
  package.check <- lapply(
    packages,
    FUN = function(x) {
      if (!require(x, character.only = TRUE)) {
        install.packages(x, repos="https://cran.rstudio.com", dependencies = TRUE)
        library(x, character.only = TRUE)
      }
    }
  )

  #source collapsed.cpp
  sourceCpp(collapsePath)
  setwd(DGEdir)
  biny = binx = sidesize

  #read files
  print("Read files")
  #bc = read.table("barcodes.tsv",header=F)$V1
  bc=fread('barcodes.tsv',header=F)$V1
  features = read.table('features.tsv',header=F)$V2
  m = readMM('matrix.mtx')
  if(any(c(length(features),length(bc)) != dim(m)))
  {
    stop('Dimension of matrix.mtx does not match with features or barcodes')
  }
  rownames(m) = features
  colnames(m) = bc

  if (nMax!='None')
  {
    m = m[,colSums(m)<=nMax&colSums(m)>0]  #remove outliers
  }


  #get spatial info
  miseq_pos = read.table(spatial)
 # print(head(miseq_pos))
  colnames(miseq_pos) = c('HDMI','lane_miseq','tile_miseq','x_miseq','y_miseq')
  miseq_pos$tile_miseq=paste(miseq_pos$lane_miseq,miseq_pos$tile_miseq,sep="_")
  tiles=(unlist(strsplit(tiles,',')))
  tiles_uniq=unique(miseq_pos$tile_miseq)
  print(tiles)
  print(tiles_uniq)
  if (layout=='HiSeq')
  { 
    lane1_tile=c(paste('1',1201:1216,sep='_'),paste('1',2201:2216,sep='_'),paste('1',1101:1116,sep='_'),paste('1',2101:2116,sep='_'))
    lane1_rtile=sapply(1:length(lane1_tile),function(i) {strsplit(lane1_tile[i],'_')[[1]][2]})
    lane1_layout=data.frame('ROW'=c(rep(1,32),rep(2,32)),'COL'=c(seq(1,32,1),seq(1,32,1)),'LANE'=1,'OTILE'=lane1_rtile)

    lane2_tile=c(paste('2',1201:1216,sep='_'),paste('2',2201:2216,sep='_'),paste('2',1101:1116,sep='_'),paste('2',2101:2116,sep='_'))
    lane2_rtile=sapply(1:length(lane2_tile),function(i) {strsplit(lane2_tile[i],'_')[[1]][2]})
    lane2_layout=data.frame('ROW'=c(rep(3,32),rep(4,32)),'COL'=c(seq(1,32,1),seq(1,32,1)),'LANE'=2,'OTILE'=lane2_rtile)

    layout=rbind(lane1_layout,lane2_layout)
    #lane12=c(paste('1',1101:1116,sep='_'),paste('1',1201:1216,sep='_'),paste('1',2101:2116,sep='_'),paste('1',2201:2216,sep='_'),paste('2',1101:1116,sep='_'),paste('2',1201:1216,sep='_'),paste('2',2101:2116,sep='_'),paste('2',2201:2216,sep='_'))
    #nrow=4
    #ncol=ceiling(length(lane12)/nrow)
  } else if  (layout == 'MiSeq')
  {
    print('In miseq')
    lane1_tile=tiles_uniq[order(tiles_uniq,decreasing=F)]
    lane1_rtile=sapply(1:length(lane1_tile),function(i) {strsplit(lane1_tile[i],'_')[[1]][2]})
    nrow=1
    ncol=length(lane1_rtile)
    layout=data.frame('ROW'=nrow,'COL'=1:ncol,'LANE'=1,'OTILE'=lane1_rtile)

  } else 
  {
    #check if path exists
   layout=read.csv(layout,row.names=1)
   colnames(layout)[4]='OTILE'
  

   print('Customized')
  }

 layout$tile=paste(layout$LANE,layout$OTILE,sep='_')
 #lane12=unique(layout$tile)
 #lane12=lane12[order(lane12)]

  if (tiles=='All')
  {
    print('cond2')
    tiles=tiles_uniq
    tiles=tiles[order(tiles)]
  } else
  { print('cond3')
    tiles = intersect(tiles,tiles_uniq)
    tiles=tiles[order(tiles)]
    print('tiles after merge')
    print(tiles)
    miseq_pos=miseq_pos[miseq_pos$tile_miseq %in% tiles,]
    if (dim(miseq_pos)[1]==0)
    {
      stop("No tiles found")
    }
  }


  setwd(outpath)
  print('merge')
  df = data.frame("HDMI" =colnames(m),"HDMIind" = 1:(dim(m)[2]))
  tile_df = merge(miseq_pos,df,by = "HDMI")
  m_tile = m[,tile_df$HDMIind]
  tile_df$UMI = colSums(m_tile)
  tile_df$tileHDMIind= match(tile_df$HDMI,colnames(m_tile))
  tile_df=tile_df[tile_df$tile_miseq %in% layout$tile,]
  write.csv(tile_df,'tile_df_raw.csv')

  print('Start collapsing')
  obj=sapply(tiles,collapseTiles, tile_df=tile_df,binx=binx,biny=biny,m_tile=m_tile)
  obj=mergeSeuratObj(obj)
  print('OBJ')
  print(obj)

  #super tile
  tile_df = obj@meta.data
  tile_df$name=colnames(obj)
  #tile_df$tile=as.factor(tile_df$tile)                                       #aggregate all tils by expanding coord
  
  #ordering=data.frame('tile'=lane12,'aggrInd'=order(as.factor(lane12))-1)
 # tile_df=merge(tile_df,ordering,by='tile')

  addson_hori = max(tile_df$Y)
  addson_verti = max(tile_df$X)
  nrow = max(layout$ROW)
  ncol = max(layout$COL)
  tile_df = merge(tile_df,layout,by='tile')
  print(head(tile_df))
  print(head(obj@meta.data))
  addson_hori = max(tile_df$Y)
  addson_verti = max(tile_df$X)
  tile_df$y_miseq_expand = tile_df$Y +  addson_hori*(tile_df$COL-1)
  #tile_df$x_miseq_expand = tile_df$X +  addson_verti*(tile_df$ROW-1)
  tile_df$x_miseq_expand = tile_df$X +  addson_verti*(nrow-tile_df$ROW)


  # if(layout=='HiSeq')
  # {
  #     tile_df$x_miseq_expand =   tile_df$X +  addson_verti*(floor((tile_df$aggrInd2/ncol)))

  # } else if (layout=='MiSeq')
  # {
  #     tile_df$x_miseq_expand =   tile_df$X 

  # } else
  # {
  #     #print('layout')
  #     #layout = read.csv(layout,row.names=1)
  #     nrow = max(layout$ROW)
  #     ncol = max(layout$COL)
  #    # layout$tile=  paste(layout[,3],layout[,4],sep="_")
  #     #layout$tile=paste(layout$LANE,layout$TILE,sep="_")
  #    #layout$ROW_COL=paste(layout$ROW,layout$COL,sep="_")
  #     #layout$aggrInd =  as.numeric(factor(layout$ROW_COL))-1
  #     print(head(layout))

  #     print(head(tile_df))
  #     #tile_df$aggrInd2 = layout$aggrInd

  #     tile_df = merge(tile_df,layout,by='tile')
  #     addson_hori = max(tile_df$Y)
  #     addson_verti = max(tile_df$X)
  #     tile_df$y_miseq_expand = tile_df$Y +  addson_hori*(tile_df$COL-1)
  #     tile_df$x_miseq_expand = tile_df$X +  addson_verti*(tile_df$ROW-1)

  # }
  

  # #rownames(tile_df)=colnames(obj)

  # #print(head(tile_df))
  tile_df=tile_df[,c('orig.ident','name','LANE','OTILE','tile','nCount_Spatial','nFeature_Spatial','X','Y','x_miseq_expand','y_miseq_expand','ROW','COL')]
  colnames(tile_df)=c('orig.ident','collapseID','lane','tile','lane_tile','nCount_Spatial','nFeature_Spatial','X','Y','X_expand','Y_expand','row','col')
  obj@meta.data=tile_df
  rownames(obj@meta.data)=colnames(obj)
  #obj@meta.data$X_expand = obj@meta.data$x_miseq_expand
  #obj@meta.data$Y_expand = obj@meta.data$y_miseq_expand
  print(obj)
  print(head(obj@meta.data))
  obj@images$image = new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = obj@meta.data[,c('Y_expand','X_expand')]
  )

  print(obj)
  m=(obj@assays$Spatial@counts)
  gene=rownames(obj)
  bc=colnames(obj)
  writeMM(obj = m, file="collapsedMatrix.mtx")
  write.csv(gene,'collapsedGenes.csv')
  write.csv(bc,'collapsedBarcodes.csv')
  saveRDS(obj,'SimpleSquareGrids.RDS')
  junk = dir(path=outpath,  pattern="Temp")
  file.remove(junk)
  junk = dir(path=outpath,  pattern="^tile.*\\Barcodes.csv")
  file.remove(junk)
  junk = dir(path=outpath,  pattern="^tile.*\\Genes.csv")
  file.remove(junk)
  junk = dir(path=outpath,  pattern="^tile.*\\Matrix.mtx")
  file.remove(junk)
                                        #   print('feature eplot')
  #png("nFeatureplot.png",width=7*2,height=6,res=300,units='in')
  #p1=VlnPlot(obj, features = "nFeature_Spatial", pt.size = 0) + NoLegend()
  #p2=VlnPlot(obj, features = "nFeature_Spatial", pt.size = 0,log=T) + NoLegend()
  #print(plot_grid(p1,p2,ncol=2))
  #dev.off()

#   if(nMax !='None')
#   {
#
#     geneCount1=median(obj@meta.data$nFeature_Spatial)  #?automatic ??????
#     print('geneCount')
#     print(geneCount1)
#
#   }
#   else
#       geneCount1=0
#
#
#   print(clustering)
  if(clustering)
  {

    geneCount1=median(obj@meta.data$nFeature_Spatial)  #?automatic ??????
   # print('geneCount')
   # print(geneCount1)
    obj_simple = subset(obj, subset = nFeature_Spatial>geneCount1)
    obj_simple = SCTransform(obj_simple, assay = "Spatial")
    obj_simple = RunPCA(obj_simple)
    obj_simple = RunUMAP(obj_simple, dims = 1:30)
    obj_simple = FindNeighbors(obj_simple, dims = 1:30)
    obj_simple = FindClusters(obj_simple)
    saveRDS(obj_simple,'SimpleSquareGridsWithClustering.RDS')
    meta=obj_simple@meta.data
    meta=cbind(meta,obj_simple@reductions$umap@cell.embeddings)
    print(nrow)
    print(ncol)
    png('SpatialDimPlot.png',width=ncol*7,height=nrow*6,res=300,units='in')
    print(ggplot(meta,aes(X_expand,Y_expand,color=seurat_clusters))+geom_point(size=1,alpha=1))
    dev.off()

    png('UMAP.png',width=ncol*7,height=6,res=300,units='in')
    print(DimPlot(obj, reduction = "umap", label = TRUE))
    dev.off()
  }

  print('Done!')

}



getSimpleGrid(DGEdir,spatial,tiles,nrow,ncol,sidesize,outpath,collapsePath,layout,nMax,clustering)






