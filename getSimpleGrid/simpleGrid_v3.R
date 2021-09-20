#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
seqscope1st=args[1]
DGEdir=args[2]
spatial=args[3]
tiles=args[4]
nrow=as.numeric(args[5])
ncol=as.numeric(args[6])
sidesize=as.numeric(args[7])
outpath=args[8]
collapsePath=args[9]
layout=args[10]
order=args[11]
nMax=args[12]
clustering=args[13]
print('nMAX')
print(nMax)
print(is.null(nMax))
print(nMax=="None")
print('tiles')
print(tiles)
print(clustering)
#tiles=as.numeric(unlist(strsplit(tiles,',')))
tiles=(unlist(strsplit(tiles,',')))



if(length(tiles)>1)
{
    print('two')
    nrow=2
    ncol=ceiling(length(tiles)/nrow)
}
if(length(tiles)==1)
{
    nrow=1
    ncol=1
}
print('tiles')
print(tiles)
print(order)
print(order=='top')
print(length(tiles))
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

  #assign('var', 100, envir=test.env)
# or simply
  test.env$len =c()
  test.env$no = c()
  grd_re = make.grid(tile_df_i$x_miseq,tile_df_i$y_miseq,tile_df_i$tileHDMIind, binx, biny, xlim, ylim,function(x) {test.env$len=c(test.env$len,length(x))})
  grd_re = make.grid(tile_df_i$x_miseq,tile_df_i$y_miseq,tile_df_i$tileHDMIind, binx, biny, xlim, ylim,function(x) {test.env$no=c(test.env$no,x)})  
  #grd_re = make.grid(tile_df_sub_wind$x_miseq,tile_df_sub_wind$y_miseq,tile_df_sub_wind$tileHDMIind, binx, biny, xlim2, ylim2,function(x) {write(length(x), file=fn1,append = T);cat(x,file=fn2,append=TRUE,sep='\n')})
  #draw the grids centers:
  # write.csv(as.numeric(colnames(grd)),paste0('grd_col',j,t,'.csv'))
  #  write.csv(as.numeric(rownames(grd)),paste0('grd_row',j,t,'.csv'))
  collapseLen = test.env$len
  collapseInd = test.env$no

  # collapseLen = cbind(collapseLen,cumsum(collapseLen$V1))
  # colnames(collapseLen) =c("len","end")
  # interv = c(0,collapseLen$end)


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
  writeMM(obj = sparse.gbm, file=paste0('tile',i,"collapsedMatrix.mtx"))
  write.csv(rownames(collapseM),paste0('tile',i,'collapsedGenes.csv'))
  write.csv(colnames(collapseM),paste0('tile',i,'collapsedBarcodes.csv'))

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

getSimpleGrid = function(seqscope1st,DGEdir,spatial,tiles,nrow,ncol,sidesize,outpath,collapsePath,layout,order,nMax,clustering)
{
  print('nMAX')


  print('inside')                                        #print(tiles)
#  if(missing(nrow)| missing(ncol))
 # {
  #    nrow=2
   #   ncol=ceiling(length(tiles)/nrow)
 # }
  print('seqscope1st')
  print(seqscope1st)
  tile_ind=tiles
  if(missing(seqscope1st))
  {
    seqscope1st="HiSeq"
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


  #install required packages
  packages = c("tidyverse","Matrix", "tictoc", "ggplot2", "ggsci","Seurat","mapplots","rlist","cowplot","dplyr","Rcpp")
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
  bc = read.table("barcodes.tsv",header=F)$V1
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
    m = m[,colSums(m)<=100&colSums(m)>0]  #remove outliers
  }


  #m = m[,colSums(m)<=100&colSums(m)>0]  #remove outliers
  # print(dim(m))

  #get spatial info
  miseq_pos = read.table(spatial)
 # print(head(miseq_pos))
  colnames(miseq_pos) = c('HDMI','lane_miseq','tile_miseq','x_miseq','y_miseq')
  miseq_pos$tile_miseq=paste(miseq_pos$lane_miseq,miseq_pos$tile_miseq,sep="_")
  print(head(miseq_pos))
  if (tiles=='All')
  {
    tiles=unique(miseq_pos$tile_miseq)
    #tile_ind='All'

  }
  print('tiles')
  print(tiles)
  tiles = intersect(tiles,unique(miseq_pos$tile_miseq))
  bottom=miseq_pos[miseq_pos$tile_miseq %in% tiles,]

  if (seqscope1st=='MiSeq')
  {
    #bottom = miseq_pos[miseq_pos$tile>2100,]
    plotwidth = plotheight=3.5
    nrow=1
    ncol=length(tiles)
  }
  if (seqscope1st=='HiSeq')
  {
    #bottom = miseq_pos
    plotheight=3.5
    plotwidth=plotheight*3
    if (length(tiles)==1)
    {
      nrow=ncol=1
    }
    else
    {
      nrow=2
      ncol=ceiling(length(tiles)/nrow)
    }
  }
  if (seqscope1st=='Custom')
  {
    if (layout=='False' & length(tiles)==1)
    {
      nrow=ncol=1

    }
    print('fix later')
  }

  print('merge')

  df = data.frame("HDMI" =colnames(m),"HDMIind" = 1:(dim(m)[2]))

  tile_df = merge(bottom,df,by = "HDMI")

  tile_df$tile=as.factor(tile_df$tile)  #aggregate all tils by expanding coord
  #x=unique(tile_df$tile)
  #ordering=data.frame('tile'=levels(tile_df$tile),'aggrInd'=order(as.factor(x))-1)
  #tile_df=merge(tile_df,ordering,by='tile')

  #tile_df$aggrInd =  as.numeric(factor(tile_df$tile_miseq))-1
  #tile_df$aggrInd = tile_df$tile_miseq - min(tile_df$tile_miseq)
  m_tile = m[,tile_df$HDMIind]
  tile_df$UMI = colSums(m_tile)
  tile_df$tileHDMIind= match(tile_df$HDMI,colnames(m_tile))



  setwd(outpath)
  write.csv(tile_df,'tile_df.csv')

  print('Start collapsing')
  obj=sapply(tiles,collapseTiles, tile_df=tile_df,binx=binx,biny=biny,m_tile=m_tile)
  obj=mergeSeuratObj(obj)
  print(head(obj@meta.data))

  tile_df = obj@meta.data
  rownames(tile_df)=colnames(obj)
  tile_df$tile=as.factor(tile_df$tile)                                       #aggregate all tils by expanding coord
  x=unique(tile_df$tile)
  if(tile_ind== 'All')
  {
    print('All 128')

    lane12=c(paste('1',1101:1116,sep='_'),paste('1',1201:1216,sep='_'),paste('1',2101:2116,sep='_'),paste('1',2201:2216,sep='_'),paste('2',1101:1116,sep='_'),paste('2',1201:1216,sep='_'),paste('2',2101:2116,sep='_'),paste('2',2201:2216,sep='_'))
  }
  if (tile_ind !='All')
  {
    lane12=tiles
  }
  #ordering=data.frame('tile'=levels(tile_df$tile),'aggrInd'=order(as.factor(x))-1)
  ordering=data.frame('tile'=lane12,'aggrInd'=order(as.factor(lane12))-1)
  tile_df=merge(tile_df,ordering,by='tile')

  #generating super tiles
  if (layout=='FALSE')
  {
    print('layout is false')
    if(seqscope1st == 'HiSeq')
    {
      print('hiseq')
      print(nrow)
      print(ncol)
      nrow=2
      ncol=ceiling(length(lane12)/nrow)
      addson_hori = max(tile_df$Y)
      addson_verti = max(tile_df$X)
      # tile_df$aggrInd =  as.numeric(factor(tile_df$tile))-1
      tile_df$aggrInd2 = tile_df$aggrInd
      tile_df$y_miseq_expand = tile_df$Y +  addson_hori*(((tile_df$aggrInd2%%ncol)))
      tile_df$x_miseq_expand =   tile_df$X +  addson_verti*(floor((tile_df$aggrInd2/ncol)))


    }
    if(seqscope1st == 'MiSeq')
    {
      addson_hori = max(tile_df$Y)
      addson_verti = max(tile_df$X)
      # tile_df$aggrInd =  as.numeric(factor(tile_df$tile))-1
      tile_df$aggrInd2 = tile_df$aggrInd
      tile_df$y_miseq_expand = tile_df$Y +  addson_hori*tile_df$aggrInd2
      tile_df$x_miseq_expand =   tile_df$X

    }
    if (seqscope1st == 'Custom')
    {
      tile_df$y_miseq_expand = tile_df$Y
      tile_df$x_miseq_expand =   tile_df$X
    }

  }

  else
    {
      print('layout')
      layout = read.csv(layout,row.names=1)
      nrow = max(layout$ROW)
      ncol = max(layout$COL)
     # layout$tile=  paste(layout[,3],layout[,4],sep="_")
      layout$tile=paste(layout$LANE,layout$TILE,sep="_")
     #layout$ROW_COL=paste(layout$ROW,layout$COL,sep="_")
      #layout$aggrInd =  as.numeric(factor(layout$ROW_COL))-1
      print(head(layout))

      print(head(tile_df))
      #tile_df$aggrInd2 = layout$aggrInd
      tile_df = merge(tile_df,layout,by='tile')
      addson_hori = max(tile_df$Y)
      addson_verti = max(tile_df$X)
      tile_df$y_miseq_expand = tile_df$Y +  addson_hori*(tile_df$COL-1)
      tile_df$x_miseq_expand = tile_df$X +  addson_verti*(tile_df$ROW-1)

    }

#
#   if (layout=='FALSE')
#   {
#     if(order=='top')
#     {
#
#       addson_hori = max(tile_df$Y)
#       addson_verti = max(tile_df$X)
#      # tile_df$aggrInd =  as.numeric(factor(tile_df$tile))-1
#       tile_df$aggrInd2 = tile_df$aggrInd
#       tile_df$y_miseq_expand = tile_df$Y +  addson_hori*(((tile_df$aggrInd2%%ncol)))
#       tile_df$x_miseq_expand =   tile_df$X +  addson_verti*(floor((tile_df$aggrInd2/ncol)))
#     }
#     if(order=='bottom')
#     {
#
#       addson_hori = max(tile_df$Y)
#       addson_verti = max(tile_df$X)
#       tile_df$aggrInd =  as.numeric(factor(tile_df$tile))-1
#
#       tile_df$aggrInd2 = 0
#       tile_df[tile_df$aggrInd<ncol,"aggrInd2"] = tile_df[tile_df$aggrInd<ncol,]$aggrInd+ncol
#       tile_df[tile_df$aggrInd>=ncol,"aggrInd2"] = tile_df[tile_df$aggrInd>=ncol,]$aggrInd-ncol
#       tile_df$y_miseq_expand = tile_df$Y +  addson_hori*(((tile_df$aggrInd2%%ncol)))
#       tile_df$x_miseq_expand =   tile_df$X +  addson_verti*(floor((tile_df$aggrInd2/ncol)))
#     }
#     #for super tile
#
#   }
#   else
#   {
#     print('layout')
#     layout = read.csv(layout,row.names=1)
#     nrow = max(layout$ROW)
#     ncol = max(layout$COL)
#    # layout$tile=  paste(layout[,3],layout[,4],sep="_")
#     layout$tile=paste(layout$LANE,layout$TILE,sep="_")
#    #layout$ROW_COL=paste(layout$ROW,layout$COL,sep="_")
#     #layout$aggrInd =  as.numeric(factor(layout$ROW_COL))-1
#     print(head(layout))
#
#     print(head(tile_df))
#     #tile_df$aggrInd2 = layout$aggrInd
#     tile_df = merge(tile_df,layout,by='tile')
#     addson_hori = max(tile_df$Y)
#     addson_verti = max(tile_df$X)
#     tile_df$y_miseq_expand = tile_df$Y +  addson_hori*(tile_df$COL-1)
#     tile_df$x_miseq_expand = tile_df$X +  addson_verti*(tile_df$ROW-1)
#
#
#   }

                                        #tile_df$orig.ident = rownames(tile_df)
  rownames(tile_df)=colnames(obj)
  print(head(tile_df))
  obj@meta.data=tile_df

  obj@meta.data$X_expand = obj@meta.data$x_miseq_expand
  obj@meta.data$Y_expand = obj@meta.data$y_miseq_expand
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

#
#                                         #clustering
#   print('feature eplot')
#   png("nFeatureplot.png",width=7*2,height=6,res=300,units='in')
#   p1=VlnPlot(obj, features = "nFeature_Spatial", pt.size = 0) + NoLegend()
#   p2=VlnPlot(obj, features = "nFeature_Spatial", pt.size = 0,log=T) + NoLegend()
#   plot_grid(p1,p2,ncol=2)
#   dev.off()
#   #png("nFeatureplot",width=7,height=6,res=300,units='in')
#   #VlnPlot(obj, features = "nFeature_Spatial", pt.size = 0,log=T) + NoLegend()
#                                         #dev.off()
#
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
#   if(clustering)
#   {
#
#    # geneCount1=median(obj@meta.data$nFeature_Spatial)  #?automatic ??????
#    # print('geneCount')
#    # print(geneCount1)
#     obj_simple = subset(obj, subset = nFeature_Spatial>geneCount1)
#     obj_simple = SCTransform(obj_simple, assay = "Spatial")
#     obj_simple = RunPCA(obj_simple)
#     obj_simple = RunUMAP(obj_simple, dims = 1:30)
#     obj_simple = FindNeighbors(obj_simple, dims = 1:30)
#     obj_simple = FindClusters(obj_simple)
#     saveRDS(obj_simple,'SimpleSquareGridsWithClustering.RDS')
#     meta=obj_simple@meta.data
#     meta=cbind(meta,obj_simple@reductions$umap@cell.embeddings)
#     print(nrow)
#     print(ncol)
#     png('SpatialDimPlot.png',width=ncol*7,height=nrow*6,res=300,units='in')
#     ggplot(meta,aes(X_expand,Y_expand,color=seurat_clusters))+geom_point(size=1,alpha=1)
#     dev.off()
#
#     #  png('SpatialDimPlot.png',width=7,height=6,res=300)
#     #SpatialDimPlot(obj_simple, stroke = 0)
#     #dev.off()
#   }

  print('Done!')

}



getSimpleGrid(seqscope1st,DGEdir,spatial,tiles,nrow,ncol,sidesize,outpath,collapsePath,layout,order,nMax,clustering)






