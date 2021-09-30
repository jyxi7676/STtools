#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
for(i in 1:5) { #-- Create objects  'r.1', 'r.2', ... 'r.6' --
  nam <- paste0("r", i)
  assign(nam, args[i])
}

####################################################################################################3
#' This function generate grids in subfields of a tile
#' @param groupid an integer indicates the id for the subfield
#' @param tile_df dataframe with spatial info of a certain tile
#' @param m_tile sparise matrix
#' @param slidestarts  a vector of the sliding steps
#' @param window sliding window size
#' @param binx x side size of the simple grid
#' @param biny y side size of the simple grid
#' @param tile The lane_ tile number
getSubGrids = function(groupid,tile_df,m_tile,slidestarts,window,binx,biny,tile)
{


  #print(groupid)
  obj=suppressWarnings(CreateSeuratObject(m_tile[,1:2]))#initialize
  tile_df_sub=tile_df
  #tile_df_sub= tile_df[tile_df$ID==groupid,]
  #sourceCpp(slidingPath) put into the main function
  #Use nested loop instead
  #obj=slidingWindow(slidestarts,tile_df_sub,window,binx,biny,simpleGrid,obj,m_tile,colnames(m_tile),rownames(m_tile),tile,groupid)
  for(j in slidestarts)
  {

    print(dim(m_tile))
    for(t in slidestarts)
    {
      obj=simpleGrid(tile_df_sub,binx,biny,window,j,t,obj,m_tile,colnames(m_tile),rownames(m_tile),tile,groupid);
    }
  }
  return(obj)

}


getGroupGrids=function(tiles,tile_df,m_tile,slidestarts,window,binx,biny)
{
  tile=tiles[i]

  submat=sapply(1:length(tiles),function(x) {tile_df_=tile_df[tile_df$tile_miseq==tiles[x],];ind=match(tile_df_$HDMI,colnames(m_tile));m_tile_=m_tile[,ind];writeMM(m_tile_,paste0('m_tile_',tiles[x]));file=paste0('group_tile_',tiles[x],'.csv');  tile_df_$tileHDMIind=1:dim(tile_df_)[1];write.csv(tile_df_,file,append=F)})
  group = unique(tile_df$tile_miseq)
  group = group[order(unique(tile_df$tile_miseq))]
  write.table(group,paste0('groupgrids_tile','.txt'),row.names=FALSE,sep="\t",col.names=F,quote=F,append=T)
}



####################################################################################################3
#' This function grids the data with user-defined binning size and collapsed HDMIs within each grid
#' @param seqscope1st Data flatform, please use "MiSeq" for now
#' @param DGEdir folder that stores barcodes.tsv, features.tsv and matrix.mtx
#' @param spatial txt file stores spatial informaiton with four columns: 'HDMI','tile_miseq','x_miseq','y_miseq'
#' @param outpath
#' @param tiles
#' @export
getSubfield = function(layout,DGEdir,spatial,outpath,tiles)
{
  options(warn=-1)
  if (missing(layout))
  {
    layout='MiSeq'
  }

  if (!dir.exists(DGEdir)){
    stop("DGEdir does not exist")
  }
  if (!dir.exists(outpath)){
    stop("outpath does not exist")
  }

  setwd(DGEdir)
  #libraries
  packages = c("R.utils","data.table","tidyverse","Matrix", "tictoc","rlist",'dplyr','purrr', "ggplot2", "ggsci","mapplots","dplyr")
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
  m = m[,colSums(m)<=100&colSums(m)>0]  #remove outliers
  dim(m)
  #get spatial info
  print('Reading spatial info')
  miseq_pos = read.table(spatial)
  colnames(miseq_pos) = c('HDMI','lane_miseq','tile_miseq','x_miseq','y_miseq')
  miseq_pos$tile_miseq=paste(miseq_pos$lane_miseq,miseq_pos$tile_miseq,sep="_")
  tiles=(unlist(strsplit(tiles,',')))
  tiles_uniq=unique(miseq_pos$tile_miseq)
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
    print(tiles)
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
  print(head(tile_df))
  print(tiles)
  print(dim(m_tile))



  # df = data.frame("HDMI" =colnames(m),"HDMIind" = 1:(dim(m)[2]))
  # setwd(outpath)
  # tile_df = merge(miseq_pos,df,by = "HDMI")  #this takes some time, just save this
  # m_tile = m[,tile_df$HDMIind]
  # print(dim( m_tile))
  # tile_df$UMI=colSums(m_tile)
  if (file.exists('groupgrids_tile.txt'))
  {
      file.remove('groupgrids_tile.txt')
  }
  print('sapply')
  out=getGroupGrids(tiles=tiles,tile_df=tile_df,m_tile=m_tile,slidestarts=slidestarts,window=window,binx=binx,biny=biny)

  #out=sapply(1:length(tiles),getGroupGrids,tiles=tiles,tile_df=tile_df,m_tile=m_tile,slidestarts=slidestarts,window=window,binx=binx,biny=biny)

 print('Finish generating subfields!')
}


getSubfield(r1,r2,r3,r4,r5)
