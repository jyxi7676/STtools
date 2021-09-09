#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
for(i in 1:5) { #-- Create objects  'r.1', 'r.2', ... 'r.6' --
  nam <- paste0("r", i)
  assign(nam, args[i])
}

#r5=as.numeric(as.vector((strsplit(r5, ",")[[1]])))
r5=(unlist(strsplit(r5,',')))
print(r5)
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


getGroupGrids=function(i,tiles,tile_df,m_tile,slidestarts,window,binx,biny)
{
  tile=tiles[i]
  tile_df_exact=tile_df[tile_df$tile_miseq == tile,]
  nsub_x=nsub_y=5
  tile_df_exact$xcut = cut(tile_df_exact$x_miseq,nsub_x)
  tile_df_exact$ycut = cut(tile_df_exact$y_miseq,nsub_y)
  tile_df_exact=tile_df_exact%>%mutate(ID = group_indices(., xcut, ycut))
  tile_df_exact$tile=tile
  
  #tile_df_exact$tile_groupid=paste0(tile_df_exact$tile,'_',tile_df_exact$ID)
  submat=sapply(1:25,function(x) {tile_df_=tile_df_exact[tile_df_exact$ID==x,];ind=match(tile_df_$HDMI,colnames(m_tile));m_tile_=m_tile[,ind];writeMM(m_tile_,paste0('m_tile_',tile,'_sub_',x))})
#insert gap filling
  tile_df_exact$tile_groupid=paste0(tile_df_exact$tile,'_',tile_df_exact$ID)
  tile_df_exact$tileHDMIind=1:dim(tile_df_exact)[1]

  file=paste0('group_tile_',tile,'.csv')
  write.csv(tile_df_exact,file,append=F)

  group = unique(tile_df_exact$tile_groupid)
  group = group[order(unique(tile_df_exact$ID))]
  write.table(group,paste0('groupgrids_tile','.txt'),row.names=FALSE,sep="\t",col.names=F,quote=F,append=T)

  #need to take the adjacent subfield into consideration
  
                                        #return(tile_df_exact)
  #obj = suppressWarnings(CreateSeuratObject(m_tile[,1:10]))
  #out=sapply(1:5,getSubGrids,tile_df=tile_df_exact,m_tile=m_tile,slidestarts=slidestarts,window=window,binx=binx,biny=biny,tile=tile)  #function applies on subfied
  #obj_tile = mergeSeuratObj(out)
  #obj_list[[as.character(tile)]] = obj_tile
  # return(obj_list)
}



####################################################################################################3
#' This function grids the data with user-defined binning size and collapsed HDMIs within each grid
#' @param seqscope1st Data flatform, please use "MiSeq" for now
#' @param DGEdir folder that stores barcodes.tsv, features.tsv and matrix.mtx
#' @param spatial txt file stores spatial informaiton with four columns: 'HDMI','tile_miseq','x_miseq','y_miseq'
#' @param outpath
#' @param tiles
#' @export
getSubfield = function(seqscope1st,DGEdir,spatial,outpath,tiles)
{
  options(warn=-1)
  if(missing(seqscope1st))
  {
    seqscope1st="MiSeq"
  }


  if (!dir.exists(DGEdir)){
    stop("DGEdir does not exist")
  }
  if (!dir.exists(outpath)){
    stop("outpath does not exist")
  }

  setwd(DGEdir)
  #libraries
  packages = c("Matrix", "tictoc","rlist",'dplyr','purrr', "ggplot2", "ggsci","mapplots","dplyr")
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
                                        #print(!any(tiles%in%miseq_pos$tile_miseq))
  print(tiles)
  print(head(miseq_pos))
  print(head(miseq_pos$tile_miseq))
#  print(table(miseq_pos$tile_miseq))
  if (all(tiles%in%miseq_pos$tile_miseq))
  {
    print('yes')
    miseq_pos=miseq_pos[miseq_pos$tile_miseq %in% tiles,]
  }
  else
  {
    stop('Not all tiles found')
  }
  df = data.frame("HDMI" =colnames(m),"HDMIind" = 1:(dim(m)[2]))
  setwd(outpath)
  tile_df = merge(miseq_pos,df,by = "HDMI")  #this takes some time, just save this
  m_tile = m[,tile_df$HDMIind]
  tile_df$UMI=colSums(m_tile)
  if (file.exists('groupgrids_tile.txt'))
  {
      file.remove('groupgrids_tile.txt')
  }
  
  out=sapply(1:length(tiles),getGroupGrids,tiles=tiles,tile_df=tile_df,m_tile=m_tile,slidestarts=slidestarts,window=window,binx=binx,biny=biny)

print('Finish generating subfields!')
}


getSubfield(r1,r2,r3,r4,r5)
