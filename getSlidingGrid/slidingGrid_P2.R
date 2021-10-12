#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
for(i in 1:11) { #-- Create objects  'r.1', 'r.2', ... 'r.6' --
  nam <- paste0("r", i)
  assign(nam, args[i])
}

r11=as.numeric(as.vector((strsplit(r11, ",")[[1]])))


####################################################################################################3
#' This function generate grids in subfields of a tile
#' @param tile_df dataframe with spatial info of a certain tile
#' @param m_tile sparise matrix
#' @param slidestarts  a vector of the sliding steps
#' @param window sliding window size
#' @param binx x side size of the simple grid
#' @param biny y side size of the simple grid
#' @param tile integer of the tile number
getSubGrids = function(tile_df,m_tile,slidestarts,window,binx,biny,tile,collapsePath)
{
                                        #  tic()
  df_miny = min(tile_df$y_miseq)
  df_maxy = max(tile_df$y_miseq)
  df_minx = min(tile_df$x_miseq)
  df_maxx = max(tile_df$x_miseq)
  obj=suppressWarnings(CreateSeuratObject(m_tile[,1:2]))#initialize note if none
  #tile_df_sub= tile_df
  #Use nested loop instead
  #obj=slidingWindow(slidestarts,tile_df_sub,window,binx,biny,simpleGridobj,m_tile,colnames(m_tile),rownames(m_tile),tile,groupid)
  for(j in slidestarts)
  {
    print(j)
    for(t in slidestarts)
    {
     print(t)
     print('start simple grids')
                                        #tic()
     #print(dim(tile_df))
     obj=simpleGrid(tile_df,binx,biny,window,j,t,obj,m_tile,colnames(m_tile),rownames(m_tile),tile,collapsePath,df_miny,df_maxy,df_minx,df_maxx);
     #toc()
    }
  }
 # toc()
  return(obj)

}

####################################################################################################3
#' This function merges a list of Seurat object
#' @param seurat_object_list a list of Seurat object
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
#' This function generate simple square grids of subfield data with user-defined binning size and collapsed HDMIs within each grid
#' @param tile_df_sub dataframe of subfied data
#' @param binx binning size of x axis
#' @param biny binning size of y axis
#' @param window sliding window size
#' @param j integer of horizontal sliding step
#' @param t interger of vertical sliding step
#' @param obj seurat object
#' @param m_tile sparse matrix, maybe made global variable?
#' @param m_bc barcodes of m_tile
#' @param m_gene barcodes of m_gene
#' @param tile integer of one tile
simpleGrid = function(tile_df_sub,binx,biny,window,j,t,obj,m_tile,m_bc,m_gene,tile,collapsePath,df_miny,df_maxy,df_minx,df_maxx)
{
   miny = min(tile_df_sub$y_miseq)
   maxy = max(tile_df_sub$y_miseq)
   minx = min(tile_df_sub$x_miseq)
   maxx = max(tile_df_sub$x_miseq)
  #print('subsetting')

  #tic()
  tile_df_sub_down = tile_df_sub[tile_df_sub$x_miseq>=(minx+window*j)&tile_df_sub$x_miseq<=(maxx-window*j),]
  #print('dim')
  #print(dim(tile_df_sub))
  #print(dim(tile_df_sub_down))
  if(dim(tile_df_sub_down)[1]==0)
  {
    return(NULL)
  }
  miny = min(tile_df_sub_down$y_miseq)
  maxy = max(tile_df_sub_down$y_miseq)
  minx = min(tile_df_sub_down$x_miseq)
  maxx = max(tile_df_sub_down$x_miseq)
  tile_df_sub_wind = tile_df_sub_down[tile_df_sub_down$y_miseq>=(miny+window*t) & tile_df_sub_down$y_miseq<=(maxy-window*t),]
  #print(dim(tile_df_sub_wind))
  if(dim(tile_df_sub_wind)[1]==0)
  {
    return(NULL)
  }
  
  xlim = c(min(tile_df_sub_wind$x_miseq),max(tile_df_sub_wind$x_miseq))
  ylim = c(min(tile_df_sub_wind$y_miseq),max(tile_df_sub_wind$y_miseq))
  print('xlim')
  print(xlim)
  print(ylim)
   print(binx)
                                        #print(xlim[0])
   gapx=xlim[2]-xlim[1]
   gapy=ylim[2]-ylim[1]
  if(gapx==0)
  {
  return(NULL)
  }
  if(gapy==0)
 {
 return(NULL)
 }
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

 #toc()
 # print('make grid')
                                        # tic()
   #print('write')
  # write.csv(tile_df_sub_wind,'tile_df_sub_wind.csv')
 print('xlim')
  print(xlim)
  print(ylim)
   print(binx)
  grd = make.grid(tile_df_sub_wind$x_miseq,tile_df_sub_wind$y_miseq,tile_df_sub_wind$UMI, binx,biny, xlim, ylim)
 #print(dim(grd))
  grd=t(grd)

  fn1=paste0('Temp_CollapsedHDMIsIndLength',tile,'_',j,t,'.csv')
  fn2=paste0('Temp_CollapsedHDMIsInd',tile,'_',j,t,'.txt')
  
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
  test.env = new.env()

  test.env$len =c()
  test.env$no = c()
  grd_re = make.grid(tile_df_sub_wind$x_miseq,tile_df_sub_wind$y_miseq,tile_df_sub_wind$tileHDMIind, binx, biny, xlim, ylim,function(x) {test.env$len=c(test.env$len,length(x))})
  grd_re = make.grid(tile_df_sub_wind$x_miseq,tile_df_sub_wind$y_miseq,tile_df_sub_wind$tileHDMIind, binx, biny, xlim, ylim,function(x) {test.env$no=c(test.env$no,x)})  

  collapseLen = test.env$len
  collapseInd = test.env$no

  collapseLen = cbind(collapseLen,cumsum(collapseLen))
  colnames(collapseLen) =c("len","end")
  collapseLen=as.data.frame(collapseLen)
  interv = c(0,collapseLen$end)
 # toc()
  #print(toc())
  #create dataframe of the assignment of collapsed grids for each HDMI
  #print('matrix subsetting')
  #tic()
  df=data.frame('HDMIind' = collapseInd,"HDMI" = m_bc[collapseInd])
  #toc()
  assign=c()
  #print('add name')
  #tic()
  out=sapply(1:(dim(collapseLen)[1]),function(x) {nrep=collapseLen$len[x];return(c(assign,rep(paste0('Collapse2_',x),nrep)))})
  assign = unlist(out)
  df$assign = assign
  #toc()
                                        #print(toc())

                                        #toc()
  tic()
  sourceCpp(collapsePath)
  print('source sucesful')
  #tic();
  collapseM = collapse(m_tile,collapseInd,interv)
  #  toc()
  #print('finished collapse!')

  
  #print(tic())
  #print('last')
  
  rownames(collapseM) = m_gene
  colnames(collapseM) = paste0("Collase_tile_",tile,"_",j,"_",t,"_",1:(length(interv)-1))
  sparse.gbm = Matrix(collapseM , sparse = T )

  toc()
  grd=t(grd)
  pos = which(!is.na(grd), TRUE)
  pos_coor = t(sapply(1:(dim(pos)[1]),function(x) {c(as.numeric(rownames(grd)[as.numeric(pos[x,1])]),as.numeric(colnames(grd)[as.numeric(pos[x,2])]) )})) #here the grd is transposed
  colnames(pos_coor) =c("X","Y")
  coord.df = data.frame("Y"=pos_coor[,2], "X"=pos_coor[,1],"tile"=tile, stringsAsFactors=FALSE)
  #fn7=paste0('tile2112_coor_df_',j,'_',t,'.csv')
  #write.csv(coord.df,fn7)

  dge1 = collapseM
  spatial1 = coord.df
  obj1 = suppressWarnings(expr=CreateSeuratObject(counts=dge1,assay='Spatial'))
  #obj1$status = "Original"
  obj1@meta.data$tile = spatial1$tile
  obj1@meta.data$X = spatial1$X
  obj1@meta.data$Y = spatial1$Y
  obj1@meta.data$interationi=j
  obj1@meta.data$interationj=t
  #print(obj1)

  if(j==0 & t==0)
  {
    obj=obj1
  }
  else
  {
    obj = merge(obj,obj1)

  }
  #print(toc())
  #toc()
  return(obj)
}






slidingWindowSub=function(collapsePath,DGEdir,outpath,window,sidesize,xargs)
{
  packages = c("Rcpp","Matrix", "tictoc", "ggplot2", "ggsci","Seurat","mapplots","rlist","cowplot","dplyr")
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

  if (!dir.exists(DGEdir)){
    stop("DGE  directory does not exist")
  }
  if (!dir.exists(outpath)){
    stop("Output path does not exist")
  }
  if(missing(sidesize))
  {
    sidesize=300
  }
  biny = binx = sidesize
  slidestarts = seq(0,(binx/window-1),1)
  setwd(DGEdir)
  features = read.table('features.tsv',header=F)$V2


  setwd(outpath)
  print('xargs')
  print(xargs)
  # pat = "(.*?_.*?)_(.*)"
  # tile=sub(pat,"\\1",xargs)
  # print('tile')
  # print(tile)
  # groupid=sub(pat,"\\2",xargs)
  # print('groupid')
  # print(groupid)
  # print('stop')
  tile=xargs
 # tile=strsplit(xargs,'_')[[1]][1]
 # groupid=strsplit(xargs,'_')[[1]][2]
  tilefile=paste0('m_tile_',tile)
  print(tilefile)
  m_tile_sub=readMM(tilefile)

  sub=read.csv(paste0('group_tile_',tile,'.csv'))
  sub_xargs=sub[sub$tile_miseq==xargs,]
  colnames(m_tile_sub)=sub_xargs$HDMI
  rownames(m_tile_sub)=features
  print('collapsing')
  sub_xargs$tileHDMIind=1:(dim(sub_xargs))[1]
  print('Pring tiles')
#  print(tile)
 # print(as.numeric(tile))
  #df_miny = min(sub_xargs$y_miseq)
  #df_maxy = max(sub_xargs$y_miseq)
  #df_minx = min(sub_xargs$x_miseq)
  #df_maxx = max(sub_xargs$x_miseq)
  #print('p2 time')
                                        #tic()
  print('subxargs')
  print(head(sub_xargs))
  print(dim(sub_xargs))
  clp=getSubGrids(sub_xargs,m_tile_sub,slidestarts,window,binx,biny,tile,collapsePath)
  #toc()
  saveRDS(clp,paste0('tile_',tile,'.RDS'))
  #out=sapply(group,getSubGrids,tile_df=tile_df_exact,m_tile=m_tile,slidestarts=slidestarts,window=window,binx=binx,biny=biny,tile=tile)  #function applies on subfied
  print('finish subfield collapsing')
}


slidingWindowSub(r1,r2,r3,as.numeric(r4),as.numeric(r5),r6)




