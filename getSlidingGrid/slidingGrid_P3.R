#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
for(i in 1:6) { #-- Create objects  'r.1', 'r.2', ... 'r.6' --
  nam <- paste0("r", i)
  assign(nam, args[i])
}
print(args)
r6=(unlist(strsplit(r6,',')))
if(length(r6)>1)
{
    print('two')
    r3=2
    r2=ceiling(length(r6)/r3)
}
if(length(r6)==1)
{
    r3=1
    r2=1
}
print('r2r3')
print(r2)
print(r3)
#print('tiles')
#print(tiles)
####################################################################################################3
#' This function merges a list of Seurat object
#' @param seurat_object_list a list of Seurat object
mergeSeuratObj=function(seurat_object_list)
{
  seurat_object_list = unlist(seurat_object_list)
  for (i in length(seurat_object_list))
  {
    seurat_object_list[[i]] = RenameCells(seurat_object_list[[i]],
                                          add.cell.id = paste0('Cell',i))
  }
  merged_combined = suppressWarnings(expr=reduce(seurat_object_list,
                                                 merge,
                                                 do.normalize = FALSE))
  return(merged_combined)
}

mergeTileSubFieldRds=function(outpath,ncol,nrow,layout,order,tiles)
{
  packages = c("Seurat","tidyverse")
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


  if (!dir.exists(outpath)){
    stop("Output path does not exist")
  }
  setwd(outpath)
  df = list.files(pattern = ".RDS") %>% map(readRDS)
  df=df[df!='SimpleSquareGrids.RDS']
  df=df[df!='SlidingSquareGrids.RDS']

  print(df)
  print('merge rds')
  obj = mergeSeuratObj(df)
  print(obj)
  print(table(obj$tile))
  print(head(obj@meta.data))
  #obj_list[[as.character(tile)]] = obj_tile
  #obj= mergeSeuratObj(obj_list)
  
  objfile='SlidingSquareGrids.RDS'
  # if (file.exists(objfile))
  #   file.remove(objfile)


  #junk = dir(path=outpath,  pattern=".RDS")
  #junk=junk[!junk %in% 'SimpleSquareGrids.RDS']
  #file.remove(junk)
  junk = dir(path=outpath,  pattern=".csv")
  file.remove(junk)
  junk = dir(path=outpath,  pattern="group_tile")
  file.remove(junk)
  junk = dir(path=outpath,  pattern="m_tile")
  file.remove(junk)
  junk = dir(path=outpath,  pattern="Temp")
  file.remove(junk)
  
 
  tile_df = obj@meta.data

  tile_df$name=colnames(obj)

  tiles_uniq = unique(tile_df$tile)
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

  #tile_df$tile=as.factor(tile_df$tile)                                       #aggregate all tils by expanding coord
  
  #ordering=data.frame('tile'=lane12,'aggrInd'=order(as.factor(lane12))-1)
 # tile_df=merge(tile_df,ordering,by='tile')

  addson_hori = max(tile_df$Y)
  addson_verti = max(tile_df$X)
  nrow = max(layout$ROW)
  ncol = max(layout$COL)
  tile_df = merge(tile_df,layout,by='tile')
  print(head(obj@meta.data))
  addson_hori = max(tile_df$Y)
  addson_verti = max(tile_df$X)
  #tile_df$y_miseq_expand = tile_df$Y +  addson_hori*(tile_df$COL-1)
  #tile_df$x_miseq_expand = tile_df$X +  addson_verti*(tile_df$ROW-1)
  tile_df$y_miseq_expand = tile_df$Y +  addson_hori*(tile_df$COL-1)
  tile_df$x_miseq_expand = tile_df$X +  addson_verti*(nrow-tile_df$ROW)
  print(head(tile_df))
  tile_df=tile_df[,c('orig.ident','name','LANE','OTILE','tile','nCount_Spatial','nFeature_Spatial','X','Y','x_miseq_expand','y_miseq_expand','ROW','COL','interationi','interationj')]
  print(head(tile_df))
  colnames(tile_df)=c('orig.ident','collapseID','lane','tile','lane_tile','nCount_Spatial','nFeature_Spatial','X','Y','X_expand','Y_expand','row','col','iter_i','iter_j')



  obj@meta.data=tile_df
  rownames(obj) = colnames(obj)

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

  saveRDS(obj,objfile)
  m=(obj@assays$Spatial@counts)
  gene=rownames(obj)
  bc=colnames(obj)
  writeMM(obj = m, file="collapsedMatrix.mtx")
  write.csv(gene,'collapsedGenes.csv')
  write.csv(bc,'collapsedBarcodes.csv')
  print('Done!')
}


mergeTileSubFieldRds(r1,as.numeric(r2),as.numeric(r3),r4,r5)

#function(seqscope1st,DGEdir,spatial,nrow,ncol,sidesize,outpath,window,collapsePath,slidingPath,tile)
#getSlidingGrid(r1,r2,r3,as.numeric(r4),as.numeric(r5),as.numeric(r6),r7,as.numeric(r8),r9,r10,as.numeric(r11))

_df) = colnames(obj)
  print(obj)

  #obj@meta.data$X_expand = obj@meta.data$x_miseq_expand
  #obj@meta.data$Y_expand = obj@meta.data$y_miseq_expand
  #print(obj)
  print(head(obj@meta.data))
  print('h3')
  obj@images$image = new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = obj@meta.data[,c('Y_expand','X_expand')]
  )
  print(obj)
  print(head(colnames(obj)))
  rownames(obj@meta.data)=colnames(obj)
  saveRDS(obj,objfile)
  m=(obj@assays$Spatial@counts)
  gene=rownames(obj)
  bc=colnames(obj)
  writeMM(m, file="collapsedMatrix.mtx")
  write.csv(gene,'collapsedGenes.csv')
  write.csv(bc,'collapsedBarcodes.csv')
  print('Done!')
}


mergeTileSubFieldRds(r1,as.numeric(r2),as.numeric(r3),r4,r5)

#function(seqscope1st,DGEdir,spatial,nrow,ncol,sidesize,outpath,window,collapsePath,slidingPath,tile)
#getSlidingGrid(r1,r2,r3,as.numeric(r4),as.numeric(r5),as.numeric(r6),r7,as.numeric(r8),r9,r10,as.numeric(r11))

