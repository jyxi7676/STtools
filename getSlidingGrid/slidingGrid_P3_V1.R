#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
for(i in 1:11) { #-- Create objects  'r.1', 'r.2', ... 'r.6' --
  nam <- paste0("r", i)
  assign(nam, args[i])
}

####################################################################################################3
#' This function merges a list of Seurat object
#' @param seurat_object_list a list of Seurat object
mergeSeuratObj=function(seurat_object_list)
{

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

mergeTileSubFieldRds=function(outpath,ncol,nrow)
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
  df = list.files(pattern = ".RDS") %>%
    map(readRDS)


  obj = mergeSeuratObj(df)

  #obj_list[[as.character(tile)]] = obj_tile
  #obj= mergeSeuratObj(obj_list)

  objfile=paste0('SlidingSquareGrids.RDS')
  # if (file.exists(objfile))
  #   file.remove(objfile)
  junk = dir(path=outpath,  pattern=".RDS")
  junk=junk[!junk %in% 'SimpleSquareGrids.RDS']
  file.remove(junk)
  junk = dir(path=outpath,  pattern=".csv")
  file.remove(junk)
  junk = dir(path=outpath,  pattern="group_tile")
  file.remove(junk)
  junk = dir(path=outpath,  pattern="m_tile")
  file.remove(junk)
  junk = dir(path=outpath,  pattern="Temp")
  file.remove(junk)
  


  tile_df = obj@meta.data
  addson_hori = max(tile_df$Y)
  addson_verti = max(tile_df$X)
  if(length(unique(tile_df$tile))==1)
  {
    obj@meta.data$X_expand=obj@meta.data$X
    obj@meta.data$Y_expand=obj@meta.data$Y
  }
  else
  {
    tile_df$aggrInd =  as.numeric(factor(tile_df$tile))-1

    tile_df$aggrInd2 = 0
    tile_df[tile_df$aggrInd<ncol,"aggrInd2"] = tile_df[tile_df$aggrInd<ncol,]$aggrInd+ncol
    tile_df[tile_df$aggrInd>=ncol,"aggrInd2"] = tile_df[tile_df$aggrInd>=ncol,]$aggrInd-ncol

    tile_df$y_miseq_expand = tile_df$Y +  addson_hori*(((tile_df$aggrInd2%%ncol)))
    tile_df$x_miseq_expand =   tile_df$X +  addson_verti*(floor((tile_df$aggrInd2/ncol)))
    #tile_df$orig.ident = rownames(tile_df)
    obj@meta.data$X_expand = tile_df$x_miseq_expand
    obj@meta.data$Y_expand = tile_df$y_miseq_expand
  }


  obj@images$image = new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = obj@meta.data[,c('Y_expand','X_expand')]
  )

  saveRDS(obj,objfile)
  print('Done!')
}


mergeTileSubFieldRds(r1,as.numeric(r2),as.numeric(r3))

#function(seqscope1st,DGEdir,spatial,nrow,ncol,sidesize,outpath,window,collapsePath,slidingPath,tile)
#getSlidingGrid(r1,r2,r3,as.numeric(r4),as.numeric(r5),as.numeric(r6),r7,as.numeric(r8),r9,r10,as.numeric(r11))

