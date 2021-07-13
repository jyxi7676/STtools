

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


aggregateSliding=function(workingdir)
{
 packages = c("purrr","Matrix", "tictoc","Seurat","rlist")
  ## add more packages to load if needed
  ## Now load or install&load all
  package.check = lapply(
    packages,
    FUN = function(x) {
      if (!require(x, character.only = TRUE)) {
        install.packages(x, repos="https://cran.rstudio.com", dependencies = TRUE)
        library(x, character.only = TRUE)
      }
    }
  )
 path=workingdir
 files = list.files(workingdir, pattern = "\\sliding.RDS$", full.names = TRUE)
 f= lapply(files, readRDS)
 obj=mergeSeuratObj(f)
 tile_df = obj@meta.data
 addson_hori = max(tile_df$Y)
 addson_verti = max(tile_df$X)
 tile_df$aggrInd =  as.numeric(factor(tile_df$tile))-1

 tile_df$aggrInd2 = 0
 tile_df[tile_df$aggrInd<ncol,"aggrInd2"] = tile_df[tile_df$aggrInd<ncol,]$aggrInd+ncol
 tile_df[tile_df$aggrInd>=ncol,"aggrInd2"] = tile_df[tile_df$aggrInd>=ncol,]$aggrInd-ncol

 tile_df$y_miseq_expand = tile_df$Y +  addson_hori*(((tile_df$aggrInd2%%ncol)))
 tile_df$x_miseq_expand =   tile_df$X +  addson_verti*(floor((tile_df$aggrInd2/ncol)))
  #tile_df$orig.ident = rownames(tile_df)
 obj@meta.data$X_expand = tile_df$x_miseq_expand
 obj@meta.data$Y_expand = tile_df$y_miseq_expand

 obj@images$image = new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = obj@meta.data[,c('Y_expand','X_expand')]
  )
 saveRDS('SlidingGrids.RDS')
 print('Done')
}
