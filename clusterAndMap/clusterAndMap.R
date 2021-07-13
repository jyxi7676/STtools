#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
for(i in 1:6) { #-- Create objects  'r1', 'r2', ... 'r6' --
  nam <- paste0("r", i)
  assign(nam, args[i])
}



#function to exist R program
exit <- function() { invokeRestart("abort") }
#' This function runs clustering and mapping of simple and sliding grids
#' @param workingdir working directory
#' @param obj1_path path to seurat object of simple square grids
#' @param obj2_path path to seurat object of sliding square grids
#' @param geneCount1 cutoff of nFeature in simple square grids
#' @param geneCount2 cutoff of nFeature in sliding square grids
#' @nFeaturePlotOnly TRUE or FALSE. If TRUE, the program end after violin plot. Otherwise, run the whole program
#' @export
runClustering=function(workingdir,obj1_path,obj2_path,geneCount1=0,geneCount2=0,nFeaturePlotOnly="FALSE",outpath)
{
  setwd(workingdir)
  #libraries
  packages = c("ggplot2", "ggsci","Seurat","rlist","cowplot","tictoc")
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

  memory.limit(size=1e13)
  if(file.exists(obj1_path)==F)
  {
    print("Obj1 does not exit")
    break
    }

  if(file.exists(obj2_path)==F)
  {
    print("Obj2 does not exit")
    break
  }

  #obj_simple=readRDS('C:/Users/jyxi/Downloads/SimpleSqureGrids.RDS')
  #obj_sliding=readRDS('C:/Users/jyxi/Downloads/SlidingSquareGrids_tile_2106.RDS')
  obj_simple=readRDS(obj1_path)
  obj_sliding=readRDS(obj2_path)

  v1=VlnPlot(obj_simple, features = "nFeature_Spatial", pt.size = 0) + NoLegend()
  v2=VlnPlot(obj_sliding, features = "nFeature_Spatial", pt.size = 0) + NoLegend()
  png("nFeatureplot",width=7*2,height=6,res=300,units='in')
  plot_grid(v1,v2,ncol=2)
  dev.off()
  if(nFeaturePlotOnly!="FALSE")
    {
      print("Done")
      exit()
  }
  else
  {
    obj_simple = subset(obj_simple, subset = nFeature_Spatial>geneCount1)
    obj_simple = SCTransform(obj_simple, assay = "Spatial")
    obj_simple = RunPCA(obj_simple)
    obj_simple = RunUMAP(obj_simple, dims = 1:30)
    obj_simple = FindNeighbors(obj_simple, dims = 1:30)
    obj_simple = FindClusters(obj_simple)

    obj_sliding = subset(obj_sliding, subset = nFeature_Spatial>geneCount2)
    obj_sliding = SCTransform(obj_sliding, assay = "Spatial")
    obj_sliding = RunPCA(obj_sliding)
    obj_sliding = RunUMAP(obj_sliding, dims = 1:30)

    #DimPlot(obj_sliding, label = TRUE)
    print('find anchor')
    anchors = FindTransferAnchors(reference = obj_simple, query = obj_sliding, normalization.method = "SCT")
    predictions.assay = TransferData(anchorset = anchors, refdata = obj_simple@active.ident, prediction.assay = TRUE,
                                     weight.reduction = obj_sliding[["pca"]], dims = 1:30)
    obj_sliding[["predictions"]] = predictions.assay
    DefaultAssay(obj_sliding) = "predictions"

    obj_sliding$predicted.id = GetTransferPredictions(obj_sliding, score.filter = 0)
    Idents(obj_sliding) = "predicted.id"

    print('plotting')
    png("dimplot.png",width=7*2,height=6,res=300,units='in')
    p1=DimPlot(obj_simple, label = TRUE)
    p2=DimPlot(obj_sliding,label = TRUE)
    plot_grid(p1,p2,ncol=2)
    dev.off()
#    print('almost')
    saveRDS(obj_simple,'simpleGrid_clus.RDS')
    saveRDS(obj_sliding,"slidingGrid_mapping.RDS")
    #insert coding here after fixing slidingWindow function
    #SpatialDimPlot(obj_sliding)


    print("Done")
  }
}


runClustering(r1,r2,r3,as.numeric(r4),as.numeric(r5),r6)


