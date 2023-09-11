
setwd(paste0("~/experimental/My projects/Flowering variability RNA-seq",
             "/desync_gene_expr"))

# Libraries used for ALL sections ==========

# plotting
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(plotly)
library(gridExtra)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(viridis)
library(ggbeeswarm)
library(ggrepel)
library(RColorBrewer)
library(VennDiagram)
library(dendsort)
library(magick) # for rasterization of heatmaps

# data input / manipulation
library(readxl)
library(xlsx)
library(matrixStats)
library(tximport)

# for saving plots as svg
library(svglite)
library(stats)
library(graphics)

#for t tests and other hypothesis tests
library(ggpubr)
library(rstatix)
library(PMCMRplus)
library(multcompView)

#for smoothing
library(fda)

# for GO / KEGG terms
library(biomaRt)
library(gprofiler2)
library(KEGGREST)
library(reshape2)

# for graph and network analysis
library(GENIE3)
library(igraph)
library(WGCNA)

# CONSTANTS ========================================

# Allows for consistent plotting throughout
no_col <- "#E65719"
yes_col <- "#19A8E6"

# The path to DynGENIE3 executable needs to be called
dyngenie3_path <- paste0(getwd(),
                         "/external_scripts/dynGENIE3_compiled/dynGENIE3.dll")

# Functions for use throughout =========

get_z_score_matrix <- function(df_in, include_na = FALSE) {
  # Input should be a data frame containing only numeric variables
  # Use include_na = TRUE to include rows with 0 standard deviation

  # Calculate z scores - i.e. transform each row by subtracting mean
  # then dividing by the standard deviation
  df_z <- (df_in - rowMeans(df_in)) / (rowSds(as.matrix(df_in)))[row(df_in)]

  # remove values with Na - i.e. anything where the standard deviation was 0
  if (!include_na) {
    df_z <- df_z %>% drop_na()
  }

  return(df_z)
}

get_01_matrix <- function(df_n, include_na = FALSE) {
  # Input should be a data frame containing only numeric variables
  # Use include_na = TRUE to include rows with 0 range

  df_as_mat <- as.matrix(df_n)

  mat_min <- rowMins(df_as_mat)
  mat_max <- rowMaxs(df_as_mat)
  mat_range <- mat_max - mat_min

  df01 <- (df_n - mat_min) / (mat_range)[row(df_n)]

  if (!include_na) {
    df01 <- df01 %>% drop_na()
  }

  return(df01)
}

# from https://stackoverflow.com/questions/47916307/specify-position-of-geom-text-by-keywords-like-top-bottom-left-right
annotation_compass <- function(label,
                               position = c('N','NE','E','SE','S','SW','W','NW'),
                               padding = grid::unit(c(0.5,0.5),"line"), ...){
  position <- match.arg(position)
  x <- switch (position,
               N = 0.5,
               NE = 1,
               E = 1,
               SE = 1,
               S = 0.5, 
               SW = 0,
               W = 0, 
               NW = 0
  )
  y <- switch (position,
               N = 1,
               NE = 1,
               E = 0.5,
               SE = 0,
               S = 0, 
               SW = 0,
               W = 0.5, 
               NW = 1
  )
  hjust <- switch (position,
                   N = 0.5,
                   NE = 1,
                   E = 1,
                   SE = 1,
                   S = 0.5, 
                   SW = 0,
                   W = 0, 
                   NW = 0
  )
  vjust <- switch (position,
                   N = 1,
                   NE = 1,
                   E = 0.5,
                   SE = 0,
                   S = 0, 
                   SW = 0,
                   W = 0.5, 
                   NW = 1
  )
  f1 <- switch (position,
                N = 0,
                NE = -1,
                E = -1,
                SE = -1,
                S = 0, 
                SW = 1,
                W = 1, 
                NW = 1
  )
  f2 <- switch (position,
                N = -1,
                NE = -1,
                E = 0,
                SE = 1,
                S = 1, 
                SW = 1,
                W = 0, 
                NW = -1
  )
  annotation_custom(grid::textGrob(label, 
                                   x = grid::unit(x,"npc") + f1*padding[1] , 
                                   y = grid::unit(y,"npc") + f2*padding[2],
                                   hjust = hjust, vjust = vjust, ...))
}