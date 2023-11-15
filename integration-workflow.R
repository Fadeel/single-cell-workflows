#!/usr/bin/env Rscript

library(Seurat)
library(yaml)
library(stringr)
library(ggplot2)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(patchwork)
set.seed(1234)

## functions
pickNumClusters <- function(object, graph = "wsnn", step = 0.05, algorithm_num = 3 , start = 0, end = 1){
  
  num.clusters = integer()
  for (res in seq(start , end, step)){
    seurat_object <- FindClusters(seurat_object,resolution = res , 
                                  graph.name = graph, algorithm = algorithm_num, verbose = FALSE)
  }
  
  i = 1
  for (res in seq(start , end, step)){
    clustering = paste0(graph,"_res.",res)
    num.clusters[i] = length( levels(seurat_object@meta.data[[clustering]]))
    i = i+ 1
  }
  
  num.clusters = as.data.frame(num.clusters)
  num.clusters$resolution <- seq(start , end, step)
  return(num.clusters)
  
}

# read command line arguments
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("A YAML parameters file is required", call.=FALSE)
} 

parameters_file <- args[1]

#parameters_file <- "params_integration.yaml"

# read parameters YAML file
params = read_yaml(parameters_file)

percent.mt = params$normalization$percent.mito
S_score = params$normalization$S_score
G2M_score = params$normalization$G2M_score

# create results folder
results_dir = paste0(params$results_dir ,params$project_name, "/")

system(sprintf("mkdir -p %s",results_dir ))

#copy parameters file to results folder
system(sprintf("cp  %s %s", parameters_file , results_dir ))

#################### read list of Seurat objects 

# add info to log file
cat(paste0("================= reading input \n"), 
    file = paste0(results_dir, "/log.txt"))

seurat_objects <- params$seurat_objects
seurat_objects <- strsplit(seurat_objects, split = ",")[[1]]

list.objs <- seurat_objects %>% lapply(readRDS)

if(length(seurat_objects) == 1){
  list.objs[[1]]$group = list.objs[[1]]@meta.data[,params$group]
  list.objs <- SplitObject(list.objs[[1]], split.by = "group")
  
  sample_ids <- unique(list.objs[[1]]$group)

}else{
  sample_ids <- strsplit(params$sample_ids, split = ",")[[1]]
}

cat(paste0("There are ", length(list.objs), " objects to be integrated\n" ), 
file = paste0(results_dir, "/log.txt"))

cat(paste0("Sample IDs:  ",paste0(sample_ids, collapse = ", ") ,"\n"), 
    file = paste0(results_dir, "/log.txt") , append = TRUE)

#################### Normalization
cat(paste0("================= Normalization \n"), 
    file = paste0(results_dir, "/log.txt") , append = TRUE)



if(!("group" %in% colnames(list.objs[[1]]@meta.data))){
  for(i in seq(length(list.objs))){
    list.objs[[i]]$group <- sample_ids[i]
    }
}

if(params$normalization$do_normalization){
  
    list.objs <- lapply(X = list.objs, FUN = function(x) {
          x <- NormalizeData(x, assay = "RNA") 
          x <- ScaleData(x, assay = "RNA", features = rownames(x)) 
       })
  
    if(params$normalization$regress_CC && "Phase" %in% colnames(list.objs[[1]]@meta.data))
    {
      list.objs <- list.objs %>% lapply(SCTransform,vars.to.regress = 
                                          c(percent.mt,S_score,G2M_score ) ) 
      
    }else{
      list.objs <- list.objs %>% lapply(SCTransform,vars.to.regress =c(percent.mt))
      }
}

#################### Integration
if(params$integration$nfeatures == "all" )
{
  features = c()
  for(i in seq(length(list.objs))){
    features <- unique(c(features,rownames(list.objs[[i]]@assays$RNA)))
  }

}else{
  features <- SelectIntegrationFeatures(object.list = list.objs, nfeatures = params$integration$nfeatures)
}

list.objs <- PrepSCTIntegration(object.list = list.objs, anchor.features = features)

list.objs <- lapply(X = list.objs, FUN = RunPCA, features = features,npcs = 30)

anchors <- FindIntegrationAnchors(object.list = list.objs, normalization.method = 
                                    params$integration$normalization.method,
                                  anchor.features = features,
                                  reduction = params$integration$reduction, 
                                  dims = 1:params$integration$dims ,k.anchor = params$integration$k.anchor )

seurat_object <- IntegrateData(anchorset = anchors, normalization.method = params$integration$normalization.method)

#################### Dimensionality reduction
cat(paste0("================= Dimensionality reduction using PCA \n"), 
    file = paste0(results_dir, "/log.txt") , append = TRUE)

DefaultAssay(seurat_object) = "RNA"
seurat_object = NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst")
seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))


DefaultAssay(seurat_object) <- "integrated"

if(params$dim_reduction$do_PCA){
  seurat_object <- RunPCA(seurat_object,npcs = 50)
  
  #pca_heatmap <- DimHeatmap(seurat_object, dims = 1:8, cells = 200, balanced = TRUE , ncol = 4,
  #                          assays = "RNA",fast = FALSE)
  elbow_plot <- ElbowPlot(seurat_object , ndims = 50)
}

#################### clustering
cat(paste0("================= Clustering \n"), 
    file = paste0(results_dir, "/log.txt") , append = TRUE)

if(params$clustering$do_clustering){
  
  seurat_object <- FindNeighbors(seurat_object, dims = 1:params$dim_reduction$n_dims,reduction = "pca",
                                 k.param = params$clustering$knn)
  
  if(params$clustering$resolution == "auto"){
    number.clusters = pickNumClusters(seurat_object, step = 0.02, end = 1, graph = 
                                        paste0(params$dim_reduction$assay,"_snn") ,algorithm_num = 1)
    t = table(number.clusters$num.clusters)
    freq_num_clusters = names(t[t==max(t)]) # most frequent number of cluster after sweeping different resolution values
    res = number.clusters$resolution[number.clusters$num.clusters == freq_num_clusters][1]
    
    p1 <- ggplot(as.data.frame(table(number.clusters$num.clusters)),
                 aes(x = Var1, y = Freq ) ) + 
      geom_point() + xlab("num. of clusters") + ylab("frequency")
    p2 <- ggplot(number.clusters, aes(x = resolution, y = num.clusters) ) + geom_point()
    num_clusters_plot <- p1 + p2 
  }else{
    res = params$clustering$resolution
  }
  
  seurat_object <- FindClusters(seurat_object,resolution = res, graph.name = 
                                  paste0(params$dim_reduction$assay,"_snn"), algorithm = 1)
  
  cat(paste0("number of cells per cluster : \n"), 
      file = paste0(results_dir, "log.txt") , append = TRUE)
  cat( table(seurat_object$seurat_clusters) , 
       file = paste0(results_dir, "log.txt") , append = TRUE)
}

#################### UMAP plot
cat(paste0("\n================= Dimensionality reduction using UMAP \n"), 
    file = paste0(results_dir, "/log.txt") , append = TRUE)

if(params$UMAP$do_UMAP){
  seurat_object = RunUMAP(seurat_object, dims = 1:params$dim_reduction$n_dims,
                          n.neighbors = params$clustering$knn)
  
  clusters_umap_plot <- DimPlot(seurat_object, reduction = "umap" , label = T, 
                                group.by = "seurat_clusters", label.size = 5)
  
  if("Phase" %in% colnames(seurat_object@meta.data)){
    cc_umap_plot <- DimPlot(seurat_object, reduction = "umap" , label = F, group.by = "Phase")
  }
  
  if("group" %in% colnames(seurat_object@meta.data)){
    sample_umap_plot <- DimPlot(seurat_object, reduction = "umap" , label = F, group.by = "group")
  }
}

#################### write Seurat object
saveRDS(seurat_object, paste0(results_dir,params$project_name , ".rds" ))


#################### cluster  markers
cat(paste0("================= Identifying cluster markers \n"), 
    file = paste0(results_dir, "/log.txt") , append = TRUE)

if(params$markers$find_markers){

    markers = FindAllMarkers(seurat_object, assay = "RNA",  min.pct = params$markers$min.pct, 
                           logfc.threshold = params$markers$logfc.threshold ,
                           only.pos = T , test.use = "LR" , latent.vars = "group" )
  
  markers = markers %>% filter(p_val_adj < params$markers$p_val_adj)
  markers = markers %>% mutate(diff.pct = pct.1 - pct.2)
  
  t = params$markers$top_markers
  cluster.specific.genes = markers %>% group_by(cluster) %>%
    slice_max(diff.pct,n = t, with_ties = FALSE)
  
  top_markers_per_clusters <- markers %>% group_by(cluster) %>%
    slice_max(avg_log2FC,n = t, with_ties = FALSE)
  
  h1 = DoHeatmap(seurat_object, features = top_markers_per_clusters$gene, assay = params$output_plots$assay) & NoLegend()  
  h2 = DoHeatmap(seurat_object, features = cluster.specific.genes$gene, assay = params$output_plots$assay) & NoLegend()  
  markers_file_name <- paste0(results_dir,"/RNA.markers_min.pct_",params$markers$min.pct,
                              "_logfc_",params$markers$logfc.threshold,".txt")
  
  write.table(markers, markers_file_name, quote = F, row.names = F)
}

#################### generate plots
cat(paste0("================= Generating output plots \n"), 
    file = paste0(results_dir, "/log.txt") , append = TRUE)

#if (params$output_plots$generate_plots)
#{
#markers_file_name <- paste0(results_dir,"/RNA.markers_min.pct_",params$markers$min.pct,
#                            "_logfc_",params$markers$logfc.threshold,".xlsx")

#markers = read.xlsx(markers_file_name)
reduction = params$output_plots$reduction
DefaultAssay(seurat_object) <- params$output_plots$assay
Idents(seurat_object) <- "seurat_clusters"
top_markers_per_clusters <- markers %>% group_by(cluster) %>%
  slice_max(avg_log2FC,n = 5, with_ties = FALSE)

clusters <- unique(markers$cluster)
p3 <- DimPlot(seurat_object, reduction = reduction, label = T, label.size = 5)  
list_plots = list()
for (i in seq(length(clusters))) 
{
  list_plots[[i]] = FeaturePlot(seurat_object, reduction = reduction, 
                                features = top_markers_per_clusters$gene
                                [top_markers_per_clusters$cluster==clusters[i]] ,
                                label = F, min.cutoff = 'q5' , ncol = 3, max.cutoff = "q95" )+
    plot_annotation(title = paste0("Cluster_",clusters[i]),
                    theme = theme(plot.title = element_text(size = 14 , color = "red"),
                                  plot.title.position = "plot") ) + p3 
  
}

## plot genes of interest
genesOfInterest <- scan(params$output_plots$genesOfinterest_file, what = "character")
genesOfInterest_plots = list()
if (length(genesOfInterest) != 0)
{
  #genesOfInterest = sort(genesOfInterest)
  for (i in seq(length(genesOfInterest))) 
  {
    if (genesOfInterest[i] %in% rownames(seurat_object)){
      genesOfInterest_plots[[i]] = FeaturePlot(seurat_object, reduction = reduction, 
                                               features = genesOfInterest[i] ,
                                               label = F, min.cutoff = 'q5',max.cutoff = "q95",
                                               pt.size = params$output_plots$pointSize)
    }
  }
}
genesOfInterest_plots[[length(genesOfInterest)+1]] <- VlnPlot(
  seurat_object, features = genesOfInterest,ncol =  3)

#Write all plots to PDF file
pdf(paste0(results_dir , "all_plots.pdf") , onefile= TRUE, width = 10 )

#pca_heatmap
elbow_plot
if(params$clustering$resolution == "auto"){
  num_clusters_plot
}

clusters_umap_plot
if("group" %in% colnames(seurat_object@meta.data)){
  sample_umap_plot 
}

h1+h2
list_plots #top_marker_umaps

if (length(genesOfInterest) != 0)
{
  genesOfInterest_plots
}

dev.off()

#}




