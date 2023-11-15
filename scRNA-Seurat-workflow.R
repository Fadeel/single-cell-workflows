#!/usr/bin/env Rscript


library(Seurat)
library(yaml)
library(stringr)
library(ggplot2)
library(plyr)
library(dplyr)
#library(openxlsx)
library(RColorBrewer)
library(patchwork)
#library(clusterProfiler)
#library("org.Hs.eg.db")
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

# read parameters YAML file
params = read_yaml(parameters_file)

# create results folder
results_dir = paste0(params$results_dir ,params$project_name, "/")

system(sprintf("mkdir -p %s",results_dir ))

#copy parameters file to results folder
system(sprintf("cp  %s %s", parameters_file , results_dir ))

#################### read input 

# add info to log file
cat(paste0("================= reading input \n"), 
    file = paste0(results_dir, "/log.txt"))

if(params$input_data == "cell_ranger"){
  data <- Read10X(data.dir = params$input_location,
                  strip.suffix = params$cell_ranger$strip.suffix)
  
  # create Seurat Object
  seurat_object = CreateSeuratObject(counts = data, 
                                     min.cells = params$cell_ranger$min.cells)
  
  # compute percent of mitochondrial reads
  if (params$species == "human")
  {
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  }else{
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
  }
  
  #generate QC plots
  QC_plots <- list()
  QC_plots[[1]] <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                           ncol = 3 , log = T , pt.size = 0)
  
  QC_plots[[2]] <-FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  QC_plots[[3]] <-FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  
  cat(paste0("number of cell barcodes = ", ncol(seurat_object), "\n"), 
      file = paste0(results_dir, "log.txt"), append = TRUE)
  cat(paste0("number of genes = ", nrow(seurat_object[["RNA"]]) , "\n") ,
      file = paste0(results_dir, "log.txt"), 
      append = TRUE)

    # write Seurat object
    #saveRDS(seurat_object, paste0(results_dir, params$project_name, "_raw.rds" ))
  
}else{
  seurat_object <- readRDS(params$input_location)
}

#################### filtering

cat(paste0("================= filtering \n"), 
    file = paste0(results_dir, "/log.txt"), append = TRUE)

if (params$filtering$do_filtering)
{
  
  if(params$input_data == "cell_ranger")
  {
    seurat_object  <- subset(
      x = seurat_object,
      subset = 
        nCount_RNA >  params$filtering$nCount_RNA_min &
        nCount_RNA < params$filtering$nCount_RNA_max &
        nFeature_RNA > params$filtering$nFeature_RNA_min &
        nFeature_RNA < params$filtering$nFeature_RNA_max &
        percent.mt < params$filtering$percent.mt)
      
      #generate QC plots
      QC_plots[[4]] <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                               ncol = 3 , log = T , pt.size = 0)
      
      QC_plots[[5]] <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
      QC_plots[[6]] <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
      
      # add info to log file
      cat(paste0("number of cell barcodes after filtering= ", ncol(seurat_object), "\n"), 
          file = paste0(results_dir, "/log.txt"), append = TRUE)
      
      cat(paste0("number of genes after filtering = ", nrow(seurat_object[["RNA"]]) , "\n") ,
          file = paste0(results_dir, "/log.txt"), 
          append = TRUE)
      #saveRDS(seurat_object, paste0(results_dir, params$project_name, "_filtered.rds" ))
      
      
  }else{
    
    ## re-analysis using certain group of cells 
    if(params$subsetting$value != "all"){
      Idents(seurat_object) <- params$subsetting$cell_label
      values = strsplit(params$subsetting$value , split = ",")[[1]]
      seurat_object  <- subset(x = seurat_object, idents = values  )
      #saveRDS(seurat_object, paste0(results_dir, params$project_name, "_filtered.rds" ))
      } 
  }
  
}



#################### cell cycle analysis

cat(paste0("================= cell cycle analysis \n"), 
    file = paste0(results_dir, "/log.txt") , append = TRUE)

DefaultAssay(seurat_object) = "RNA"
seurat_object = NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst")
seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))

if (params$cell_cycle$do_analysis)
{
  if(params$species == "mouse"){
    mouse_cc  = readRDS("mouse_CC.genes.rds")
    s.genes = mouse_cc$s.genes
    g2m.genes = mouse_cc$g2m.genes
  }else{
    s.genes = cc.genes$s.genes
    g2m.genes = cc.genes$g2m.genes
  }
  
  seurat_object <- CellCycleScoring(seurat_object, s.features = s.genes, 
                                    g2m.features = g2m.genes, set.ident = TRUE)
  seurat_object = RunPCA(seurat_object , features = c(s.genes, g2m.genes),
                         reduction.name = "pca.cc")
  CC_pca_plot <- DimPlot(seurat_object, reduction = "pca.cc" , pt.size =0.5)
}



#################### Normalization

if(! ("percent.mt" %in% colnames(seurat_object@meta.data))){
  if (params$species == "human")
    {
      seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-", assay = "RNA")
    }else{
      seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-", assay = "RNA")
    }
  
  
}

cat(paste0("================= Normalization \n"), 
    file = paste0(results_dir, "/log.txt") , append = TRUE)

if(params$normalization$do_normalization){
  if(params$cell_cycle$do_regression == "diff"){
    # regressing diff. between S.score and G2M.score
    seurat_object$CC.Difference = seurat_object$S.Score - seurat_object$G2M.Score
  
    # run SCTransform
    seurat_object <- SCTransform(seurat_object, vars.to.regress = 
                                   c("percent.mt","CC.Difference"), verbose = FALSE)
    
    seurat_object = RunPCA(seurat_object , features = c(s.genes, g2m.genes),
                           reduction.name = "pca.cc2")
    
    CC_regressed_pca_plot <- DimPlot(seurat_object, reduction = "pca.cc2" , pt.size = 0.5)
    
  }else if(params$cell_cycle$do_regression == "all"){
    
    seurat_object <- SCTransform(seurat_object, vars.to.regress = 
                                   c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)
    
    seurat_object = RunPCA(seurat_object , features = c(s.genes, g2m.genes),
                           reduction.name = "pca.cc2")
    CC_regressed_pca_plot <- DimPlot(seurat_object, reduction = "pca.cc2" , pt.size = 0.5)
  }else{
    
    if(params$normalization$regress_mito)
    {
        seurat_object <- SCTransform(seurat_object, vars.to.regress = 
                                       c("percent.mt"), verbose = FALSE)
    }else{
      seurat_object <- SCTransform(seurat_object, verbose = FALSE)
    }

  }
  
  ## top genes
  top = head(VariableFeatures(seurat_object , assay = "SCT"), 10)
  plot1 <- VariableFeaturePlot(seurat_object, assay = "SCT" )
  top_var_plot <- LabelPoints(plot = plot1, points = top, repel = TRUE)
}

#################### Dimensionality reduction
cat(paste0("================= Dimensionality reduction using PCA \n"), 
    file = paste0(results_dir, "/log.txt") , append = TRUE)

if(params$dim_reduction$do_PCA){
  seurat_object <- RunPCA(seurat_object, reduction.name = "pca" ,npcs = 50, assay = params$dim_reduction$assay)
  
  pca_heatmap <- DimHeatmap(seurat_object, dims = 1:8, cells = 200, balanced = TRUE , ncol = 4,
                            assays = params$dim_reduction$assay,fast = FALSE)
  elbow_plot <- ElbowPlot(seurat_object , ndims = 50)
}

#################### clustering
cat(paste0("================= Clustering \n"), 
    file = paste0(results_dir, "/log.txt") , append = TRUE)

if(params$clustering$do_clustering){
  
seurat_object <- FindNeighbors(seurat_object, dims = 1:params$dim_reduction$n_dims, k.param = params$clustering$knn)

  if(params$clustering$resolution == "auto"){
    number.clusters = pickNumClusters(seurat_object, step = 0.02, end = 1, graph = paste0(params$dim_reduction$assay,"_snn") ,algorithm_num = 1)
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

  seurat_object <- FindClusters(seurat_object,resolution = res, graph.name = paste0(params$dim_reduction$assay,"_snn"), algorithm = 1)
  
  cat(paste0("number of cells per cluster : \n"), 
      file = paste0(results_dir, "log.txt") , append = TRUE)
  cat( table(seurat_object$seurat_clusters) , 
      file = paste0(results_dir, "log.txt") , append = TRUE)
}

#################### UMAP plot
cat(paste0("\n================= Dimensionality reduction using UMAP \n"), 
    file = paste0(results_dir, "/log.txt") , append = TRUE)

if(params$UMAP$do_UMAP){
  seurat_object = RunUMAP(seurat_object, dims = 1:params$dim_reduction$n_dims, n.neighbors = params$clustering$knn)
  
  clusters_umap_plot <- DimPlot(seurat_object, reduction = "umap" , label = T, 
                                group.by = "seurat_clusters", label.size = 5)
  
  if("Phase" %in% colnames(seurat_object@meta.data)){
    cc_umap_plot <- DimPlot(seurat_object, reduction = "umap" , label = F, group.by = "Phase")
  }
  
  if("sample" %in% colnames(seurat_object@meta.data)){
    sample_umap_plot <- DimPlot(seurat_object, reduction = "umap" , label = F, group.by = "sample")
  }
}

saveRDS(seurat_object, paste0(results_dir,params$project_name , ".rds" ))


#################### cluster  markers
cat(paste0("================= Identifying cluster markers \n"), 
    file = paste0(results_dir, "/log.txt") , append = TRUE)

if(params$markers$find_markers){
  markers = FindAllMarkers(seurat_object, assay = "RNA",  min.pct = params$markers$min.pct, 
                           logfc.threshold = params$markers$logfc.threshold , only.pos = T )
  markers = markers %>% filter(p_val_adj < params$markers$p_val_adj)
  markers = markers %>% mutate(diff.pct = pct.1 - pct.2)
  
  t = params$markers$top_markers
  cluster.specific.genes = markers %>% group_by(cluster) %>%
    slice_max(diff.pct,n = t, with_ties = FALSE)
  
  top_markers_per_clusters <- markers %>% group_by(cluster) %>%
    slice_max(avg_log2FC,n = t, with_ties = FALSE)
  
  h1 = DoHeatmap(seurat_object, features = top_markers_per_clusters$gene, assay = "RNA") & NoLegend()  
  h2 = DoHeatmap(seurat_object, features = cluster.specific.genes$gene, assay = "RNA") & NoLegend()  
  markers_file_name <- paste0(results_dir,"/RNA.markers_min.pct_",params$markers$min.pct,
                              "_logfc_",params$markers$logfc.threshold,".txt")
  
  write.table(markers, markers_file_name, quote = F, row.names = F)
}

if(params$markers$GO_enrich){
  
  clusters <- unique(markers$cluster)
  list_GO_plots <- list()
  for(i in seq(length(clusters)))
  {
  genes = markers$gene[markers$cluster == clusters[i]]

  if (length(genes) > 0)
  {
    #ENSEMBL_ids <- mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column="ENSEMBL")
    
    # go_enrich <- enrichGO(gene = ENSEMBL_ids,
    #                       OrgDb = "org.Hs.eg.db", 
    #                       keyType = 'ENSEMBL',
    #                       readable = T,
    #                       ont = "BP",
    #                       pvalueCutoff = 0.2, 
    #                       qvalueCutoff = 0.5)
    
    
    #write.xlsx(go_enrich@result, file = paste0(results_dir, "/cluster_",clusters[i],"_GO.xlsx"))
    
    ## use simplify to remove redundant terms
    #list_GO_plots[[i]] <- barplot(go_enrich,drop = FALSE,showCategory = 10, 
    #                              title = "GO Biological Pathways",font.size = 10)
  }
  }
  #dotplot(go_enrich, title = "Enriched Pathways", showCategory = 15) 
}


#################### generate plots
cat(paste0("================= Generating output plots \n"), 
    file = paste0(results_dir, "/log.txt") , append = TRUE)

#if (params$output_plots$generate_plots)
#{
  markers_file_name <- paste0(results_dir,"/RNA.markers_min.pct_",params$markers$min.pct,
                              "_logfc_",params$markers$logfc.threshold,".xlsx")
  
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
                                  label = F, min.cutoff = 'q1' , ncol = 3 )+
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
                                                 label = F, min.cutoff = 'q5',max.cutoff = 'q95',
                                                 pt.size = params$output_plots$pointSize)
      }
    }
  }
  genesOfInterest_plots[[length(genesOfInterest)+1]] <- VlnPlot(
    seurat_object, features = genesOfInterest,ncol =  3)

  #Write all plots to PDF file
  pdf(paste0(results_dir,params$project_name  , "_plots.pdf") , onefile= TRUE, width = 10 )
  
  if(params$input_data == "cell_ranger")
  {
    QC_plots
  }
  
  if(params$cell_cycle$do_analysis)
  {
    CC_pca_plot
  }
  
  if(params$cell_cycle$do_regression %in% c("diff","all"))
  {
    CC_regressed_pca_plot
  }
  
  if(params$normalization$do_normalization)
  {
    top_var_plot
  }
  
  pca_heatmap
  elbow_plot
  if(params$clustering$resolution == "auto"){
    num_clusters_plot
  }
  clusters_umap_plot
  if("Phase" %in% colnames(seurat_object@meta.data)){
    cc_umap_plot  
  }
  
  if("sample" %in% colnames(seurat_object@meta.data)){
    sample_umap_plot 
  }
  
  h1+h2
  list_plots #top_marker_umaps
  if(params$markers$GO_enrich)
    {
    print(length(list_GO_plots))
    list_GO_plots
    }
  
  if (length(genesOfInterest) != 0)
  {
  genesOfInterest_plots
  }
  #cluster_annot_umap
  
  dev.off()

#}




