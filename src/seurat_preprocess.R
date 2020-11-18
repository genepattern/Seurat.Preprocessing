# Requirements:
# - Seurat 3.1.5
# - scater

# Seurat wrapper
# Imputs: 1- matrix.mtx
#         2- genes.tsv
#         3- barcodes.tsv

################################################################################
#Function Definitions
################################################################################

## Load the RDS file

load_rds <- function(input.file){
  pbmc=NULL

  if (file.exists(input.file)){
  	pbmc = readRDS(input.file)
  } else {
  Print('Input file could not be found!')
  }
  return(pbmc)
}

### Load 10x data function
setupR <- function(tenx_data_dir){
    write("Loading libraries...", stdout())
    suppressMessages(library(Seurat))
    suppressMessages(library(scater))
    fig_height=450
    fig_width=800
    # Load the PBMC dataset
    write("Loading the dataset...", stdout())

    ## READ DATA HERE
    ## FIRST UNZIP THE FILE
    # UNZIP TBI
    print("About to read")
    print(tenx_data_dir)
    if (grepl('http', tenx_data_dir, fixed=TRUE)){
      download.file(tenx_data_dir, paste('/temp/',basename(tenx_data_dir),sep=''))
      print("File downloaded to ")
      print(paste('/temp/',basename(tenx_data_dir),sep=''))
      print(list.files('/temp/'))
      if (grepl('.tar', tenx_data_dir, fixed=TRUE)){
          print('Untarring')
          untar(paste('/temp/',basename(tenx_data_dir),sep=''),exdir='/temp/10xdata/')
        }else if(grepl('.zio', tenx_data_dir, fixed=TRUE)){
          print('Unzipping')
          unzip(paste('/temp/',basename(tenx_data_dir),sep=''),exdir='/temp/10xdata/')
        }
      print('File extracted to /temp/10xdata/')
      print(list.files('/temp/'))
      print(list.files('/temp/10xdata/'))
    } else {
      if (grepl('.tar', tenx_data_dir, fixed=TRUE)){
          print('Untarring')
          untar(tenx_data_dir, extir='/temp/10xdata/')
        }else if(grepl('.zio', tenx_data_dir, fixed=TRUE)){
          print('Unzipping')
          unzip(tenx_data_dir, extir='/temp/10xdata/')
        }
    }

    #suppressMessages(pbmc.data <- Read10X(data.dir = "data/pbmc3k/filtered_gene_bc_matrices/hg19/"))
    suppressMessages(pbmc.data <- Read10X(data.dir = '/temp/10xdata/filtered_gene_bc_matrices/hg19/'))

    ## READ SPARSE COUNTS, IGNRED FOR NOW
    #raw_counts <- readSparseCounts(file="https://datasets.genepattern.org/data/module_support_files/Conos/HNSCC_noribo.txt")
    #hnscc <- CreateSeuratObject(counts = raw_counts, project = "HNSCC")

    # Initialize the Seurat object with the raw (non-normalized data).
    suppressMessages(pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200))
    write("Done", stdout())
    #     return(hnscc)
    return(pbmc)
}

# mitochondrial QC metrics
set_mito_qc <- function(colName, pat) {

    if (colName == 'PARAMETER_LEFT_INTENTIONALLY_BLANK'){
        print("set_mito_qc won't be called because colName is empty")
        print(colName)
    } else if (pat == 'PARAMETER_LEFT_INTENTIONALLY_BLANK'){
        print("set_mito_qc won't be called because pat is empty")
        print(pat)
    } else {
        write("Calculating the frequency of mitochondrial genes...", stdout())
        pattern <- paste("^", trimws(pat, which = "both"), sep="")

        # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
        pbmc[[colName]] <- PercentageFeatureSet(pbmc, pattern = pattern)
        write("Done!", stdout())
    }

    return(pbmc)
}

#### Violin Plots
tripleViolin <- function(first, second, third){

    feats <- c(first, second, third)
    plot(VlnPlot(pbmc, features = feats, ncol = 3, combine=TRUE), fig.height=5, fig.width=15)
    return("")
}

# FILTER DATA

my_subset <- function(min_n_features, max_n_features, max_percent_mitochondrial){
    write("About to filter data...", stdout())
    pbmc <- subset(pbmc, subset = nFeature_RNA > min_n_features & nFeature_RNA < max_n_features & percent.mt < max_percent_mitochondrial)
    write('filtering done!', stdout())
    return(pbmc)
}

# NORMALIZE DATA

norm_pbmc <- function(meth, scale){
    write("Normalizing data...", stdout())
    invisible(pbmc <- NormalizeData(pbmc, normalization.method = meth, scale.factor = scale, verbose = F))
    write('Normalization done!', stdout())
    return(pbmc)
}

# FEATURE SELECTION

feat_sel_plot <- function(meth, nFeat, nLabel){
    write("Identifying variable features...", stdout())
    invisible(capture.output(pbmc <- FindVariableFeatures(pbmc, selection.method = meth, nfeatures = nFeat,
                                                         verbose=F)))
    write("Done!", stdout())

    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(pbmc), nLabel)

    # plot variable features with and without labels
    invisible(capture.output(plot1 <- VariableFeaturePlot(pbmc)))
    invisible(capture.output(plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)))
    print(plot2)
    #plot(CombinePlots(plots = list(plot1, plot2)))
    return(pbmc)
}

## SCALE DATA

myscale <- function(pbmc){
    write("Scaling data...", stdout())
    all.genes <- rownames(pbmc)
    invisible(capture.output(pbmc <- ScaleData(pbmc, features = all.genes, verbose = F)))
    write('done!', stdout())
    return(pbmc)
}

## PCA

mypca <-function(pbmc){
    feats <- VariableFeatures(object = pbmc, verbose = F)
    pbmc <-RunPCA(pbmc, features = feats, nfeatures.print=5, verbose = F)

    return(pbmc)
}

## VISUALIZE PCA

vdp <- function(p1){

    plot(DimPlot(pbmc, reduction = "pca"))
    return("")
}

## More plots
vdl <- function(nDims){
    dim_range = c(1,strtoi(nDims))
    print(VizDimLoadings(pbmc, dims = dim_range, reduction = "pca"))
    return("")
}


## ELBOW PLOT

ebp <- function(){
    plot(ElbowPlot(pbmc))
    return(pbmc)
}

## HEAT MAP
vdhm <- function(nd,c){
    if (nd == 1){
        dim_range = 1
    } else {
        dim_range = c(1:strtoi(nd))
    }

    print(DimHeatmap(pbmc, dims = dim_range, cells = c, balanced = TRUE))
    return(pbmc)
}

save_it <- function(fileName){
    saveRDS(pbmc, file = fileName)
    print("Saved file!")
    return(pbmc)
}

################################################################################
#Parse Parameters
################################################################################
print('==========================================================')
print("Loading library: optparse")
library("optparse")

# Parse input arguments
parser = OptionParser()
# parameter types: 'character', 'integer', 'logical', 'double', or 'complex'
# ====================================
# Paramter for setupR
parser <- add_option(parser, c("--tenx_data_dir"), help = "List of files to load.")
parser <- add_option(parser, c("--input_rds"), help = "RDS file created by Seurat.QC.")
# ====================================
# PARAMETERS for set_mito_qc
parser <- add_option(parser, c("--column_name"), type='character', default='PARAMETER_LEFT_INTENTIONALLY_BLANK', help = "column name of percent mitochondrial genes [often times it's called percent.mt].")
parser <- add_option(parser, c("--pattern"),type='character', default='PARAMETER_LEFT_INTENTIONALLY_BLANK', help = "what pattern to use to label mitochondrial genes [often times this is annotated as MT-].")
# ====================================
# Parameters for the violin plot, tripleViolin
parser <- add_option(parser, c("--first_feature"),type='character',default='nFeature_RNA', help = "First feature to plot as a violin plot.")
parser <- add_option(parser, c("--second_feature"),type='character', default='nCount_RNA', help = "Second feature to plot as a violin plot.")
parser <- add_option(parser, c("--third_feature"),type='character', default='percent.mt', help = "Third feature to plot as a violin plot.")
# ====================================
# Parameters for filter data, my_subset
parser <- add_option(parser, c("--min_n_features"),type='integer', help = "Min number of Genes that need to be expressed in a cell to be included.")
parser <- add_option(parser, c("--max_n_features"),type='integer', help = "Max number of Genes that can to be expressed in a cell to be included.")
parser <- add_option(parser, c("--max_percent_mitochondrial"), type='integer', default=5, help = "Maximum percent of Genes that are labeled as mitochondrial a cell to be included.")
# ====================================
# Parameters for normalize data, norm_pbmc
parser <- add_option(parser, c("--norm_method"),type='character',default='LogNormalize', help = "Method for normalization.")
parser <- add_option(parser, c("--scale_factor"),type='integer',default=10000, help = "Scaling to be applied after normalization.")
# ====================================
# Parameters for Feature Selection, feat_sel_plot
parser <- add_option(parser, c("--feat_sel_method"),type='character',default='vst', help = "Method for feature selection. You should probably not change this unless you really know what you are doing.")
parser <- add_option(parser, c("--num_features"),type='integer',default=2000, help = "Number of top features color during feature selection.")
parser <- add_option(parser, c("--num_to_label"),type='integer',default=10, help = "Number of top features to label.")
# ====================================
# Parameter for Vizualize Dimention Loadings, vdl
parser <- add_option(parser, c("--vdl_num_dims"),type='integer',default=2, help = "Number of PCA dimensions to visualize.")
# ====================================
# Parameters for Heat Map, vdhm
parser <- add_option(parser, c("--vdhm_num_dims"),type='integer',default=15, help = "Number of dimensions for the dimensional reduction heatmap.")
parser <- add_option(parser, c("--cells"),type='integer',default=500, help = "Number of top cells to plot.")
# ====================================
#parameter for save_it
parser <- add_option(parser, c("--file_name"),type='character',default='seurat_preprocessed_dataset', help = "Basename of the file to be saved.")
# ====================================


print('==========================================================')
args <- parse_args(parser)
print('Parameters used:')
print(args)
print('==========================================================')


################################################################################
#Begin Running the functions
################################################################################


#-------------------------------------------------------------------------------
# Ignore these calls, these were here when this script did it all in one shot
# you should call Seurat.QC first
#-------------------------------------------------------------------------------
# Call the setupR function
#suppressMessages(pbmc <- setupR(args$tenx_data_dir))
# call the set_mito_qc function
#suppressMessages(pbmc <- set_mito_qc(args$column_name, args$pattern))
#tripleViolin(args$first_feature, args$second_feature, args$third_feature)
#-------------------------------------------------------------------------------
# End of ignore
#-------------------------------------------------------------------------------

pbmc <- load_rds(args$input_rds)

pbmc <- my_subset(args$min_n_features, args$max_n_features, args$max_percent_mitochondrial)

pbmc <- norm_pbmc(args$norm_method, args$scale_factor)

pbmc <- feat_sel_plot(args$feat_sel_method, args$num_features, args$num_to_label)

pbmc <- myscale(pbmc)

write("Performing PCA...", stdout())

# We should add Paramters for PCA
pbmc <- mypca(pbmc)
write("Done!", stdout())

vdp()

vdl(args$vdl_num_dims)

ebp()

vdhm(args$vdhm_num_dims, args$cells)

save_it(paste(args$file_name,'.rds',sep=''))

#sessionInfo()
