print("about to install pagoda2")
install.packages(c("devtools", "BiocManager"))
BiocManager::install(c("AnnotationDbi", "BiocGenerics", "GO.db", "pcaMethods"))
devtools::install_github("hms-dbmi/pagoda2")
library('pagoda2')

print("about to install conos")
devtools::install_github("hms-dbmi/conos")


#install.packages("BiocManager", repos = "http://cran.us.r-project.org")

#package="SingleCellExperiment"
#BiocManager::install(package)

#print("About to install scater")
#BiocManager::install("scater")



#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager", repos = "http://cran.us.r-project.org")

#print("about to install'SingleCellExperiment', 'SummarizedExperiment', 'DelayedArray', 'DelayedMatrixStats', 'BiocNeighbors', 'BiocSingular', 'BiocParallel', 'beachmat'")
#BiocManager::install(c('SingleCellExperiment', 'SummarizedExperiment', 'DelayedArray', 'DelayedMatrixStats', 'BiocNeighbors', 'BiocSingular', 'BiocParallel', 'beachmat'))
#print("about to install scater")
#BiocManager::install("scater", dependencies=TRUE)
