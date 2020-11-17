docker run -v $PWD:/LOCAL -w $PWD/Job_1 -t genepattern/seurat-suite:2.4 Rscript --no-save --quiet --slave --no-restore  /LOCAL/seurat_preprocess.R\
 --tenx_data_dir 'https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz' \
 --column_name "percent.mt" --pattern 'MT-'\
 --first_feature 'nFeature_RNA' --second_feature 'nCount_RNA' --third_feature 'percent.mt'\
 --min_n_features 200 --max_n_features 2500 --max_percent_mitochondrial 5\
 --norm_method 'LogNormalize' --scale_factor 10000\
 --num_features 2000 --num_to_label 10\
 --vdl_num_dims 2\
 --vdhm_num_dims 15 --cells 500\
 --file_name "test_run"
