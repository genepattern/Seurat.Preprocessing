docker run -v $PWD:/LOCAL -w /LOCAL/Job_1 -t genepattern/seurat-suite:4.0.3 Rscript --no-save --quiet --slave --no-restore  /LOCAL/src/seurat_preprocess.R\
 --input_rds '/LOCAL/data/test_run.rds' \
 --column_name "percent.mt" --pattern 'MT-'\
 --first_feature 'nFeature_RNA' --second_feature 'nCount_RNA' --third_feature 'percent.mt'\
 --min_n_features 2 --max_n_features 6000 --max_percent_mitochondrial 25\
 --norm_method 'LogNormalize' --scale_factor 10000\
 --num_features 2000 --num_to_label 10\
 --vdl_num_dims 2\
 --vdhm_num_dims 15 --cells 500\
 --file_name "test_run"\
 --keep_scale_data "TRUE"
