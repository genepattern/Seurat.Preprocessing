# Seurat.Preprocessing (v1)
---
**Description**: GenePattern module which implements the preprocessing steps for Seurat. You may need to run this module multiple times if you want to change the filtering step.

**Author**: Edwin Juárez

**Contact**: [Forum Link](https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/genepattern-help)

**Algorithm Version**: Seurat 3.2.0

---

## Summary
---

The `Seurat.Preprocessing` Module aims to provide a way to integrate the multiple stemps in the preprocessing of single-cell datasets. The resultant preprocessed dataset  can then be used for downstream analyses and visualizations (such as Seurat.Clustering).


## References
---
[Satija Lab](https://satijalab.org)

[Seurat](https://satijalab.org/seurat/)

### Technical Notes


## Parameters
---

| Name | Description |
-------|--------------
| tenx_data_dir         | `.tar.gz` or `.zip` file input that contains the  raw single cell data -- currently only 10x data is supported.|
| column_name            | 	column name of percent mitochondrial genes. Note: not all datasets have this column, those who do often times name it percent.mt].|
| pattern        | 	what pattern to use to label mitochondrial genes [often times this is annotated as MT-].|
| file_name      | 	Basename of the file to be saved.|
|**QC Parameters**||
| first_feature  | [For QC plots] First feature to plot as a violin plot [typically one of the columns in the matrix.mtx file. Sometimes named features.tsv].|
|second_feature|	[For QC plots] Second feature to plot as a violin plot [typically one of the columns in the matrix.mtx file. Sometimes named features.tsv].|
|third_feature|	[For QC plots] Third feature to plot as a violin plot [typically one of the columns in the matrix.mtx file. Sometimes named features.tsv]. Leave blank if you don't want this to be plotted.|
|min_n_features| [For filtering] Min number of Genes that need to be expressed in a cell to be included|
|max_n_features| [For filtering] Max number of Genes that can to be expressed in a cell to be included.|
|max_percent_mitochondrial|[For filtering] Maximum percent of Genes that are labeled as mitochondrial a cell to be included.|
| **Normalization & Dimension Reduction parameters** | advanced parameters|
|norm_method|	Method for normalization.|
|scale_factor|	Scaling to be applied after normalization.|
|feat_sel_method|Method for feature selection. You should probably not change this unless you really know what you are doing.|
|num_features|	Number of top features color during feature selection.|
|num_to_label|Number of top features to label.|
|vdl_num_dims|Number of PCA dimensions to visualize.|
|vdhm_num_dims|	Number of dimensions for the dimensional reduction heatmap.|
|cells|Number of top cells to plot.|


## Output Files
---

1. `<file_name>.rds` [typically `seurat_preprocessed_dataset.rds`]
    - The `.rds` file can be used on another one of GenePattern's Seurat suite modules, such as the `Seurat.Clustering` module.
2. `<your_output_file_name>.pdf` [typically `Rplots.pdf`]
    - The `.pdf` file which contains plots of the preprocessing steps.


## License
---

`Seurat.Preprocessing` is distributed under a modified BSD license available at https://github.com/genepattern/Seurat.Preprocessing/blob/develop/LICENSE


## Platform Dependencies
---

| Task Type | CPU Type | Operating System | Language |
------------|----------|------------------|----------|
|           |  Any     | Any              | R 4.0.2  |


## Version Comments
---

| Version | Release Date | Description                                 |
----------|--------------|---------------------------------------------|
| 1       | 2020-11-16          | Initial Release of `Seurat.Preprocessing` |
