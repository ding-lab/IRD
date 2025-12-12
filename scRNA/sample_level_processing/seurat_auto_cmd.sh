
sample=BM2
output_path=/diskmnt/Projects/myeloma_scRNA_analysis/MMY_IRD/${sample}/Seurat
tag=noPollock

Rscript /diskmnt/Projects/Users/jwang/scripts/IRD_Seurat.R \
-i /diskmnt/Projects/myeloma_scRNA_analysis/MMY_IRD/${sample}/BAM_cellranger_hg38/${sample}/outs/raw_feature_bc_matrix/ \
-s ${sample} \
-o ${output_path} \
-t ${tag} \
--scrublet /diskmnt/Projects/myeloma_scRNA_analysis/MMY_IRD/${sample}/Seurat/ \
--mito_max 20 \
--ncount_min 1000 \
--ncount_max 15000 \
--nfeature_min 200 \
--nfeature_max 10000 \

