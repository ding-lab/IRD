########################################################
# run using cellchat environment
########################################################
#library(harmony)
library(RColorBrewer)
library(Seurat)
library(ggplot2)
library(reshape2)
library(viridis)
library(dplyr)
library(Matrix)
library(cowplot)
library(getopt)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(scProportionTest)
library(enrichR)
library(ggpubr)
library(DESeq2)
library(rstatix)
options(future.globals.maxSize = 500 * 1024^3) 

#---------Gene Annotations---------------------------------------
gs <- readRDS('/diskmnt/Projects/myeloma_scRNA_analysis/MMY_IRD/analysis/resources/gene_sets.rds')
annotations <- gs$annotations
s.genes <- gs$s.genes
g2m.genes <- gs$g2m.genes
ig_genes <- gs$ig_genes #length=436
hemogenes <- gs$hemo_genes
protein_coding_genes <- gs$protein_coding_genes
#--- palettes
colors <- readRDS('/diskmnt/Projects/myeloma_scRNA_analysis/MMY_IRD/analysis/resources/colors.rds')
#----------------------------------------------
# function for filtering out mitochondrial, hemoglobin, Ig, ribosomal genes from a list of genes
#----------------------------------------------
filter_genes <- function(genes){
  de_fi1 <- genes[!genes%in% ig_genes & !genes %in% hemogenes]
  de_fim1 <- de_fi1[!str_starts(de_fi1, 'MT-')]
  de_fimr1 <- de_fim1[!str_starts(de_fim1, 'RPL') & !str_starts(de_fim1, 'RPS')]
  genesf <- de_fimr1
  return(genesf)
}
# return only protein-coding genes (excl. mito, hemo, Ig, ribo)
filter_genes_p <- function(genes){
de_fi1 <- genes[!genes%in% ig_genes & !genes %in% hemogenes & genes %in% protein_coding_genes]
  de_fim1 <- de_fi1[!str_starts(de_fi1, 'MT-')]
  de_fimr1 <- de_fim1[!str_starts(de_fim1, 'RPL') & !str_starts(de_fim1, 'RPS')]
  genesf <- de_fimr1
  return(genesf)
}

#-----------------------
setwd("/diskmnt/Projects/myeloma_scRNA_analysis/MMY_IRD/analysis/MMRF_longitudinal")

#483 sample object in nat can submission
mmrf <- readRDS("/diskmnt/Projects/MMRF_analysis/Integration/20240516_All_CENTER_B4/Integration/All_center_B1_B4/All_compartments/SeuratObj_in_483_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_v1.rds")
obj <- mmrf

mmrf_meta <- readRDS("/diskmnt/Projects/MMRF_analysis/Integration/20240516_All_CENTER_B4/Integration/All_center_B1_B4/All_compartments/Meta_data_SeuratObj_in_483_samples_Human_Ref_SM_CB_LogNorm_PC25_Harmony_v2.rds")

# filter for longitudinal pairs with remission/PT collection_event plus one other timepoint
clin_pt <- mmrf_meta %>% filter(
    d_amm_tx_asct_1st=="yes",
    collection_event %in% c( "100 Post ASCT" , "Post-transplant", "Remission/Response"),
)
pt_id <- unique(clin_pt$public_id)

# make histogram of study day of visit for post-ASCT / remission timepoints
visitdays <- data.frame(value=as.integer(unique(clin_pt$VISITDY)))
# Plot histogram
p <- ggplot(visitdays, aes(x = value)) +
  geom_histogram(binwidth = 50, fill = "black", color = "white") +
  labs( x = "Study Day of Visit for PT/Remission Collection", y = "Count") + scale_y_continuous(limits=c(0,8), breaks=seq(0,8,by=2)) +
  theme_minimal()
pdf('histogram_days_after_dx_mmrfLongitudinal.pdf', width=4, height=4)
print(p)
dev.off()

# get all samples from the cases that have a post-ASCT/ remission timepoint
clinf <- mmrf_meta %>% filter(
    public_id %in% pt_id,
    collection_event %in% c( "100 Post ASCT" , "Post-transplant", "Remission/Response", "Baseline", "Relapse/Progression")
) #unique(clinf$sample_id) is 97, public_id is 38
rownames(clinf) <- clinf$barcode_full

objf <- subset(obj, subset = sample_id %in% unique(clinf$sample_id) )

# add meta data
meta_add <- clinf[rownames(objf@meta.data), c('public_id', 'collection_event')]
objf <- AddMetaData(objf, meta_add)

# get levels for ordering mmrf celltypes
ctlevels <- read.table('/diskmnt/Projects/myeloma_scRNA_analysis/MMY_IRD/analysis/MMRF_longitudinal/mmrf_obj_celltype_levels.txt', sep='\t', header=T, quote='', comment.char='')



#---------------
# see how cell type abundance changes over time in non-PC fraction
#----------------
pcs <- unique(ctlevels$ct[grepl("^Pc_", ctlevels$ct)])
bs <- unique(ctlevels$ct[grepl("^B_", ctlevels$ct)])
tnks <- ctlevels$ct[29:58]
cd8s <- ctlevels$ct[40:53]
cd4s <- ctlevels$ct[29:39]

# non PC only
tme <- objc_meta %>% filter(!celltype_f %in% pcs)

# Count number of cells per cell type per collection_event per public ID
tme$celltype_f <- droplevels(tme$celltype_f)
tme$sample_id <- paste0(tme$public_id,"|", tme$collection_event)
st <- as.data.frame(table(tme$sample_id))
names(st) <- c("Sample", "nCells_Sample")
sca <- as.data.frame(table(tme$sample_id, tme$celltype_f))
colnames(sca) <- c('Sample', 'CT', 'Freq')
scas <- merge(sca, st, by='Sample')
scas$Prop_of_Sample <- scas$Freq/scas$nCells_Sample*100
scas <- scas %>%
  separate(Sample, into = c("public_id", "collection_event"), sep = "\\|")

scas$CT <- factor(scas$CT, levels=ctlevels$ct)
scas$collection_event <- factor(scas$collection_event, levels=c('Baseline', 'Post-transplant', 'Remission/Response', 'Relapse/Progression' ))

# Perform statistical test and adjust p-values
my_comparisons <- list(c('Baseline', 'Remission/Response'))
stat.test <- scas %>% group_by(CT) %>% 
    wilcox_test(Prop_of_Sample ~ collection_event, comparisons = my_comparisons, p.adjust.method = "fdr") %>%
    add_significance("p.adj")

# Filter for desired CTs
selected_CTs <- bs
scac <- scas %>% filter(CT %in% selected_CTs)
stat.test.filtered <- stat.test %>% filter(CT %in% selected_CTs)
stat.test.filtered <- stat.test.filtered %>% add_xy_position(x='collection_event')
# Plot all at once with faceting
p <- ggplot(scac, aes(x=collection_event, y=Prop_of_Sample, group=public_id)) +  geom_line(color = "gray", size = 0.25) + 
  geom_point(size = 2) +  stat_pvalue_manual(stat.test.filtered, label = "p.adj.signif") + theme_classic() +
  labs(x = "Collection Event", y = "% of Sample Excl. PC") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_wrap(~ CT, scales = "free_y", ncol=4)

pdf("mmrf_Bcell_proportions_overTime.pdf", width = 8, height = 8)
print(p)
dev.off()



tnk_obj <- subset(objc, subset= celltype_f %in% tnks)


#-------------
# get cd4/cd8 ratio in IRD samples
m$id <- paste0(m$sample_id, "|", m$Timepoint)
sca <- as.data.frame(table(m$id, m$subset))
colnames(sca) <- c('Sample', 'CT', 'Freq')
cd4_cd8_ratio <- sca %>%
  filter(CT %in% c("CD4T", "CD8T")) %>%
  select(Sample, CT, Freq) %>%
  pivot_wider(names_from = CT, values_from = Freq) %>%
  mutate(CD4_CD8_ratio = CD4T / CD8T)
cd4_cd8_ratio <- cd4_cd8_ratio %>%
  separate(Sample, into = c("sample", "timepoint"), sep = "\\|")
cd4_cd8_ratio$CD4_CD8_ratio <- log2(cd4_cd8_ratio$CD4_CD8_ratio)
cd4_cd8_ratio <- cd4_cd8_ratio %>% filter(timepoint %in% c('NDMM', 'Post-ASCT'))
cd4_cd8_ratio$timepoint <- factor(cd4_cd8_ratio$timepoint, levels=c("NDMM", "Post-ASCT"))
# Plot all at once with faceting
p <- ggboxplot(cd4_cd8_ratio, x = "timepoint", y = "CD4_CD8_ratio", outlier.shape = NA) +
  stat_compare_means(method="wilcox.test", label='p.signif') +
  geom_jitter(aes(x = timepoint, y = CD4_CD8_ratio),
              position = position_jitter(0.2), size = 2) +
  labs(y = "Log2 CD4:CD8 T Cell Ratio") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
pdf("IRD_cd4cd8ratio_overTime.pdf", width=3, height=4)
print(p)
dev.off()

