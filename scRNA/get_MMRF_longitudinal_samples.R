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

#--------------------
# run SCT on object
objf <- NormalizeData(objf, normalization.method='LogNormalize', assay='RNA')
objf <- CellCycleScoring(objf, s.features = s.genes, g2m.features = g2m.genes, assay='RNA', set.ident = TRUE)
objf$CC.Difference <- objf$S.Score - objf$G2M.Score
objf <- SCTransform(objf, assay='RNA', vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', "CC.Difference"), return.only.var.genes = F, conserve.memory=TRUE)
saveRDS(objf, "mmrf_longitudinal.rds", compress=F)

#--------------------
# label transfer IRD
ird <- readRDS("/diskmnt/Projects/myeloma_scRNA_analysis/MMY_IRD/analysis/organized/Overview/IRD_WU_ref_cleaned.rds")

set.seed(666)
query.obj <- objf
ref.obj <- ird
DefaultAssay(ref.obj) <- 'SCT'
ref.obj <- RunUMAP(ref.obj, reduction = "pca", dims=1:50, return.model = TRUE)
DefaultAssay(query.obj) <- 'SCT'
anchors <- FindTransferAnchors(
  reference = ref.obj,
  query = query.obj,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:50
)
query.obj <- MapQuery(
  anchorset = anchors,
  query = query.obj,
  reference = ref.obj,
  refdata = list(
    celltype.l1 = 'subset'
  ),
  reference.reduction = "pca",
  reduction.model = "umap"
)
colnames(query.obj@meta.data)[colnames(query.obj@meta.data) == "predicted.celltype.l1"] = "IRD_subset"
colnames(query.obj@meta.data)[colnames(query.obj@meta.data) == "predicted.celltype.l1.score"] = "IRD_subset_score"
objf <- query.obj
saveRDS(objf,"mmrf_longitudinal.rds", compress=F)

Idents(objf) <- 'IRD_subset'
objf <- RenameIdents(objf, 
"PC_Normal"='PC',
"PC_Neoplastic"='PC')
objf$IRD_subset <- Idents(objf)
objf$IRD_subset_f <- factor(objf$IRD_subset, levels=names(colors$subset))
objf$celltype_f <- factor(objf$celltype_subclusters_integration_Yizhe_v1, levels=ctlevels$ct)
saveRDS(objf,"mmrf_longitudinal.rds", compress=F)

#objf <- readRDS('mmrf_longitudinal.rds')

# remove samples that only have remission/PT and no baseline or relapse
single_timepoint_ids <- objf@meta.data %>%
  distinct(public_id, collection_event) %>%
  count(public_id) %>%
  filter(n == 1) %>%
  pull(public_id)

objf <- subset(objf, subset = !public_id %in% single_timepoint_ids)
saveRDS(objf,"mmrf_longitudinal.rds", compress=F)

#-----
# plot how IRD subset compares to MMRF celltyping
#-----
df <- as.data.frame(table(objf$IRD_subset_f, objf$celltype_f))
names(df) <- c('IRD_subset', 'MMRF_celltype', 'Freq')
p <- ggplot(subset(df, Freq>0), aes(x = MMRF_celltype, y = IRD_subset, size=Freq)) +
  geom_point(color = "black") +
  labs(title = "MMRF annotation vs. IRD label", x = "MMRF Celltype", y = "IRD Subset") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
pdf('mmrfAnnot_vs_IRDlabel.pdf', width=18)
print(p)
dev.off()

# remove doublets 
db <- ctlevels$ct[85:100] #all doublet clusters
objc <- subset(objf, subset = !celltype_f %in% db)
saveRDS(objc,"mmrf_longitudinal_noDoublets.rds", compress=F)
# objc <- readRDS("mmrf_longitudinal_noDoublets.rds")

# visualize and annotate
imm <- c(
'TNFRSF17','SLAMF7','SDC1', # PCs
"RAG1", "VPREB1", "CD19","MS4A1","CD79A", 'IGHM', 'IGHD', 'IGHA1', 'IGHG1', #B
'CD3D','CD3E','CD4','IL7R','TCF7', 'FOXP3', #CD4
'CD8A','CD8B','NKG7', #CD8
'FCGR3A','GNLY','KLRC1', 'KLRD1', #NK
'CD14','S100A9', 'SELL', #CD14
'MS4A7', 'TNFRSF1B', #CD16
'FCER1A', 'CD1C', 'CLEC9A', # DCs
'MPO','AZU1','ELANE', # neutrophils
'PF4', 'PPBP', 'THPO',#megakaryocytes
'HBB', 'HBA1', 'HBA2', 'GATA1', 'KLF1', 'CA3', # ery
'CXCL12', 'LEPR', 'KITLG' #MSC
)

p1 <- DotPlot(objc, features=imm, group.by='IRD_subset_f', assay='RNA') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf("IRDlabels_marker_dotplot.pdf", width=12, height=8)
print(p1)
dev.off()

p2 <- DotPlot(objc, features=imm, group.by='celltype_f', assay='RNA') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf("MMRFct_marker_dotplot.pdf", width=12, height=18)
print(p2)
dev.off()

#--------------
# make presence tile plot of sample availability
#--------------
objc_meta <- objc@meta.data[, c('public_id', 'collection_event', 'IRD_subset_f', 'celltype_f')]
objc_meta$collection_event[objc_meta$collection_event == "100 Post ASCT"] <- "Post-transplant"
#write.table(objc_meta, "mmrf_longitudinal_noDoublets_metadata.txt", sep="\t", quote=F)
#objc_meta <- read.table("/diskmnt/Projects/myeloma_scRNA_analysis/MMY_IRD/analysis/MMRF_longitudinal/mmrf_longitudinal_noDoublets_metadata.txt", sep="\t", header=T)
objc_meta$collection_event <- factor(objc_meta$collection_event, levels=c('Baseline', 'Post-transplant', 'Remission/Response', 'Relapse/Progression'))

presence_df <- objc_meta %>%
  distinct(public_id, collection_event) %>%
  mutate(has_data = TRUE)  
# All combinations of public_id and collection_event
all_combinations <- expand.grid(
  public_id = unique(cell_counts$public_id),
  collection_event = unique(cell_counts$collection_event)
)
# Join to presence_objc_meta to identify missing combos
plot_df <- left_join(all_combinations, presence_df, by = c("public_id", "collection_event")) %>%
  mutate(has_data = ifelse(is.na(has_data), FALSE, TRUE))

p <- ggplot(plot_df, aes(x = collection_event, y = public_id, fill = has_data)) +
  geom_tile(color = "grey80") +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "white")) +
  labs(title = "Timepoint Availability per Sample", x = "Timepoint", y = "Sample ID") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf("longitudinal_sample_availability_tileplot.pdf", width=3, height=6)
print(p)
dev.off()

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


# get samples where remission response proportion is higher than baseline, for each b subtype
higher_rr_ids <- list()
for(btype in bs){
  bcounts_sub <- scas %>% filter(CT ==btype & collection_event %in% c("Baseline", "Remission/Response"))
  n_ids_with_both_events <- sum(table(bcounts_sub$public_id) == 2)
  print(btype)
  print(n_ids_with_both_events)
  higher_rr_ids[[paste0(btype)]] <- bcounts_sub %>%
    select(public_id, collection_event, Prop_of_Sample) %>%
    pivot_wider(names_from = collection_event, values_from = Prop_of_Sample) %>%
    filter(`Remission/Response` > Baseline) %>%
    pull(public_id)
}

# Get vector lengths
lengths_vec <- sapply(higher_rr_ids, length)
# Calculate proportion by dividing by 32, the number of samples with data for both baseline and remission timepoints
props <- lengths_vec / 32
# Make data frame with rownames and a column "prop"
df_prop <- data.frame(prop = props)
# Row names are automatically the names of the list elements

pdf("proportion_of_b_subtype_where_remissionprop_higherThan_baseline.pdf", width=4)
df_prop$celltype <- factor(rownames(df_prop), levels=ctlevels$ct)
p <- ggplot(df_prop, aes(x = celltype, y = prop)) +
  geom_col(fill = "steelblue") +       # bar plot
  ylim(0, 1) +                         # proportions between 0 and 1
  labs(title = "Proportion per Celltype",
       x = "Cell Type",
       y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p)
dev.off()


#---------------
# plot B cell markers in B cells only
#---------------
b <- subset(objc, subset=celltype_f %in% bs)
ird_bobj <- readRDS('/diskmnt/Projects/myeloma_scRNA_analysis/MMY_IRD/analysis/organized/B/B_allsamples_normalized_for_DEG_filtered.rds')

bmarkers <- c('CD99', 'DNTT', 'ADA', 'RAG1', 'VPREB1', 'IGLL1', #preB dev
"TUBA1B", "TUBA1B", "STMN1", "HMGB2", "MKI67", # cycling
'CD19','SOX4', 'NEIL1', #immB dev
'ACSM3', 'HSP90AA1', 'HSP90AB1', 'JUND', 'FOS', # immb DEGs
'CD24', 'CD38', 'MME', #transitional B
'MS4A1', 'TNFRSF13C','CD52', 'BANK1', 'IGHM', 'IGHD', # nvB dev
'TCL1A', 'LTB', 'SELL', 'PLPP5', 'RGS2', 'CD83', 'IFI44L', 'IFITM1', 'PARP14', #nv
'IGHA1', 'IGHG1', 'S100A6')

p_mmrf <- DotPlot(b, features=unique(bmarkers), group.by='celltype_f', assay='RNA') + 
scale_color_gradient(low = "white", high = "black") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf("b_marker_dotplot.pdf", width=11, height=3)
p_ird <- DotPlot(ird_bobj, features=unique(bmarkers), group.by='ctm', assay='RNA') + 
scale_color_gradient(low = "white", high = "black") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf("b_marker_dotplot.pdf", width=11, height=3.5)
print(p_mmrf)
print(p_ird)
dev.off()

bi <- subset(objc, subset=IRD_subset_f %in% c('PreB', 'ImmatureB', 'NaiveB', 'MemoryB'))
p <- DotPlot(bi, features=unique(bmarkers), group.by='IRD_subset_f', assay='RNA') + 
scale_color_gradient(low = "white", high = "black") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf("b_marker_IRDlabels_dotplot.pdf", width=11, height=3)
print(p)
dev.off()

# -----------------------
# for IRD: get changes in proportion of non-PC fraction across B cells with transitional b cells separated

allm <- read.table("/diskmnt/Projects/myeloma_scRNA_analysis/MMY_IRD/analysis/organized/Overview/overview_metadata_simplified.txt", sep="\t", header=T, comment.char='', quote='')

m <- allm %>% filter(lineage != 'PC')
m[which(m$ctm %in% c('NvB_TCL1A', 'NvB_RGS2')), 'subset'] <- 'TransB'
st <- as.data.frame(table(m$sample_id, m$Timepoint))
st <- st[which(st$Freq!=0),]
names(st) <- c("Sample", "Timepoint", "nCells_Sample")
sca <- as.data.frame(table(m$sample_id, m$subset))
colnames(sca) <- c('Sample', 'CT', 'Freq')

scas <- merge(sca, st, by='Sample')
scas$Prop_of_Sample <- scas$Freq/scas$nCells_Sample*100
scas <- scas %>% filter(Timepoint %in% c('NBM', 'NDMM', 'Post-ASCT'))
scas$Timepoint <- factor(scas$Timepoint, levels=c('NBM', 'NDMM', 'Post-ASCT'))
scac$CT <- factor(scac$CT, levels= c('PreB', 'ImmatureB', 'TransB', 'NaiveB', 'MemoryB'))

my_comparisons=list(c('NBM', 'NDMM'), c('NDMM', 'Post-ASCT'), c('NBM', 'Post-ASCT'))
stat.test <- scas %>% group_by(CT) %>% 
    wilcox_test(Prop_of_Sample ~ Timepoint, comparisons = my_comparisons, p.adjust.method = "fdr") %>%
    add_significance("p.adj")

# Filter for desired CTs
selected_CTs <- c('PreB', 'ImmatureB', 'TransB', 'NaiveB', 'MemoryB')
scac <- scas %>% filter(CT %in% selected_CTs)
stat.test.filtered <- stat.test %>% filter(CT %in% selected_CTs)
stat.test.filtered <- stat.test.filtered %>% add_xy_position(x='Timepoint')
# Plot all at once with faceting
p <- ggboxplot(scac, x = "Timepoint", y = "Prop_of_Sample", outlier.shape = NA) +
  stat_pvalue_manual(stat.test.filtered, label = "p.adj.signif") +
  geom_jitter(aes(x = Timepoint, y = Prop_of_Sample),
              position = position_jitter(0.2), size = 1) +
  labs(y = "% of Sample Excl. PC") +
  theme_classic() + facet_wrap(~ CT, scales = "free_y", ncol=3) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

pdf("B_withTransB_subset_proportion_of_nonPCs_boxplot_IRD_facet.pdf", width = 6, height = 6)
print(p)
dev.off()


#----------------------

# compare pseudobulk gene expression in samples where B subtype was higher in remission vs. samples where it was lower
ab <- c('TNFRSF13C', 'TNFRSF13B', 'TNFRSF17') #april/baff
nfkb <- c('NFKB1', 'NFKB2', 'NFKBIA', 'REL', 'RELA', 'RELB', 'MYC', 'BCL2', 'CD40', 'IRF4')

# first check expression changes in NDMM vs remission for IRD B cells --------------
bref <- readRDS('/diskmnt/Projects/myeloma_scRNA_analysis/MMY_IRD/analysis/organized/B/B_mapping_reference.rds')
#pseudobulk ref for raw counts
addm <- ird_bobj@meta.data[rownames(bref@meta.data), c('ctm'), drop=F]
breff <- subset(bref, cells=rownames(addm))
breff <-AddMetaData(breff,addm)

Idents(breff) <- 'ctm'
breff <- RenameIdents(breff,
"PreB_Cycling"="PreB", "ImmB_Cycling"="ImmB", "ImmB_ACSM3"="ImmB",   "ImmB_HSP90"  = "ImmB",
"NvB_RGS2" = "TrnsB" , "NvB_TCL1A"="TrnsB",    "NvB_LTB" = "NvB",      "NvB_CD83"  = "NvB",   "NvB_BANK1"   = "NvB",
 "NvB_IFN" ="NvB",      "MemoryB" ="MemB")
breff$ctm_simple <- factor(Idents(breff), levels=c('PreB', 'ImmB', 'TrnsB', "NvB", 'MemB'))

breff$id <- paste0(breff$sample_id, "|", breff$Timepoint, "|", breff$ctm_simple)
psi <- AggregateExpression(
  object = breff,
  group.by = "id",     # or "sample_id + celltype" if needed
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)
pci <- psi$RNA  # extract pseudobulk counts matrix
pci <- pci[rownames(pci) %in% c( nfkb), ]

# Sample metadata — must match column names of pseudobulk matrix
imeta <- data.frame(id = colnames(pci))
rownames(imeta) <- imeta$id
imeta <- imeta %>% separate(id, into = c("sample", "time", 'ct'), sep = "\\|")
imeta$ct <- str_replace_all(imeta$ct, "-", "_")
imeta$sample <- str_replace_all(imeta$sample, "-", "_")

ird_avg_norm_exp <- list()
ird_res_list <- list() # compare baseline pseudobulk expression vs remission

for (ctm in levels(breff$ctm_simple)) {
  print(ctm)
  # Subset pseudobulk counts and metadata for this cell type
  ctm_pc <- pci[, rownames(imeta[which(imeta$ct == ctm & imeta$time %in% c('NDMM', 'Post-ASCT')), ])]
  ctm_meta <- imeta[which(imeta$ct==ctm & imeta$time %in% c('NDMM', 'Post-ASCT')),]
  # Create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = ctm_pc,
    colData = ctm_meta,
    design = ~ time)
  
  dds$time <- relevel(dds$time, ref = "NDMM")
  
  # Filter: keep genes with at least 3 counts in at least 2 samples
  keep_genes <- rowSums(counts(dds) >= 3) >= 2
  dds <- dds[keep_genes, ]
  
  # Filter out samples with zero total counts
  dds <- dds[, colSums(counts(dds)) > 0]
  
  # Normalize and run DESeq
  dds <- estimateSizeFactors(dds, type = "poscounts")
  normalized_counts <- counts(dds, normalized = TRUE)
  
  # Store average normalized expression
  average_expression <- setNames(data.frame(rowMeans(normalized_counts)), ctm)
  ird_avg_norm_exp[[ctm]] <- average_expression
  
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("time", "Post-ASCT", "NDMM"))
  ird_res_list[[ctm]] <- res
}

#plot
# Extract log2FoldChange vectors and combine into a single data frame
log2fc_wide <- sapply(ird_res_list, function(df) df$log2FoldChange)
rownames(log2fc_wide) <- rownames(ird_res_list[[1]])
mat <- as.matrix(log2fc_wide)
# Do the same for padj values
padj_wide <- sapply(ird_res_list, function(df) df$padj)
# Set row names (they should all be the same)
rownames(padj_wide) <- rownames(ird_res_list[[1]])

asterisk_table <- apply(padj_wide, c(1, 2), function(p) {
  if (p < 0.001) {
    "***"
  } else if (p < 0.01) {
    "**"
  } else if (p < 0.05) {
    "*"
  } else {
    ""
  }
})

hm <- Heatmap(mat, name = "log2FC", row_names_side = "left",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(as.matrix(asterisk_table)[i, j], x, y, gp = gpar(fontsize = 12))}, 
  col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE, show_row_names = TRUE)

pdf('irdB_withTrnsB_nfkbGenExpChanges_PTvsNDMM.pdf', width=5)
draw(hm)
dev.off()

# now check for MMRF B cells ----------------------------
objc$id <- paste0(objc$public_id, "|", objc$collection_event, "|", objc$celltype_f)
pseudo <- AggregateExpression(
  object = objc,
  group.by = "id",     # or "sample_id + celltype" if needed
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)

pseudo_counts <- pseudo$RNA  # extract pseudobulk counts matrix
pseudo_counts <- pseudo_counts[rownames(pseudo_counts) %in% c(ab, nfkb), ]

# Sample metadata — must match column names of pseudobulk matrix
sample_metadata <- data.frame(
  id = colnames(pseudo_counts)
)
rownames(sample_metadata) <- sample_metadata$id
sample_metadata <- sample_metadata %>%
  separate(id, into = c("public_id", "collection_event", 'celltype'), sep = "\\|")
sample_metadata$celltype <- str_replace_all(sample_metadata$celltype, "-", "_")
sample_metadata$public_id <- str_replace_all(sample_metadata$public_id, "-", "_")

avg_norm_exp <- list()
res_list_more <- list() # compare baseline pseudobulk expression vs remission
res_list_less <- list()
for(b_subtype in c('B_preLg', 'B_preSm', 'B_imm', 'B_trns', 'B_naive', 'B_Unswitched_mem', "B_mem")){
  print(b_subtype)
  b_pseudo_counts <- pseudo_counts[, rownames(sample_metadata[which(sample_metadata$celltype == b_subtype),]) ]
  dds <- DESeqDataSetFromMatrix(countData = b_pseudo_counts,
    colData = sample_metadata[which(sample_metadata$celltype == b_subtype),],
    design = ~ 1)  # no design needed for simple normalization
  dds <- dds[rowSums(counts(dds)) > 0, colSums(counts(dds)) > 0]
  dds <- estimateSizeFactors(dds,type = "poscounts")
  normalized_counts <- counts(dds, normalized=TRUE)
  average_expression <- setNames(data.frame(rowMeans(normalized_counts)), b_subtype)
  avg_norm_exp[[b_subtype]] <- average_expression

  more_in_rem <- rownames(sample_metadata[which(sample_metadata$celltype == b_subtype & sample_metadata$public_id %in% higher_rr_ids[[b_subtype]]),])
  less_in_rem <- rownames(sample_metadata[which(sample_metadata$celltype == b_subtype & ! sample_metadata$public_id %in% higher_rr_ids[[b_subtype]]),])
  b_pseudo_counts_more <-  pseudo_counts[, more_in_rem]
  b_pseudo_counts_less <-  pseudo_counts[, less_in_rem]
  b_metadata_more <- data.frame(id = colnames(b_pseudo_counts_more))
  rownames(b_metadata_more) <- b_metadata_more$id
  b_metadata_more <- b_metadata_more %>%
    separate(id, into = c("public_id", "collection_event", 'celltype'), sep = "\\|")
  b_metadata_less <- data.frame(id = colnames(b_pseudo_counts_less))
  rownames(b_metadata_less) <- b_metadata_less$id
  b_metadata_less <- b_metadata_less %>%
    separate(id, into = c("public_id", "collection_event", 'celltype'), sep = "\\|")
  dds_more <- DESeqDataSetFromMatrix(
    countData = b_pseudo_counts_more,
    colData = b_metadata_more,
    design = ~ collection_event
  )
  dds_more$collection_event <- relevel(dds_more$collection_event, ref = "Baseline")
  dds_less <- DESeqDataSetFromMatrix(
    countData = b_pseudo_counts_less,
    colData = b_metadata_less,
    design = ~ collection_event
  )
  dds_less$collection_event <- relevel(dds_less$collection_event, ref = "Baseline")
  keep_genes_more <- rowSums(counts(dds_more) >= 3) >= 2  # At least 5 counts in 2+ samples
  dds_more <- dds_more[keep_genes_more, ]
  sample_zeros_more <- colSums(counts(dds_more)) == 0
  dds_more <- dds_more[, !sample_zeros_more]
  dds_more <- estimateSizeFactors(dds_more, type = "poscounts")
  dds_more <- DESeq(dds_more)
  res_more <- results(dds_more, contrast=c("collection_event", "Remission/Response", "Baseline"))
  res_more_df <- as.data.frame(res_more[order(res_more$pvalue), ])
  res_more_df$nonPCpropComparedToBaseline <- 'higher'
  res_list_more[[b_subtype]] <- res_more_df
  keep_genes_less <- rowSums(counts(dds_less) >= 3) >= 2  # At least 5 counts in 2+ samples
  dds_less <- dds_less[keep_genes_less, ]
  sample_zeros_less <- colSums(counts(dds_less)) == 0
  dds_less <- dds_less[, !sample_zeros_less]
  dds_less <- estimateSizeFactors(dds_less, type = "poscounts")
  dds_less <- DESeq(dds_less)
  res_less <- results(dds_less, contrast=c("collection_event", "Remission/Response", "Baseline"))
  res_less_df <- as.data.frame(res_less[order(res_less$pvalue), ])
  res_less_df$nonPCpropComparedToBaseline <- 'lower'
  res_list_less[[b_subtype]] <- res_less_df
}

# plot changes -------------------------
# get maximum log2 average normalized expression of all genes of interest
max_values <- sapply(avg_norm_exp, function(df) max(unlist(df), na.rm = TRUE))
overall_max <- max(max_values) 
log2(overall_max) # 7.112

pdf('Bsubtype_geneExp_changes_splitBy_proportionHigherOrLower_relativeToNDMM.pdf', width=5)
for(b_subtype in c('B_preLg', 'B_preSm', 'B_imm', 'B_trns', 'B_naive', 'B_Unswitched_mem', "B_mem")){
more <- res_list_more[[b_subtype]]
more$gene <- rownames(more)
less <- res_list_less[[b_subtype]]
less$gene <- rownames(less)
all <- rbind(more, less)
alldf <- all %>% mutate(sig_label = ifelse(padj < 0.05, "*", ""))

# Only keep relevant columns
df_clean <- alldf %>% select(gene, nonPCpropComparedToBaseline, log2FoldChange, padj)
rownames(df_clean) <- NULL
# Define conditions
conditions <- c("higher", "lower")
# Create a complete grid of all genes × conditions
full_grid <- expand.grid(
  gene = unique(c(ab, nfkb)),
  nonPCpropComparedToBaseline = conditions,
  stringsAsFactors = FALSE)
# Join existing data to the full grid
df_complete <- full_grid %>%
  left_join(df_clean, by = c("gene", "nonPCpropComparedToBaseline")) %>%
  arrange(gene, nonPCpropComparedToBaseline)

# Make a matrix for heatmap
heatmap_matrix <- df_complete %>%
  select(gene, nonPCpropComparedToBaseline, log2FoldChange) %>%
  pivot_wider(names_from = nonPCpropComparedToBaseline, values_from = log2FoldChange) %>%
  column_to_rownames("gene") %>% as.matrix()

mat1 <- heatmap_matrix[intersect(nfkb, rownames(heatmap_matrix)),]
mat2 <- heatmap_matrix[intersect(ab, rownames(heatmap_matrix)),]

# Create a matrix of labels ("*" where padj < 0.05, "" otherwise)
sig_label_matrix <- df_complete %>%
  mutate(
    sig_label = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ "")) %>%
  select(gene, nonPCpropComparedToBaseline, sig_label) %>%
  pivot_wider(names_from = nonPCpropComparedToBaseline, values_from = sig_label) %>%
  column_to_rownames("gene") %>%  as.matrix()
sig1 <- sig_label_matrix[nfkb,]
sig2 <- sig_label_matrix[ab,]

# get log2 mean normalized expression across all samples where the gene is sufficiently detected:
avgexp <- avg_norm_exp[[b_subtype]]
avgexp1 <- log2(avgexp[rownames(mat1),] +1)
names(avgexp1) <- rownames(mat1)
row_anno1 <- rowAnnotation(
  "log2(Mean Norm. Exp)" = anno_barplot(avgexp1, 
                                  border = FALSE, ylim=c(0,7.2),
                                  gp = gpar(fill = "skyblue")))
avgexp2 <- log2(avgexp[rownames(mat2),] +1)
names(avgexp2) <- rownames(mat2)
row_anno2 <- rowAnnotation(
  "log2(Mean Norm. Exp)" = anno_barplot(avgexp2, 
                                  border = FALSE, ylim=c(0,7.2),
                                  gp = gpar(fill = "skyblue")))
hm1 <- Heatmap( mat1, name = "log2FC", right_annotation = row_anno1,row_names_side = "left",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sig1[i, j], x, y, gp = gpar(fontsize = 12))}, 
  col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = TRUE, show_row_names = TRUE)
hm2 <- Heatmap(mat2, name = "log2FC",  right_annotation = row_anno2, row_names_side = "left",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sig2[i, j], x, y, gp = gpar(fontsize = 12))},
  col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), cluster_rows = FALSE,
  cluster_columns = FALSE, show_column_names = TRUE, show_row_names = TRUE)
draw(hm1 %v% hm2,
   column_title=b_subtype,
   column_title_gp=grid::gpar(fontsize=16))
}
dev.off()

#---------------------------
# look at other cell types
#---------------------------
# look at marker genes
tnk_obj <- subset(objc, subset= celltype_f %in% tnks)

tm <- unique(c('CD3D', 'CD3E', 'CD4', 'SELL','CCR7', 'TCF7', 'FOXP1', # cd4 tnv0
'ANXA1', 'S100A11', 'TIMP1', 'IL7R', 'CD40LG', #cd4 tmem
'FOXP3', 'CTLA4','IL2RA', # treg
'CD8A', 'CD8B', 'GZMA', 'GZMB', 'GZMH', 'GZMK', 'TIGIT', 'LAG3', 'EOMES',
'TGFB1', 'LEPROTL1','ZEB1','ZEB2', # CD4, CD8 TGFB
'ISG15', 'IFI44L', 'IFI6', 'IFITM1', 'STAT1', #CD4 IFN
'IFIT2', 'IFIT3', 'CCL5', 'TNF', 'IFNG', # CD8 IFN
'KLRB1', 'DPP4', 'SLC4A10', #MAIT
'TRGC1', 'TRDC', #gdT
'NCAM1', 'FCGR3A', 'PRF1', 'GZMA', 'GZMB', 'GNLY', 'SPON2', 'CST7',  # CD56dim
'KLRC1', 'SELL', 'CD69', # CD56bright?
'STMN1', 'MKI67', 'TOP2A'))

p <- DotPlot(tnk_obj, features=unique(tm), group.by='celltype_f', assay='RNA') + 
scale_color_gradient(low = "white", high = "black") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf("tnk_marker_dotplot.pdf", width=11, height=10)
print(p)
dev.off()

# proportion changes over time
tnk_meta <- tnk_obj@meta.data
tnk_meta$sample_id <- paste0(tnk_meta$public_id,"|",tnk_meta$collection_event)
st <- as.data.frame(table(tnk_meta$sample_id))
names(st) <- c("Sample", "nCells_Sample")
sca <- as.data.frame(table(tnk_meta$sample_id, tnk_meta$celltype_f))
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
selected_CTs <- tnks
stat.test.filtered <- stat.test %>% filter(CT %in% selected_CTs & p.adj.signif !='ns')
stat.test.filtered <- stat.test.filtered %>% add_xy_position(x='collection_event')
scac <- scas %>% filter(CT %in% unique(as.character(stat.test.filtered$CT)) & collection_event %in% c('Baseline', 'Remission/Response') )

# Plot all at once with faceting
p <- ggpaired(scac, x='collection_event', y='Prop_of_Sample', id = "public_id", line.color = "gray", line.size = 0.5, dot.size = 3,  palette = "jco") + 
  stat_pvalue_manual(stat.test.filtered, label = "p.adj.signif") + facet_wrap(~ CT, scales ="free_y", ncol=4) +
  labs(x = "Collection Event", y = "% of total T/NK") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("mmrf_TNK_proportions_overTime.pdf", width = 8, height = 12)
print(p)
dev.off()


# simplify t cell annotations

objc$ct_simple <- as.character(objc$celltype_f)
objc@meta.data[which(objc@meta.data$ct_simple %in% cd4s), 'ct_simple'] <- 'CD4T'
objc@meta.data[which(objc@meta.data$ct_simple %in% cd8s), 'ct_simple'] <- 'CD8T'

# get CD4:CD8 ratios across timepoints
objc_meta <- objc@meta.data
objc_meta$sample_id <- paste0(objc_meta$public_id,"|",objc_meta$collection_event)

sca <- as.data.frame(table(objc_meta$sample_id, objc_meta$ct_simple))
colnames(sca) <- c('Sample', 'CT', 'Freq')

cd4_cd8_ratio <- sca %>%
  filter(CT %in% c("CD4T", "CD8T")) %>%
  select(Sample, CT, Freq) %>%
  pivot_wider(names_from = CT, values_from = Freq) %>%
  mutate(CD4_CD8_ratio = CD4T / CD8T)

cd4_cd8_ratio <- cd4_cd8_ratio %>%
  separate(Sample, into = c("public_id", "collection_event"), sep = "\\|")

cd4_cd8_ratio$CD4_CD8_ratio <- log2(cd4_cd8_ratio$CD4_CD8_ratio)
cd4_cd8_ratio$collection_event <- factor(cd4_cd8_ratio$collection_event, levels=c("Baseline", "Post-transplant", "Remission/Response", "Relapse/Progression" ))
# Plot all at once with faceting
p <- ggpaired(cd4_cd8_ratio, x='collection_event', y='CD4_CD8_ratio', id = "public_id", line.color = "gray", line.size = 0.5, dot.size = 3,  palette = "jco") + stat_compare_means(method='wilcox.test', label='p.signif', comparisons=list(c('Baseline', 'Remission/Response'))) +
  labs(x = "Collection Event", y = "CD4:CD8 Ratio") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("mmrf_cd4cd8ratio_overTime.pdf", width=3, height=4)
print(p)
dev.off()

# get samples with higher vs lower ratio at PT
# get samples where the ratio increases or decreases
tratio <- cd4_cd8_ratio[which(cd4_cd8_ratio$collection_event %in% c('Baseline', 'Remission/Response')), c('public_id', 'collection_event', 'CD4_CD8_ratio')]

wide_df <- tratio %>%
  pivot_wider(
    names_from = collection_event,
    values_from = CD4_CD8_ratio
  ) %>% na.omit()

wide_df$higher <- wide_df$Baseline < wide_df$`Remission/Response`
wide_df <- wide_df %>% arrange(wide_df$higher)
higher_samples <- wide_df %>% filter(higher==TRUE) %>% select(public_id)
lower_samples <- wide_df %>% filter(higher==FALSE) %>% select(public_id)
objc@meta.data$T48ratio_PTvsNDMM <- 'Unk'
objc@meta.data[which(objc@meta.data$public_id %in% higher_samples$public_id), 'T48ratio_PTvsNDMM'] <- 'Higher'
objc@meta.data[which(objc@meta.data$public_id %in% lower_samples$public_id), 'T48ratio_PTvsNDMM'] <- 'Lower'

# check T cell gene expression in higher vs lower groups
core_markers <- c('CD2', 'CD3D', 'CD3E', 'CD4', 'CD8A')
costim <- c('CD27', 'CD28', 'SLAMF1', 'SLAMF7', 'TRAC', 'PTPRC', 'IL7R', 'FOXP3')
cytox <- c('PRF1', 'GZMA', 'GZMB', 'GZMK', 'NKG7')
immreg <- c('TGFB1', 'PDCD1', 'CTLA4', 'LAG3', 'HAVCR2', 'EOMES', 'TIGIT', 'FAS', 'BTNL9')
chemok <- c('CCR7', 'CXCR4', 'CCL5', 'CYTIP', 'KLRB1')
activ <- c('CD69', 'CD247', 'DUSP2', 'IL2RA', 'RGS16', 'PRDM1', 'SPI1', 'SPIB')
ifn <- c('IFI6', 'IFI44L', 'ISG15')

objc$ratiogroup_ct <- paste0(objc$T48ratio_PTvsNDMM, "|", objc$ct_simple, "|", objc$collection_event)
Idents(objc) <- 'ratiogroup_ct'
# A positive avg_log2FC means the gene is upregulated in higher|CD8+ T-Cell compared to lower|CD8+ T-Cell.
cd8_de <- FindMarkers(
  object = objc,
  ident.1 = "Higher|CD8T|Remission/Response",   # Replace with your actual group names
  ident.2 = "Lower|CD8T|Remission/Response",
  features = unique(c(immreg, cytox, activ)),
  logfc.threshold = 0,  # set to 0 to get all genes in the list
  min.pct = 0.05,        # only test genes expressed in >10% of cells
  test.use = "MAST", # or "MAST", "DESeq2", "LR", etc.
)
cd8_de <- cd8_de %>% arrange(avg_log2FC) %>% filter(p_val_adj < 0.05)
cd8_de$pct_diff <- (cd8_de$pct.1 - cd8_de$pct.2) * 100
cd8_de$gene <- rownames(cd8_de)

cd8_de_n <- FindMarkers(
  object = objc,
  ident.1 = "Higher|CD8T|Baseline",   # Replace with your actual group names
  ident.2 = "Lower|CD8T|Baseline",
  features = unique(c(immreg, cytox, activ)),
  logfc.threshold = 0,  # set to 0 to get all genes in the list
  min.pct = 0.05,        # only test genes expressed in >10% of cells
  test.use = "MAST", # or "MAST", "DESeq2", "LR", etc.
)
cd8_de_n <- cd8_de_n %>% arrange(avg_log2FC) %>% filter(p_val_adj < 0.05)
cd8_de_n$pct_diff <- (cd8_de_n$pct.1 - cd8_de_n$pct.2) * 100
cd8_de_n$gene <- rownames(cd8_de_n)

cd4_de <- FindMarkers(
  object = objc,
  ident.1 = "Higher|CD4T|Remission/Response",   # Replace with your actual group names
  ident.2 = "Lower|CD4T|Remission/Response",
  features = unique(c(immreg, cytox, activ)),
  logfc.threshold = 0,  # set to 0 to get all genes in the list
  min.pct = 0.05,        # only test genes expressed in >10% of cells
  test.use = "MAST", # or "MAST", "DESeq2", "LR", etc.
)
cd4_de <- cd4_de %>% arrange(avg_log2FC) %>% filter(p_val_adj < 0.05)
cd4_de$pct_diff <- (cd4_de$pct.1 - cd4_de$pct.2) * 100
cd4_de$gene <- rownames(cd4_de)

cd4_de_n <- FindMarkers(
  object = objc,
  ident.1 = "Higher|CD4T|Baseline",   # Replace with your actual group names
  ident.2 = "Lower|CD4T|Baseline",
  features = unique(c(immreg, cytox, activ)),
  logfc.threshold = 0,  # set to 0 to get all genes in the list
  min.pct = 0.05,        # only test genes expressed in >10% of cells
  test.use = "MAST", # or "MAST", "DESeq2", "LR", etc.
)
cd4_de_n <- cd4_de_n %>% arrange(avg_log2FC) %>% filter(p_val_adj < 0.05)
cd4_de_n$pct_diff <- (cd4_de_n$pct.1 - cd4_de_n$pct.2) * 100
cd4_de_n$gene <- rownames(cd4_de_n)


# plot log2FC and change in % expressed 

cd8_de$gene_f <- factor(cd8_de$gene, levels=cd8_de$gene)
cd8_de_n$gene_f <- factor(cd8_de_n$gene, levels=cd8_de_n$gene)
p1 <- ggplot(cd8_de, aes(x = gene_f, y = pct_diff, fill = avg_log2FC)) +
  geom_bar(stat = "identity") +                # Use 'identity' to plot actual values
  scale_fill_gradient2(                        # Blue-white-red gradient centered at zero
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,limits=c(-2.2, 1.3),
    name = "avg_log2FC"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)  # Rotate x labels if many
  ) +
  labs(
    y = "Difference in % Expression",
    x = "Gene", title="CD8+ T Cells: Remission"
  )
p2 <- ggplot(cd8_de_n, aes(x = gene_f, y = pct_diff, fill = avg_log2FC)) +
  geom_bar(stat = "identity") +                # Use 'identity' to plot actual values
  scale_fill_gradient2(                        # Blue-white-red gradient centered at zero
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,limits=c(-2.2, 1.3),
    name = "avg_log2FC"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)  # Rotate x labels if many
  ) +
  labs(
    y = "Difference in % Expression",
    x = "Gene", title="CD8+ T Cells: Baseline"
  ) 

pdf('paired_cd8T_DEGs_higherVslowerTratio_remissionVsbaseline.pdf', width=8, height=3)
print(plot_grid(p2, p1))
dev.off()




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


# compare pseudobulk gene expression of immune suppressive / exhaustion markers in T cells
m <- c('TGFB1', 'EOMES', 'TIGIT', 'LAG3', 'HAVCR2', 'TOX', 'PRDM1', 'GZMK')

objc$id <- paste0(objc$public_id, "|", objc$collection_event, "|", objc$ct_simple)
pseudo <- AggregateExpression(
  object = objc,
  group.by = "id",     # or "sample_id + celltype" if needed
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE)
pseudo_counts <- pseudo$RNA  # extract pseudobulk counts matrix
pseudo_counts <- pseudo_counts[rownames(pseudo_counts) %in% m, ]

# Sample metadata — must match column names of pseudobulk matrix
sample_metadata <- data.frame(
  id = colnames(pseudo_counts)
)
rownames(sample_metadata) <- sample_metadata$id
sample_metadata <- sample_metadata %>%
  separate(id, into = c("public_id", "collection_event", 'celltype'), sep = "\\|")
sample_metadata$celltype <- str_replace_all(sample_metadata$celltype, "-", "_")
sample_metadata$public_id <- str_replace_all(sample_metadata$public_id, "-", "_")

avg_norm_exp <- list()
res_list <- list() # compare baseline pseudobulk expression vs remission

for(t_subtype in c('CD4T', 'CD8T')){
  print(t_subtype)
  t_pseudo_counts <- pseudo_counts[, rownames(sample_metadata[which(sample_metadata$celltype == t_subtype),]) ]
  dds <- DESeqDataSetFromMatrix(countData = t_pseudo_counts,
    colData = sample_metadata[which(sample_metadata$celltype == t_subtype),],
    design = ~ 1)  # no design needed for simple normalization
  dds <- dds[rowSums(counts(dds)) > 0, colSums(counts(dds)) > 0]
  dds <- estimateSizeFactors(dds,type = "poscounts")
  normalized_counts <- counts(dds, normalized=TRUE)
  average_expression <- setNames(data.frame(rowMeans(normalized_counts)), t_subtype)
  avg_norm_exp[[t_subtype]] <- average_expression

  t_metadata <- data.frame(id = colnames(t_pseudo_counts))
  rownames(t_metadata) <- t_metadata$id
  t_metadata <- t_metadata %>%
    separate(id, into = c("public_id", "collection_event", 'celltype'), sep = "\\|")
  dds <- DESeqDataSetFromMatrix(
    countData = t_pseudo_counts,
    colData = t_metadata,
    design = ~ collection_event)

  dds$collection_event <- relevel(dds$collection_event, ref = "Baseline")
  keep_genes <- rowSums(counts(dds) >= 3) >= 2  # At least 5 counts in 2+ samples
  dds <- dds[keep_genes, ]
  sample_zeros <- colSums(counts(dds)) == 0
  dds <- dds[, !sample_zeros]
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("collection_event", "Remission/Response", "Baseline"))
  res_df <- as.data.frame(res[order(res$pvalue), ])
  res_list[[t_subtype]] <- res_df
}


# map MMRF T cells to IRD T cell UMAP
tref <- readRDS('/diskmnt/Projects/myeloma_scRNA_analysis/MMY_IRD/analysis/organized/TNK/TNK_mapping_reference_reannotated.rds')

set.seed(666)
query.obj <- tnk_obj
ref.obj <- tref

DefaultAssay(ref.obj) <- 'integrated'
DefaultAssay(query.obj) <- 'SCT'
anchors <- FindTransferAnchors(
  reference = ref.obj,
  query = query.obj,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:50
)
query.obj <- MapQuery(
  anchorset = anchors,
  query = query.obj,
  reference = ref.obj,
  refdata = list(
    celltype.l1 = 'ctm'
  ),
  reference.reduction = "pca",
  reduction.model = "umap"
)

colnames(query.obj@meta.data)[colnames(query.obj@meta.data) == "predicted.celltype.l1"] = "IRD_ctm"
colnames(query.obj@meta.data)[colnames(query.obj@meta.data) == "predicted.celltype.l1.score"] = "IRD_ctm_score"
tnk_obj <- query.obj
saveRDS(tnk_obj,"mmrf_tnk.rds", compress=F)

# plot T/NK markers again with mapped IRD ctm annotation
tnk_obj$IRD_ctm <- factor(tnk_obj$IRD_ctm, levels=levels(tref$ctm))
p <- DotPlot(tnk_obj, features=unique(tm), group.by='IRD_ctm', assay='RNA') + 
scale_color_gradient(low = "white", high = "black") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf("tnk_marker_dotplot_IRDmappedCtm.pdf", width=12, height=5)
print(p)
dev.off()

# Count number of cells per cell type per collection_event per public ID for IRD ctm annotation
tnk_meta <- tnk_obj@meta.data
tnk_meta$sample_id <- paste0(tnk_meta$public_id,"|",tnk_meta$collection_event)
st <- as.data.frame(table(tnk_meta$sample_id))
names(st) <- c("Sample", "nCells_Sample")
sca <- as.data.frame(table(tnk_meta$sample_id, tnk_meta$IRD_ctm))
colnames(sca) <- c('Sample', 'CT', 'Freq')

scas <- merge(sca, st, by='Sample')
scas$Prop_of_Sample <- scas$Freq/scas$nCells_Sample*100
scas <- scas %>%
  separate(Sample, into = c("public_id", "collection_event"), sep = "\\|")

scas$collection_event <- factor(scas$collection_event, levels=c('Baseline', 'Post-transplant', 'Remission/Response', 'Relapse/Progression' ))
scac <- scas %>% filter()
# Perform statistical test and adjust p-values
my_comparisons <- list(c('Baseline', 'Remission/Response'))
stat.test <- scas %>% group_by(CT) %>% 
    wilcox_test(Prop_of_Sample ~ collection_event, comparisons = my_comparisons, p.adjust.method = "fdr") %>%
    add_significance("p.adj")

# Filter for desired CTs
selected_CTs <- c('CD4_TGFB1', 'CD8_TGFB1')
stat.test.filtered <- stat.test %>% filter(CT %in% selected_CTs & p.adj.signif !='ns')
stat.test.filtered <- stat.test.filtered %>% add_xy_position(x='collection_event')
scac <- scas %>% filter(CT %in% selected_CTs & collection_event %in% c('Baseline', 'Remission/Response'))

# Plot all at once with faceting
p <- ggpaired(scac, x='collection_event', y='Prop_of_Sample', id = "public_id", line.color = "gray", line.size = 0.5, dot.size = 3,  palette = "jco") + stat_compare_means(method='wilcox.test', label='p.signif') + facet_wrap(~ CT, scales ="free_y", ncol=4) +
  labs(x = "Collection Event", y = "% of total T/NK") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("mmrf_tnk_IRD_tgfb1_proportions_overTime.pdf", width=3, height=4)
print(p)
dev.off()


pdf("mmrf_tnk_IRD_tgfb1_proportions_overTime.pdf", width=3, height=4)
for(ct in c('CD4_TGFB1', 'CD8_TGFB1')){
  ctf <- scas_f[which(scas_f$CT==ct), c('public_id', 'collection_event', 'Prop_of_Sample')]
  ctf <- ctf %>%
    group_by(public_id) %>%
    filter(all(c("Baseline", "Remission/Response") %in% collection_event)) %>%
    ungroup()
  ctf$time <- as.character(ctf$collection_event)
  ctf[which(ctf$collection_event=='Remission/Response'),'time'] <- 'Remission'
  p <- ggpaired(ctf, 
                x='time', 
                y='Prop_of_Sample',
                id = "public_id",
                line.color = "gray", 
                line.size = 0.5, dot.size = 4, 
                palette = "jco") + 
    stat_compare_means(method = "wilcox.test", p.adjust.method = "fdr", paired = TRUE, label = "p.signif") +
    labs(title = ct, x = "Collection Event", y = "% of total T/NK") #+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(p)
}
dev.off()

#------------------
# Look at BAFF/APRIL expression in myeloid cells
mye <- ctlevels$ct[59:72]

ablig <- c('TNFSF13', 'TNFSF13B')
objc@meta.data[which(objc@meta.data$ct_simple %in% mye), 'ct_simple'] <- 'Myeloid'


objc$id <- paste0(objc$public_id, "|", objc$collection_event, "|", objc$ct_simple)
pseudo <- AggregateExpression(
  object = objc,
  group.by = "id",     # or "sample_id + celltype" if needed
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)

pseudo_counts <- pseudo$RNA  # extract pseudobulk counts matrix
pseudo_counts <- pseudo_counts[rownames(pseudo_counts) %in% c(ablig), ]


# Sample metadata — must match column names of pseudobulk matrix
pmeta <- data.frame(id = colnames(pseudo_counts))
rownames(pmeta) <- pmeta$id
pmeta <- pmeta %>% separate(id, into = c("public_id", "collection_event", 'ct'), sep = "\\|")
pmeta$ct <- str_replace_all(pmeta$ct, "-", "_")
pmeta$public_id <- str_replace_all(pmeta$public_id, "-", "_")

avg_norm_exp <- list()
res_list <- list() # compare baseline pseudobulk expression vs remission

# Subset pseudobulk counts and metadata for myeloid
ctm_pc <- pseudo_counts[, rownames(pmeta[which(pmeta$ct == "Myeloid" & pmeta$collection_event %in% c('Baseline', 'Remission/Response')), ])]

expr_long <- as.data.frame(ctm_pc) %>%
  rownames_to_column("Gene") %>%             # Make gene names a column
  pivot_longer(
    cols = -Gene,                            # All other columns are expression values
    names_to = "SampleInfo", 
    values_to = "Expression"
  ) %>%
  separate(SampleInfo, into = c("SampleID", "Condition", "CellType"), sep = "\\|")

# Plot all at once with faceting
p <- ggpaired(expr_long, x='Condition', y='Expression', id = "SampleID", line.color = "gray", line.size = 0.5, dot.size = 3,  palette = "jco") + 
  stat_compare_means(method='wilcox.test', label='p.signif') + facet_wrap(~ Gene, scales ="free_y", ncol=4) +
  labs(x = "Collection Event", y = "Pseudobulk Expression") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("mmrf_myeloid_aprilBaff_overTime.pdf", width = 8, height = 12)
print(p)
dev.off()

ctm_meta <- pmeta[colnames(ctm_pc),]
# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = ctm_pc,
  colData = ctm_meta,
  design = ~ collection_event)

dds$collection_event <- relevel(dds$collection_event, ref = "Baseline")

# Filter: keep genes with at least 3 counts in at least 2 samples
keep_genes <- rowSums(counts(dds) >= 3) >= 2
dds <- dds[keep_genes, ]

# Filter out samples with zero total counts
dds <- dds[, colSums(counts(dds)) > 0]

# Normalize and run DESeq
dds <- estimateSizeFactors(dds, type = "poscounts")
normalized_counts <- counts(dds, normalized = TRUE)

# Store average normalized expression
average_expression <- setNames(data.frame(rowMeans(normalized_counts)), 'Myeloid')

dds <- DESeq(dds, fitType = "mean")
res <- results(dds, contrast = c("collection_event", "Remission/Response", "Baseline"))



#plot
# Extract log2FoldChange vectors and combine into a single data frame
log2fc_wide <- sapply(ird_res_list, function(df) df$log2FoldChange)
rownames(log2fc_wide) <- rownames(ird_res_list[[1]])
mat <- as.matrix(log2fc_wide)
# Do the same for padj values
padj_wide <- sapply(ird_res_list, function(df) df$padj)
# Set row names (they should all be the same)
rownames(padj_wide) <- rownames(ird_res_list[[1]])

asterisk_table <- apply(padj_wide, c(1, 2), function(p) {
  if (p < 0.001) {
    "***"
  } else if (p < 0.01) {
    "**"
  } else if (p < 0.05) {
    "*"
  } else {
    ""
  }
})

hm <- Heatmap(mat, name = "log2FC", row_names_side = "left",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(as.matrix(asterisk_table)[i, j], x, y, gp = gpar(fontsize = 12))}, 
  col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE, show_row_names = TRUE)

pdf('irdB_withTrnsB_nfkbGenExpChanges_PTvsNDMM.pdf', width=5)
draw(hm)
dev.off()