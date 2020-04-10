
diff_expr_analysis <- function(
  read_counts_df,
  metadata_df,
  gene_df,
  
  group_str,
  numerator_level_str,
  denominator_level_str,
  coef,

  design_formula,
  pv_adj_threshold,
  dir_output
) {
  metadata_df$Line <- factor(metadata_df$Line)
  metadata_df$Condition <- factor(metadata_df$Condition)
  
  # edgeR
  edgeR.DGElist <- DGEList(
    counts = read_counts_df,
    samples = metadata_df,
    group = metadata_df[, group_str]
  )
  
  design <- model.matrix(
    data = edgeR.DGElist$samples,
    design_formula
  )
  
  # Filtering
  to_keep_list  <- rowSums(read_counts_df >= 5) >= (min(table(metadata_df[, group_str])) / 2)
  edgeR.DGElist <- edgeR.DGElist[to_keep_list, , keep.lib.sizes=FALSE]

  edgeR.DGElist <- calcNormFactors(edgeR.DGElist, method = "TMM")
  edgeR.DGElist <- estimateDisp(edgeR.DGElist, design)
  
  edger_ql_fit <- glmQLFit(edgeR.DGElist, design)
  edger_qlf <- glmQLFTest(edger_ql_fit, coef = coef)
  
  edgeR.ql.res.sorted <- topTags(
    edger_qlf,
    n = Inf , # To retrieve all genes
    p.value = 1,
    sort.by = "PValue", adjust.method = "BH"
  )
  
  # DESeq2
  DESeq.ds <- DESeqDataSetFromMatrix(
    countData = read_counts_df,
    colData = metadata_df,
    design = design_formula
  )
  
  DESeq.ds <- DESeq.ds[to_keep_list, ]
  DESeq.ds <- DESeq(DESeq.ds)
  
  coef_deseq2 <- paste(group_str, numerator_level_str, 'vs', denominator_level_str, sep = '_')
  
  DESeq.res <- results(
    object = DESeq.ds,
    contrast = c(group_str, numerator_level_str, denominator_level_str),
    independentFiltering = TRUE,
    alpha = pv_adj_threshold,
    pAdjustMethod = "BH",
  )
  DESeq.res_lfcShrink <- lfcShrink(
    res = DESeq.res,
    dds = DESeq.ds,
    coef = coef_deseq2,
    type = "apeglm",
  )
  
  # limma-voom
  voomTransformed <- voom(edgeR.DGElist, design, plot = FALSE)
  voomed.fitted <- lmFit(voomTransformed, design = design)
  
  voomed.fitted <- eBayes(voomed.fitted)
  
  limmaV.res.sorted <- topTable(
    voomed.fitted, coef = coef,
    number = Inf,
    p.value = 1,
    sort.by = "P", adjust.method = "BH"
  )
  
  # Results
  dir_xxx = file.path(
    dir_output,
    paste0('des_', gsub("\\ \\+\\ ", "", as.character(design_formula)[2]))
  )

  DE_list <- list(
    LimmaV = rownames(subset(limmaV.res.sorted, adj.P.Val <= pv_adj_threshold)),
    edgeR = rownames(subset(edgeR.ql.res.sorted$table, FDR <= pv_adj_threshold)),
    DESeq2 = rownames(subset(DESeq.res_lfcShrink, padj <= pv_adj_threshold))
  )
  DEG_All <- Reduce(union, DE_list)
  
  DEG_All.expFCpadj <- data.frame(
    LimmaV_AveExpr = limmaV.res.sorted[DEG_All, ]$AveExpr,
    edgeR_logCPM = edgeR.ql.res.sorted[DEG_All, ]$table$logCPM,
    DESeq2_baseMean = DESeq.res_lfcShrink[DEG_All, ]$baseMean,
    
    LimmaV_FC = limmaV.res.sorted[DEG_All, ]$logFC,
    edgeR_FC = edgeR.ql.res.sorted[DEG_All, ]$table$logFC,
    DESeq2_FC = DESeq.res_lfcShrink[DEG_All, ]$log2FoldChange,
    
    LimmaV_padj = limmaV.res.sorted[DEG_All, ]$adj.P.Val,
    edgeR_padj = edgeR.ql.res.sorted[DEG_All, ]$table$FDR,
    DESeq2_padj = DESeq.res_lfcShrink[DEG_All, ]$padj,
    
    row.names = DEG_All
  )
  
  den_sample_list <- metadata_df %>% rownames_to_column(var = 'sample') %>% filter(!!sym(group_str) == denominator_level_str) %>% pull(sample)
  num_sample_list <- metadata_df %>% rownames_to_column(var = 'sample') %>% filter(!!sym(group_str) == numerator_level_str) %>% pull(sample)
  
  DEG_All.expFCpadj <- merge(
    DEG_All.expFCpadj %>% rownames_to_column(var = 'gene_id'),
    read_counts_df %>% rownames_to_column(var = 'gene_id') %>% select(all_of(c('gene_id', den_sample_list, num_sample_list))),
    by = 'gene_id'
  ) %>% column_to_rownames(var = 'gene_id')
  

  dir.create(dir_xxx, recursive = TRUE, showWarnings = FALSE)
  
  png(file.path(dir_xxx, paste0(denominator_level_str, 'Vs', numerator_level_str, '_edgeR-DESeq2-limmaVoom_venn.png')))
  gplots::venn(DE_list)
  dev.off()
  
  dir.create(dir_xxx, recursive = TRUE, showWarnings = FALSE)
  
  # Get the values that appear in two or more lists
  rl <- rle(sort(unlist(DE_list)))
  DEG_atLeast_i_samples <- rl$values[rl$lengths >= 2]
  
  x_df <- DEG_All.expFCpadj[DEG_atLeast_i_samples,]
  DEG_atLeast_i_samples.under <- rownames(x_df[x_df$DESeq2_FC < 0, ])
  DEG_atLeast_i_samples.over <- rownames(x_df[x_df$DESeq2_FC > 0, ])
  
  print(paste('All:', length(DEG_atLeast_i_samples), '=', length(DEG_atLeast_i_samples.under), '(under) +', length(DEG_atLeast_i_samples.over), '(over)'))
  print(dir_xxx)
  

  merge(
    ensembl_annotation_df[DEG_atLeast_i_samples, ] %>% rownames_to_column(var = 'gene_id'),
    DEG_All.expFCpadj[DEG_atLeast_i_samples, ] %>% rownames_to_column(var = 'gene_id'),
    by = 'gene_id'
  ) %>% arrange(desc(DESeq2_FC)) %>%
    write_tsv(
      file.path(
        dir_xxx,
        paste0(denominator_level_str, 'Vs', numerator_level_str, '_genes.tsv')
      ),
      quote = FALSE
    )

  merge(
    ensembl_annotation_df[DEG_atLeast_i_samples.under, ] %>% rownames_to_column(var = 'gene_id'),
    DEG_All.expFCpadj[DEG_atLeast_i_samples.under, ] %>% rownames_to_column(var = 'gene_id'),
    by = 'gene_id'
  ) %>% arrange(desc(DESeq2_FC)) %>%
    write_tsv(
      file.path(
        dir_xxx,
        paste0(denominator_level_str, 'Vs', numerator_level_str, '_genes.under.tsv')
      ),
      quote = FALSE
    )
  
  merge(
    ensembl_annotation_df[DEG_atLeast_i_samples.over, ] %>% rownames_to_column(var = 'gene_id'),
    DEG_All.expFCpadj[DEG_atLeast_i_samples.over, ] %>% rownames_to_column(var = 'gene_id'),
    by = 'gene_id'
  ) %>% arrange(desc(DESeq2_FC)) %>%
    write_tsv(
      file.path(
        dir_xxx,
        paste0(denominator_level_str, 'Vs', numerator_level_str, '_genes.over.tsv')
      ),
      quote = FALSE
    )
}


library(DESeq2)
library(edgeR)
library(limma)

library(tidyverse)


path_ensembl_annotation = '../data/ensembl_annotation.tsv'
path_metadata = '../data/metadata.txt'

# The matrix is on the VM for now. When it will be shared through other means, it will be specified the directory where take it
path_read_count = 'gene_count_matrix.csv'

ensembl_annotation_df <- read_csv(path_ensembl_annotation) %>% column_to_rownames(var = 'gene_id')

read_counts_df <- read_csv(file = path_read_count) %>%
  column_to_rownames(var = 'gene_id')


metadata_df <- read.table(path_metadata, header = T, row.names = 1, sep = ',')
metadata_df$Condition <- factor(metadata_df$Condition, levels = c("Control", "Infected"))





group_str = 'Condition'
numerator_level_str = 'Infected'
denominator_level_str = 'Control'
coef = 'ConditionInfected'
pv_adj_threshold = 0.05


# All
diff_expr_analysis(
  read_counts_df = read_counts_df,
  metadata_df = metadata_df,

  gene_df = ensembl_annotation_df,
  group_str = group_str,
  numerator_level_str = numerator_level_str,
  denominator_level_str = denominator_level_str,
  coef = coef,
  
  design_formula = ~ Line + Batch + Condition,

  pv_adj_threshold = pv_adj_threshold,
  
  dir_output = 'All'
)

# A549_All
sample_to_take <- metadata_df %>% rownames_to_column(var = 'sample') %>% filter(Line == 'A549') %>% pull(sample)
diff_expr_analysis(
  read_counts_df = read_counts_df %>% select(all_of(sample_to_take)),
  metadata_df = metadata_df %>% rownames_to_column(var = 'sample') %>% filter(sample %in% sample_to_take) %>% column_to_rownames(var = 'sample'),
  
  gene_df = ensembl_annotation_df,
  group_str = group_str,
  numerator_level_str = numerator_level_str,
  denominator_level_str = denominator_level_str,
  coef = coef,
  
  design_formula = ~ Batch + Condition,
  
  pv_adj_threshold = pv_adj_threshold,
  
  dir_output = 'A549_All'
)

# NHBE_WA1
sample_to_take <- metadata_df %>% rownames_to_column(var = 'sample') %>% filter(Line == 'NHBE') %>% pull(sample)
diff_expr_analysis(
  read_counts_df = read_counts_df %>% select(all_of(sample_to_take)),
  metadata_df = metadata_df %>% rownames_to_column(var = 'sample') %>% filter(sample %in% sample_to_take) %>% column_to_rownames(var = 'sample'),
  
  gene_df = ensembl_annotation_df,
  group_str = group_str,
  numerator_level_str = numerator_level_str,
  denominator_level_str = denominator_level_str,
  coef = coef,
  
  design_formula = ~ Condition,
  
  pv_adj_threshold = pv_adj_threshold,
  
  dir_output = 'NHBE_WA1'
)

# A549_H1N1
sample_to_take <- metadata_df %>% rownames_to_column(var = 'sample') %>% filter(Line == 'A549' & Batch == 'H1N1') %>% pull(sample)
diff_expr_analysis(
  read_counts_df = read_counts_df %>% select(all_of(sample_to_take)),
  metadata_df = metadata_df %>% rownames_to_column(var = 'sample') %>% filter(sample %in% sample_to_take) %>% column_to_rownames(var = 'sample'),
  
  gene_df = ensembl_annotation_df,
  group_str = group_str,
  numerator_level_str = numerator_level_str,
  denominator_level_str = denominator_level_str,
  coef = coef,
  
  design_formula = ~ Condition,
  
  pv_adj_threshold = pv_adj_threshold,
  
  dir_output = 'A549_H1N1'
)

# A549_WA1
sample_to_take <- metadata_df %>% rownames_to_column(var = 'sample') %>% filter(Line == 'A549' & Batch == 'WA1') %>% pull(sample)
diff_expr_analysis(
  read_counts_df = read_counts_df %>% select(all_of(sample_to_take)),
  metadata_df = metadata_df %>% rownames_to_column(var = 'sample') %>% filter(sample %in% sample_to_take) %>% column_to_rownames(var = 'sample'),
  
  gene_df = ensembl_annotation_df,
  group_str = group_str,
  numerator_level_str = numerator_level_str,
  denominator_level_str = denominator_level_str,
  coef = coef,
  
  design_formula = ~ Condition,
  
  pv_adj_threshold = pv_adj_threshold,
  
  dir_output = 'A549_WA1'
)

# A549_A2
sample_to_take <- metadata_df %>% rownames_to_column(var = 'sample') %>% filter(Line == 'A549' & Batch == 'A2') %>% pull(sample)
diff_expr_analysis(
  read_counts_df = read_counts_df %>% select(all_of(sample_to_take)),
  metadata_df = metadata_df %>% rownames_to_column(var = 'sample') %>% filter(sample %in% sample_to_take) %>% column_to_rownames(var = 'sample'),
  
  gene_df = ensembl_annotation_df,
  group_str = group_str,
  numerator_level_str = numerator_level_str,
  denominator_level_str = denominator_level_str,
  coef = coef,
  
  design_formula = ~ Condition,
  
  pv_adj_threshold = pv_adj_threshold,
  
  dir_output = 'A549_A2'
)

# A549andNHBE_WA1
sample_to_take <- metadata_df %>% rownames_to_column(var = 'sample') %>% filter(Batch == 'WA1') %>% pull(sample)
diff_expr_analysis(
  read_counts_df = read_counts_df %>% select(all_of(sample_to_take)),
  metadata_df = metadata_df %>% rownames_to_column(var = 'sample') %>% filter(sample %in% sample_to_take) %>% column_to_rownames(var = 'sample'),
  
  gene_df = ensembl_annotation_df,
  group_str = group_str,
  numerator_level_str = numerator_level_str,
  denominator_level_str = denominator_level_str,
  coef = coef,
  
  design_formula = ~ Line + Condition,
  
  pv_adj_threshold = pv_adj_threshold,
  
  dir_output = 'A549andNHBE_WA1'
)