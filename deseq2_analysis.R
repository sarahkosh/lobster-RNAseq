## Differential expression of lobster gill and muscle RNA-seq (DESeq2)
##
## For each tissue, compares control against each treatment (warming,
## acidification, multistressor).
##
## Input:  per-sample featureCounts files (quants_trimmed/featurecounts_<id>.txt)
##         metadata.csv with columns: sample, treatment, tissue_type
## Output: results/ (significant-gene tables, volcano plots, PCA plots)
##
## Fold-change direction: contrasts are control vs. treatment, so a positive
## log2 fold change means HIGHER expression in the control.

library(DESeq2)
library(apeglm)
library(ggplot2)
library(readr)

## ---- Settings ----
quant_dir <- "quants_trimmed"
meta_file <- "metadata.csv"
out_dir   <- "results"
dir.create(out_dir, showWarnings = FALSE)

padj_cutoff <- 0.05
lfc_cutoff  <- 1

# Minimum total count per gene (summed across samples), applied before testing.
# Thresholds differ by tissue, as in the original analysis.
min_total_count <- c(gill = 10, muscle = 5)

sample_ids <- list(
  gill   = c("g003", "g005", "g007", "g009", "g014", "g020", "g022", "g032"),
  muscle = c("m003", "m005", "m007", "m009", "m014", "m015", "m020", "m022", "m032")
)

# Treatment -> suffix used in output file names
treatments <- c(warming = "warm", acidification = "acid", multistressor = "multi")

vals5 <- c(warming = "#55b0db", acidification = "#bf389f",
           multistressor = "#7272cc", control = "#519c3d")

## ---- Helpers ----
# Read per-sample featureCounts files into one gene x sample count matrix
read_counts <- function(ids) {
  per_sample <- lapply(ids, function(id) {
    f  <- read.table(file.path(quant_dir, paste0("featurecounts_", id, ".txt")),
                     header = TRUE)
    df <- data.frame(Geneid = f$Geneid, count = f[[7]])  # column 7 = read counts
    names(df)[2] <- id
    df
  })
  merged <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), per_sample)
  stopifnot(!anyNA(merged))
  rownames(merged) <- gsub("-", ".", merged$Geneid)
  merged$Geneid <- NULL
  merged
}

meta <- as.data.frame(read_csv(meta_file, show_col_types = FALSE))
rownames(meta) <- meta$sample

# Metadata rows in the same order as the count columns
get_coldata <- function(ids) {
  stopifnot(all(ids %in% meta$sample))
  coldata <- meta[ids, ]
  coldata$treatment <- relevel(factor(coldata$treatment), ref = "control")
  coldata
}

run_deseq <- function(tissue) {
  ids     <- sample_ids[[tissue]]
  counts  <- read_counts(ids)
  coldata <- get_coldata(ids)
  stopifnot(identical(colnames(counts), rownames(coldata)))
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                                design = ~ treatment)
  dds <- dds[rowSums(counts(dds)) >= min_total_count[[tissue]], ]
  DESeq(dds)
}

# Results, significant-gene tables, and volcano plots for each treatment
export_results <- function(dds, tissue) {
  prefix <- substr(tissue, 1, 1)
  for (trt in names(treatments)) {
    message("\n", tissue, ": control vs. ", trt)
    res <- results(dds, contrast = c("treatment", "control", trt))
    print(summary(res))

    sig <- subset(res, padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff)
    write.csv(sig, file.path(out_dir,
              sprintf("%s_sig_genes_ctrlV%s.csv", prefix, treatments[[trt]])))

    p <- ggplot(as.data.frame(sig), aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point() +
      geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "red") +
      geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "blue") +
      labs(x = "Log2 Fold Change", y = "-log10(adjusted p-value)") +
      theme_minimal() +
      ggtitle(paste0("Control vs. ", trt, " - ", tissue))
    ggsave(file.path(out_dir, sprintf("%s_volcano_ctrlV%s.png", prefix, treatments[[trt]])),
           p, width = 5, height = 4)
  }
}

## ---- Differential expression ----
m_dds <- run_deseq("muscle")
g_dds <- run_deseq("gill")

export_results(m_dds, "muscle")
export_results(g_dds, "gill")

# Optional: shrunken log2 fold changes for muscle, warming. Note the sign:
# this coefficient is warming vs. control, the opposite of the contrasts above.
res_shrunk <- lfcShrink(m_dds, coef = "treatment_warming_vs_control", type = "apeglm")
png(file.path(out_dir, "m_MA_shrunk_warming.png"), width = 600, height = 500)
plotMA(res_shrunk, ylim = c(-2, 2), main = "Muscle: warming vs. control (shrunken)")
dev.off()

## ---- PCA (rlog-transformed counts) ----
pca_plot <- function(rld, title, subtitle = NULL, ntop = 5000) {
  d   <- plotPCA(rld, intgroup = "treatment", ntop = ntop, returnData = TRUE)
  pct <- round(100 * attr(d, "percentVar"))
  d$tissue <- ifelse(substr(d$name, 1, 1) == "m", "muscle", "gill")
  ggplot(d, aes(x = PC1, y = PC2, color = treatment, shape = tissue)) +
    geom_point(alpha = 0.9, size = 5) +
    scale_color_manual(values = vals5, name = "") +
    labs(x = paste0("PC1: ", pct[1], "% variance"),
         y = paste0("PC2: ", pct[2], "% variance"),
         title = title, subtitle = subtitle) +
    theme_minimal()
}

m_rld <- rlog(m_dds, blind = FALSE)
g_rld <- rlog(g_dds, blind = FALSE)

ggsave(file.path(out_dir, "m_PCA.png"), pca_plot(m_rld, "PCA", "Muscle"), width = 6, height = 5)
ggsave(file.path(out_dir, "g_PCA.png"), pca_plot(g_rld, "PCA", "Gill"),   width = 6, height = 5)

# PCA with m009, m014, and m015 removed. For visualization only: these samples
# are still included in the differential expression analysis above.
m_keep <- c("m003", "m005", "m007", "m020", "m022", "m032")
ggsave(file.path(out_dir, "m_PCA_no_009_014_015.png"),
       pca_plot(m_rld[, m_rld$sample %in% m_keep], "PCA",
                "Muscle without m009, m014, m015"), width = 6, height = 5)

# Both tissues together
all_ids <- c(sample_ids$muscle, sample_ids$gill)
all_dds <- DESeqDataSetFromMatrix(countData = read_counts(all_ids),
                                  colData = get_coldata(all_ids), design = ~ 1)
all_dds <- all_dds[rowSums(counts(all_dds)) >= 10, ]
all_rld <- rlog(DESeq(all_dds), blind = FALSE)

ggsave(file.path(out_dir, "all_PCA.png"),
       pca_plot(all_rld, "PCA", "All tissue types"), width = 6, height = 5)

# Subset with the 009, 014, and 015 samples removed (visualization only)
all_keep <- c(m_keep, "g003", "g005", "g007", "g020", "g022", "g032")
ggsave(file.path(out_dir, "all_PCA_subset.png"),
       pca_plot(all_rld[, all_rld$sample %in% all_keep], "PCA",
                "All tissue types, subset of samples"), width = 6, height = 5)
