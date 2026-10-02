# Differential expression with DESeq2 for every contrast defined in config.yaml.
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

cfg <- snakemake@params[["diffexp"]]
contrast_names <- snakemake@params[["contrasts"]]

counts <- as.matrix(read.delim(snakemake@input[["counts"]], row.names = 1, check.names = FALSE))
samples <- read.delim(snakemake@input[["samples"]], colClasses = "character", check.names = FALSE)
rownames(samples) <- samples$sample
samples <- samples[colnames(counts), , drop = FALSE]

design <- as.formula(cfg$design)
for (v in all.vars(design)) {
  if (!v %in% names(samples)) stop("Design variable '", v, "' is not a column of the sample sheet")
  samples[[v]] <- factor(samples[[v]], levels = unique(samples[[v]]))
}

dds <- DESeqDataSetFromMatrix(countData = counts, colData = samples, design = design)
min_count <- if (is.null(cfg$min_count)) 10 else cfg$min_count
min_samples <- if (is.null(cfg$min_samples)) 1 else min(cfg$min_samples, ncol(dds))
keep <- rowSums(counts(dds) >= min_count) >= min_samples
message(sprintf("Pre-filtering: %d / %d genes kept", sum(keep), nrow(dds)))
dds <- DESeq(dds[keep, ])
saveRDS(dds, snakemake@output[["rds"]])

norm <- counts(dds, normalized = TRUE)
write.table(cbind(gene_id = rownames(norm), as.data.frame(round(norm, 3))),
            snakemake@output[["normalized"]], sep = "\t", quote = FALSE, row.names = FALSE)

# --- Sample-level QC -------------------------------------------------------
vsd <- if (nrow(dds) >= 1000) vst(dds, blind = TRUE) else varianceStabilizingTransformation(dds, blind = TRUE)
first_var <- all.vars(design)[length(all.vars(design))]
pca <- plotPCA(vsd, intgroup = first_var, returnData = TRUE)
pv <- round(100 * attr(pca, "percentVar"), 1)
p <- ggplot(pca, aes(PC1, PC2, colour = .data[[first_var]], label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  labs(x = sprintf("PC1: %s%% variance", pv[1]), y = sprintf("PC2: %s%% variance", pv[2])) +
  theme_bw()
ggsave(snakemake@output[["pca"]], p, width = 6, height = 5, dpi = 150)

d <- as.matrix(dist(t(assay(vsd))))
pheatmap(d, clustering_distance_rows = dist(t(assay(vsd))),
         clustering_distance_cols = dist(t(assay(vsd))),
         annotation_col = samples[, first_var, drop = FALSE],
         main = "Sample-to-sample distances (VST)",
         filename = snakemake@output[["heatmap"]], width = 6, height = 5)

# --- Contrasts -----------------------------------------------------------
colors <- c(Up = "#C0392B", Down = "#2471A3", "Not significant" = "grey70")
for (i in seq_along(contrast_names)) {
  name <- contrast_names[[i]]
  ct <- unlist(cfg$contrasts[[name]])
  res <- results(dds, contrast = ct, alpha = cfg$alpha)
  shr <- lfcShrink(dds, contrast = ct, res = res, type = "normal")
  df <- data.frame(gene_id = rownames(res), baseMean = res$baseMean,
                   log2FoldChange = res$log2FoldChange,
                   log2FoldChange_shrunk = shr$log2FoldChange,
                   lfcSE = res$lfcSE, stat = res$stat,
                   pvalue = res$pvalue, padj = res$padj)
  df$status <- ifelse(!is.na(df$padj) & df$padj < cfg$alpha & df$log2FoldChange_shrunk > cfg$lfc_threshold, "Up",
               ifelse(!is.na(df$padj) & df$padj < cfg$alpha & df$log2FoldChange_shrunk < -cfg$lfc_threshold, "Down",
                      "Not significant"))
  df <- df[order(df$padj, na.last = TRUE), ]
  write.table(df, snakemake@output[["results"]][[i]], sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("%s: %d up, %d down (padj < %s, |shrunken LFC| > %s)", name,
                  sum(df$status == "Up"), sum(df$status == "Down"), cfg$alpha, cfg$lfc_threshold))

  title <- sprintf("%s vs %s", ct[2], ct[3])
  v <- df[!is.na(df$padj), ]
  pv <- ggplot(v, aes(log2FoldChange_shrunk, -log10(padj), colour = status)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_colour_manual(values = colors) +
    geom_vline(xintercept = c(-1, 1) * cfg$lfc_threshold, linetype = "dashed") +
    geom_hline(yintercept = -log10(cfg$alpha), linetype = "dashed") +
    labs(title = paste("Volcano plot:", title), x = "log2 fold change (shrunken)",
         y = "-log10 adjusted p-value", colour = NULL) +
    theme_bw()
  ggsave(snakemake@output[["volcano"]][[i]], pv, width = 6, height = 5, dpi = 150)

  pm <- ggplot(df[df$baseMean > 0, ], aes(baseMean, log2FoldChange_shrunk, colour = status)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_x_log10() +
    scale_colour_manual(values = colors) +
    geom_hline(yintercept = 0) +
    labs(title = paste("MA plot:", title), x = "Mean of normalised counts",
         y = "log2 fold change (shrunken)", colour = NULL) +
    theme_bw()
  ggsave(snakemake@output[["ma"]][[i]], pm, width = 6, height = 5, dpi = 150)
}

writeLines(capture.output(sessionInfo()))
