#!/usr/bin/env Rscript
# Fig 4c reproduction: ASD variant prioritization curves
# Uses: Table S2 (An et al. 2018), asd.FLARE.txt, ASD.FLARE-fb.txt, SFARI gene list

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(readxl)
})

base      <- "."
sfari_f   <- "./SFARI-Gene_genes_05-01-2026release_06-27-2026export.csv"
tables2_f <- file.path(base, "aat6576_table-s2.xlsx")
flare_f   <- file.path(base, "asd.FLARE.txt")
annot_f   <- file.path(base, "ASD.FLARE-fb.txt")
out_pdf   <- file.path(base, "fig4c_multi.pdf")

cat("Loading Table S2 (An et al. 2018)...\n")
s2_raw <- read_excel(tables2_f, sheet = 1, skip = 1, col_names = TRUE)
# Row 1 after skip is the real header
names(s2_raw)[1:13] <- c("Chr","Pos","Ref","Alt","Type","Fam","Pheno",
                          "SampleID","Consequence","SYMBOL","Gene","DISTANCE","NEAREST")
s2 <- as.data.frame(s2_raw[, c("Chr","Pos","Ref","Alt","Pheno","SampleID","SYMBOL","NEAREST")])
s2$snp_id     <- paste0(s2$Chr, ":", s2$Pos, ":", s2$Ref, ":", s2$Alt)
s2$GeneSymbol <- ifelse(s2$SYMBOL != ".", s2$SYMBOL, s2$NEAREST)
s2$Pheno      <- factor(s2$Pheno, levels = c("control", "case"))
cat("Table S2:", nrow(s2), "variants\n")

cat("Loading FLARE scores...\n")
flare <- fread(flare_f, data.table = FALSE)
cat("FLARE:", nrow(flare), "variants\n")

cat("Loading feature annotations...\n")
annot_cols <- c("snp_id", "gene_distance_1.log10", "abs_logfc.mean.c11.trevino_2021")
annot <- fread(annot_f, data.table = FALSE, select = annot_cols)

cat("Loading SFARI gene list...\n")
sfari     <- fread(sfari_f, data.table = FALSE)
sfari_syn <- subset(sfari, syndromic == 1 & `number-of-reports` >= 10)[, "gene-symbol"]
cat("Syndromic SFARI genes (>=10 reports):", length(sfari_syn), "\n")

# Merge everything on snp_id
cat("Merging datasets...\n")
df <- merge(flare, s2[, c("snp_id","Pheno","GeneSymbol","SampleID")], by = "snp_id")
df <- merge(df,   annot, by = "snp_id")
cat("After merge:", nrow(df), "variants\n")

# SFARI filter: keep variants near syndromic genes
df$sfari <- df$GeneSymbol %in% sfari_syn
tmp <- subset(df, sfari & !is.na(Pheno))
cat("Near syndromic SFARI genes:", nrow(tmp),
    "| cases:", sum(tmp$Pheno == "case"),
    "| controls:", sum(tmp$Pheno == "control"), "\n")

bg <- mean(tmp$Pheno == "case")
cat("Background case proportion:", round(bg, 6), "\n")

# Outlier prioritization curve
outlier_curve <- function(score, data, decreasing = TRUE) {
  s   <- data[order(data[[score]], decreasing = decreasing), ]
  rng <- 3:min(1000, nrow(s))
  rbindlist(lapply(rng, function(j) {
    tab <- table(factor(s$Pheno[1:j], levels = c("control", "case")))
    bt  <- binom.test(tab["case"], j, bg)
    data.table(n    = j,
               prop = as.numeric(tab["case"] / j),
               l    = bt$conf.int[1],
               h    = bt$conf.int[2])
  }))
}

predictors <- c("FLARE_fb", "FLARE_heart",
                "abs_logfc.mean.c11.trevino_2021",
                "phylop", "gene_distance_1.log10")
pred_present <- intersect(predictors, names(tmp))
cat("Predictors available:", paste(pred_present, collapse = ", "), "\n")

cat("Computing prioritization curves...\n")
res_lst <- lapply(pred_present, function(p) {
  dec <- (p != "gene_distance_1.log10")  # gene_distance: higher = closer, still sort desc
  r   <- outlier_curve(p, tmp, decreasing = dec)
  r$predictor <- p
  r
})
res <- rbindlist(res_lst)

# Ribbon data: CI only for FLARE_fb
ci_data <- res[res$predictor == "FLARE_fb", ]

color_map <- c(
  FLARE_fb                         = "#E0CA70",
  FLARE_heart                      = "#B30606",
  abs_logfc.mean.c11.trevino_2021  = "#555599",
  phylop                           = "#3CC2B2",
  gene_distance_1.log10            = "black"
)
label_map <- c(
  FLARE_fb                         = "FLARE: fetal brain",
  FLARE_heart                      = "FLARE: heart",
  abs_logfc.mean.c11.trevino_2021  = "ChromBPNet: early RG",
  phylop                           = "PhyloP",
  gene_distance_1.log10            = expression(log[10]*"(TSS distance)")
)

g <- ggplot(res, aes(x = log10(n), y = prop, col = predictor)) +
  geom_ribbon(data    = ci_data,
              mapping = aes(x = log10(n), ymin = l, ymax = h),
              fill    = "#E3E0E0", alpha = 0.4, inherit.aes = FALSE) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = bg, col = "red", linetype = "dashed") +
  geom_vline(xintercept = log10(16), linetype = "dotted") +
  scale_x_continuous(breaks = 1:3, labels = c("10", "100", "1000")) +
  scale_color_manual(values = color_map[pred_present],
                     labels = label_map[pred_present]) +
  labs(x   = "Outlier variant rank\n(near syndromic ASD genes)",
       y   = "Proportion of case mutations",
       col = NULL) +
  theme_classic() +
  theme(plot.title       = element_text(hjust = 0.5),
        legend.position  = "right",
        legend.text.align = 0)

pdf(out_pdf, width = 5.4 * 1.3, height = 2.96 * 1.3)
print(g)
dev.off()
cat("Wrote:", out_pdf, "\n")
cat("Top 16 case fraction (FLARE_fb):", sum(head(tmp[order(tmp$FLARE_fb, decreasing=TRUE),]$Pheno, 16)=="case"), "/ 16\n")
