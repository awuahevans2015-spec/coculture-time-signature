# scripts/04a_merge_open_targets.R -- run from the project root (Raw_Dat2)
fs <- list.files("data/genesets/open_targets_raw", pattern = "\\.tsv$", full.names = TRUE)
stopifnot(length(fs) >= 1)
tabs <- lapply(fs, function(f) {
  d <- read.delim(f, check.names = FALSE)
  stopifnot(all(c("symbol", "globalScore") %in% colnames(d)))   # catches wrong tables
  cat(basename(f), ":", nrow(d), "rows\n")
  data.frame(symbol = d$symbol, score = as.numeric(d$globalScore),
             disease = sub("\\.tsv$", "", basename(f)))
})
all <- do.call(rbind, tabs)
all <- all[!is.na(all$symbol) & all$symbol != "" & !is.na(all$score), ]
merged <- aggregate(score ~ symbol, data = all, FUN = max)
write.csv(merged, "data/genesets/open_targets_hf.csv", row.names = FALSE)
cat("Files used:", basename(fs), "\n")
cat("Genes:", nrow(merged), "| with score >= 0.30:", sum(merged$score >= 0.30), "\n")


d <- read.delim("data/genesets/open_targets_raw/heart_failure.tsv")
sum(duplicated(d$symbol)); sum(is.na(d$symbol) | d$symbol == ""); sum(is.na(d$globalScore))



# scripts/04_build_panel.R -- run from the project root (Raw_Dat2)
set.seed(42)
library(msigdbr); library(dplyr)
stopifnot(file.exists("data/processed/filtered_gene_background.csv"),
          file.exists("data/genesets/open_targets_hf.csv"))
bg <- read.csv("data/processed/filtered_gene_background.csv")[[1]]

get_sets <- function(coll, sub = NULL) {
  if (is.null(sub)) msigdbr(species = "Homo sapiens", collection = coll)
  else msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub)
}
kegg <- tryCatch(get_sets("C2", "CP:KEGG_LEGACY"), error = function(e) get_sets("C2", "CP:KEGG"))
all_sets <- bind_rows(get_sets("H"), get_sets("C2", "CP:REACTOME"), kegg, get_sets("C5", "GO:BP"))
wanted <- c("HALLMARK_ANGIOGENESIS", "HALLMARK_HYPOXIA", "HALLMARK_TGF_BETA_SIGNALING",
            "HALLMARK_NOTCH_SIGNALING", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "REACTOME_SIGNALING_BY_PDGF",
            "REACTOME_SIGNALING_BY_NOTCH", "REACTOME_SIGNALING_BY_VEGF", "REACTOME_COLLAGEN_FORMATION",
            "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS", "REACTOME_NON_INTEGRIN_MEMBRANE_ECM_INTERACTIONS",
            "KEGG_ECM_RECEPTOR_INTERACTION", "KEGG_FOCAL_ADHESION", "KEGG_NOTCH_SIGNALING_PATHWAY",
            "KEGG_VEGF_SIGNALING_PATHWAY", "GOBP_ANGIOGENESIS", "GOBP_CELL_CELL_JUNCTION_ORGANIZATION")
missing_sets <- setdiff(wanted, all_sets$gs_name)
cat("Sets not found:", missing_sets, "\n")
pathway_genes <- all_sets %>% filter(gs_name %in% wanted) %>% pull(gene_symbol) %>% unique()

mechano  <- c("YAP1", "WWTR1", "TEAD1", "TEAD4", "CCN1", "CCN2", "ANKRD1", "AMOTL2", "LATS1", "LATS2")
coupling <- c("PDGFB", "PDGFRB", "ANGPT1", "ANGPT2", "TEK", "NOTCH3", "JAG1", "DLL4", "COL4A1", "COL4A2",
              "LAMA4", "LAMB1", "LAMC1", "NID1", "TGFB1", "TGFBR2")
CUTOFF <- 0.30
ot <- read.csv("data/genesets/open_targets_hf.csv")
disease <- ot$symbol[ot$score >= CUTOFF]
markers <- c("PECAM1", "CDH5", "KDR", "PDGFRB", "CSPG4", "RGS5")
core <- intersect(unique(c(pathway_genes, mechano, coupling, disease)), bg)
core <- setdiff(core, markers)
with_markers <- union(core, intersect(markers, bg))

d <- format(Sys.Date(), "%Y-%m-%d")
writeLines(core, paste0("data/genesets/MVD_panel_core_", d, ".txt"))
writeLines(with_markers, paste0("data/genesets/MVD_panel_with_markers_", d, ".txt"))
writeLines(pathway_genes, paste0("data/genesets/pathway_genes_", d, ".txt"))
writeLines(intersect(disease, bg), paste0("data/genesets/disease_genes_", d, ".txt"))
cat("Panel sizes: core", length(core), "| with markers", length(with_markers),
    "| pathway genes", length(pathway_genes), "| disease genes in background",
    length(intersect(disease, bg)), "\n")

writeLines(c(paste("Panel built", d),
             paste("msigdbr version:", as.character(packageVersion("msigdbr"))),
             paste("Open Targets: UI v26.9.1, API 26.9.0; score cutoff", CUTOFF),
             "Open Targets diseases used: heart failure (MONDO_0005252); interstitial cardiac fibrosis (HP_0031329).",
             "Open Targets download date: EDIT",
             paste("MSigDB sets not found:", ifelse(length(missing_sets) == 0, "none", paste(missing_sets, collapse = ", ")))
), "docs/panel_provenance.md")
writeLines(capture.output(sessionInfo()), "docs/sessionInfo_04.txt")


cat("Sets not found:", missing_sets, "\n")
cat("Panel sizes: core", length(core), "| with markers", length(with_markers),
    "| pathway genes", length(pathway_genes), "| disease genes in background",
    length(intersect(disease, bg)), "\n")
list.files("data/genesets")

cat("Sets not found:", missing_sets, "\n")


sort(unique(grep("^GOBP_.*ANGIOGENESIS", all_sets$gs_name, value = TRUE)))
