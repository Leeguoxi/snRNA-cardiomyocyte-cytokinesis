library(Seurat)
library(dplyr)
library(clusterProfiler)

#############Using a stepwise filtering strategy based on cell cycle status, developmental stage, and cell subpopulations, we identified cardiomyocyte subgroups that are likely to complete cytokinesis after entering the cell cycle (Cytokinesis+), likely to fail cytokinesis after entering the cell cycle (Cytokinesis−), and remain cell-cycle silent (Silent).
mouse_cell_cycle_genes <- readRDS("mouse_cell_cycle_genes.rds")  # Load predefined mouse cell-cycle gene sets
cm<- readRDS("cm.rds")                                            # Load Seurat object

s.genes <- mouse_cell_cycle_genes$s.genes                         # Extract S-phase gene list
g2m.genes <- mouse_cell_cycle_genes$g2m.genes                     # Extract G2/M-phase gene list
cm<- CellCycleScoring(                                             # Score each cell for cell-cycle phase
  object = cm,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE                                                # Do not overwrite current identities
)

cc_slect <- subset(cm, subset = Phase%in% c("G2M","G1"))          # Keep only cells in G2M or G1 phase

cc_age_slect<- subset(cc_slect, subset =age %in%  c("E18.5","P4","P5","P6","P14","P28"))  # Keep selected developmental ages

cc_age_celltype_slect<-subset(cc_age_slect, subset = celltype %in% c("a","c","e"))         # Keep only celltypes a, c, and e

cm_new=cc_age_celltype_slect                                       # Create a new working object
cm_new=RenameIdents(cm_new,"a"="Cytokinesis +")                   # Rename identity a
cm_new=RenameIdents(cm_new,"c"="Cytokinesis -")                   # Rename identity c
cm_new=RenameIdents(cm_new,"e"="Silent")                          # Rename identity e
Idents(cm_new)=factor(Idents(cm_new),levels = c("Cytokinesis +","Cytokinesis -","Silent")) # Set identity order

########Explore gene sets that are progressively upregulated and downregulated across these three cell types.

# 1) Differential expression: Cytokinesis - vs Cytokinesis +
m1 <- FindMarkers(
  cm_new,
  ident.1 = "Cytokinesis -",
  ident.2 = "Cytokinesis +",
  logfc.threshold = 0,   # keep all genes regardless of minimum logFC
  min.pct = 0.25,        # gene must be expressed in at least 25% of cells in either group
  test.use = "wilcox"    # Wilcoxon rank-sum test
)

# 2) Differential expression: Silent vs Cytokinesis -
m2 <- FindMarkers(
  cm_new,
  ident.1 = "Silent",
  ident.2 = "Cytokinesis -",
  logfc.threshold = 0,
  min.pct = 0.25,
  test.use = "wilcox"
)

# 3) Differential expression: Silent vs Cytokinesis +
m3 <- FindMarkers(
  cm_new,
  ident.1 = "Silent",
  ident.2 = "Cytokinesis +",
  logfc.threshold = 0,
  min.pct = 0.25,
  test.use = "wilcox"
)

# 4) Apply Benjamini-Hochberg correction to raw p-values for each comparison
m1$FDR_BH <- p.adjust(m1$p_val, method = "BH")
m2$FDR_BH <- p.adjust(m2$p_val, method = "BH")
m3$FDR_BH <- p.adjust(m3$p_val, method = "BH")

# 5) Keep significantly upregulated genes in each comparison
#    fc_col should be the logFC column name (e.g., "avg_log2FC" in Seurat v4/v5)
up1 <- rownames(m1)[m1$FDR_BH < 0.05 & m1[[fc_col]] > 0]
up2 <- rownames(m2)[m2$FDR_BH < 0.05 & m2[[fc_col]] > 0]
up3 <- rownames(m3)[m3$FDR_BH < 0.05 & m3[[fc_col]] > 0]

# 6) Candidate trend genes = intersection of upregulated sets from all three comparisons
genes_trend <- Reduce(intersect, list(up1, up2, up3))

# 7) Compute average expression per group for candidate genes
#    (use genes present in cm_new RNA feature set)
avg_exp <- AverageExpression(
  cm_new,
  features = intersect(genes_trend, rownames(cm_new)),
  group.by = "celltype_fix",   # grouping metadata column
  assays = "RNA",
  slot = "data"                # log-normalized expression
)$RNA

# 8) Ensure columns are ordered as progression: Cytokinesis + -> Cytokinesis - -> Silent
avg_exp <- avg_exp[, c("Cytokinesis +", "Cytokinesis -", "Silent"), drop = FALSE]

# 9) Enforce minimum stepwise increase between adjacent states
delta <- 0.10
genes_mono_delta <- rownames(avg_exp)[
  (avg_exp[, 2] - avg_exp[, 1]) > delta &
  (avg_exp[, 3] - avg_exp[, 2]) > delta
]







