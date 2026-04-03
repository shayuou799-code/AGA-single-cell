# ============================================================================
# Statistical Analysis
# Project: AGA single-cell RNA-seq (Nature Communications)
# Description: Statistical tests used for figure panels and organ culture
#              experiments (complementing GraphPad Prism analyses)
# ============================================================================

set.seed(42)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

# ============================================================================
# Part 1: Wilcoxon Signed-Rank Test (Paired Patient-Level Comparisons)
# Used for: Fig 1i (immune cell proportions), module score comparisons
# ============================================================================

# --- 1. Calculate per-patient cell proportions ---
# Example: Mast cell proportion per patient (paired Vertex vs Occipital)
cell_counts <- AGA@meta.data %>%
  group_by(PatientID, Type, level1_identity) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(PatientID, Type) %>%
  mutate(total = sum(n),
         proportion = n / total) %>%
  ungroup()

mast_prop <- cell_counts %>%
  filter(level1_identity == "Mast")

# --- 2. Paired Wilcoxon signed-rank test ---
# Paired because Vertex and Occipital come from the same patient
vertex_prop  <- mast_prop %>% filter(Type == "Vertex") %>% arrange(PatientID)
occ_prop     <- mast_prop %>% filter(Type == "Occipital") %>% arrange(PatientID)

wilcox_result <- wilcox.test(vertex_prop$proportion,
                              occ_prop$proportion,
                              paired = TRUE,
                              alternative = "two.sided")
print(wilcox_result)

# --- 3. Paired boxplot visualization ---
ggpaired(mast_prop, x = "Type", y = "proportion",
         id = "PatientID",
         color = "Type", line.color = "gray70",
         palette = c("Occipital" = "#4DBBD5", "Vertex" = "#E64B35")) +
  stat_compare_means(method = "wilcox.test", paired = TRUE) +
  ylab("Mast Cell Proportion") +
  ggtitle("Paired Mast Cell Proportion: Vertex vs Occipital")

# ============================================================================
# Part 2: Module Score Comparison (Wilcoxon Rank-Sum with Bonferroni)
# Used for: Fig 1d (inflammation scores), Fig 2a (MP2/MP9 scores)
# ============================================================================

# --- 4. Example: Compare inflammation score between regions ---
# AddModuleScore was already run upstream; scores are in metadata
inflammation_data <- AGA@meta.data %>%
  select(level1_identity, Type, chronic_inflammation_score)

# Wilcoxon rank-sum test per cell type with Bonferroni correction
stat_results <- inflammation_data %>%
  group_by(level1_identity) %>%
  wilcox_test(chronic_inflammation_score ~ Type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

print(stat_results)

# ============================================================================
# Part 3: Spearman Correlation
# Used for: Mast cell abundance vs HFSC stress scores
# ============================================================================

# --- 5. Patient-level aggregation for correlation ---
patient_summary <- AGA@meta.data %>%
  group_by(PatientID, Type) %>%
  summarise(
    mast_proportion = sum(level1_identity == "Mast") / n(),
    mean_hfsc_stress = mean(HFSC_stress_score[level1_identity %in%
                              c("HFSC1", "HFSC2")], na.rm = TRUE),
    .groups = "drop"
  )

# Spearman correlation
cor_test <- cor.test(patient_summary$mast_proportion,
                     patient_summary$mean_hfsc_stress,
                     method = "spearman")
print(cor_test)

ggscatter(patient_summary,
          x = "mast_proportion", y = "mean_hfsc_stress",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Mast Cell Proportion",
          ylab = "Mean HFSC Stress Score")

# ============================================================================
# Part 4: Ex Vivo Organ Culture Statistics
# Used for: Fig 2j (TNF-α shaft elongation), Fig 6h (AREG shaft elongation)
# Note: These analyses were also performed in GraphPad Prism (v9.0)
# ============================================================================

# --- 6. Two-tailed Student's t-test (hair shaft elongation) ---
# Example data structure (replace with actual measurements):
# control_elongation <- c(...)  # mm, n=12 follicles
# tnf_elongation     <- c(...)  # mm, n=12 follicles

t_test_result <- t.test(control_elongation, tnf_elongation,
                         paired = FALSE,
                         alternative = "two.sided",
                         var.equal = FALSE)  # Welch's t-test
print(t_test_result)

# --- 7. Fisher's exact test (hair cycle stage proportions) ---
# Example: Anagen vs Catagen counts in control vs TNF-treated groups
# Used for: Fig 2k, Fig 6i
cycle_table <- matrix(c(
  10, 2,   # Control: anagen=10, catagen=2
   5, 7    # TNF-treated: anagen=5, catagen=7
), nrow = 2, byrow = TRUE,
dimnames = list(
  Group = c("Control", "TNF"),
  Stage = c("Anagen", "Catagen")
))

fisher_result <- fisher.test(cycle_table)
print(fisher_result)
# Reports odds ratio and 95% CI

# ============================================================================
# Part 5: Mouse Model Statistics
# Used for: Fig 3d-g quantification
# ============================================================================

# --- 8. Mast cell count comparison (Toluidine blue quantification) ---
# Two-tailed Student's t-test
# control_counts <- c(...)  # mast cells per HPF, n = [mice per group]
# model_counts   <- c(...)

t_test_mast <- t.test(control_counts, model_counts,
                       paired = FALSE,
                       alternative = "two.sided")
print(t_test_mast)

# --- 9. IF quantification (fluorescence intensity) ---
# MFI comparison: TNF-α, AREG, c-FOS, Ki67
# t_test_mfi <- t.test(control_mfi, model_mfi, ...)

# ============================================================================
# Summary of statistical methods
# ============================================================================
# | Test                        | Application                    | Figure |
# |-----------------------------|--------------------------------|--------|
# | Wilcoxon rank-sum + Bonf.   | DEG, module score comparisons  | 1d,2g  |
# | Wilcoxon signed-rank        | Paired V vs O proportions      | 1i     |
# | Spearman correlation        | Mast abundance vs HFSC stress  | --     |
# | Two-tailed Student's t-test | Organ culture, mouse model     | 2j,6h  |
# | Fisher's exact test         | Hair cycle stage proportions   | 2k,6i  |
# ============================================================================

# sessionInfo()  # Uncomment to print session information
