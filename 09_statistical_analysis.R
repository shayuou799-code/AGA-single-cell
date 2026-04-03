# ============================================================================
set.seed(42)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

# ============================================================================
# Part 1: Wilcoxon Signed-Rank Test (Paired Patient-Level Comparisons)
# ============================================================================

#  Paired Wilcoxon signed-rank test ---
# Paired because Vertex and Occipital come from the same patient
vertex_prop  <- mast_prop %>% filter(Type == "Vertex") %>% arrange(PatientID)
occ_prop     <- mast_prop %>% filter(Type == "Occipital") %>% arrange(PatientID)

wilcox_result <- wilcox.test(vertex_prop$proportion,
                              occ_prop$proportion,
                              paired = TRUE,
                              alternative = "two.sided")
print(wilcox_result)

#  Paired boxplot visualization ---
ggpaired(mast_prop, x = "Type", y = "proportion",
         id = "PatientID",
         color = "Type", line.color = "gray70",
         palette = c("Occipital" = "#4DBBD5", "Vertex" = "#E64B35")) +
  stat_compare_means(method = "wilcox.test", paired = TRUE) +
  ylab("Mast Cell Proportion") +
  ggtitle("Paired Mast Cell Proportion: Vertex vs Occipital")

# ============================================================================
# Part 2: Module Score Comparison (Wilcoxon Rank-Sum with Bonferroni)
# ============================================================================

#   Compare inflammation score between regions ---
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
# Spearman Correlation
# ============================================================================

# Patient-level aggregation for correlation ---
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


# sessionInfo()  # Uncomment to print session information
