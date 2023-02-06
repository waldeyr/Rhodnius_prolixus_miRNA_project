# Este script (maduros novos):
# 1) lê os resultados do mirDeep2 para cada amostra
# 2) filtra os resultados segundo alguns critérios usados em trabalhos anteriores
# 3) cria um diagrama de Venn
# 4) salva os resultados após filtro

library("readr")
library("dplyr")

# Takes this script dir as absolute path
basepath = dirname(rstudioapi::getSourceEditorContext()$path)

# SALIVARY GLAND
## Good ones
GLA_mature <- read_csv(file.path(basepath, "csv_files", "GLA_predicted_mature_novel.csv"))
GLA_mature['Sample'] = 'GLA'
GLA_mature <- GLA_mature %>%
  filter(`significant randfold p-value` == "yes") %>%
  filter(`miRDeep2 score` >= 1) %>%
  filter(`mature read count` >= 10) %>%
  filter(`star read count` >= 1)
write.table(
  GLA_mature, 
  file.path(basepath, "csv_files", "GLA_predicted_mature_novel_filtered.csv"), 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)

## Hypothetical ones
GLA_mature_hypothetical <- read_csv(file.path(basepath, "csv_files", "GLA_predicted_mature_novel.csv"))
GLA_mature_hypothetical['Sample'] = 'GLA'
GLA_mature_hypothetical <- GLA_mature_hypothetical %>%
  filter(`significant randfold p-value` != "yes") %>%
  filter(`miRDeep2 score` < 1) %>%
  filter(`mature read count` < 10) %>%
  filter(`star read count` < 1)
write.table(
  GLA_mature_hypothetical, 
  file.path(basepath, "csv_files", "GLA_predicted_mature_novel_filtered_hypothetical.csv"), 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)

# GUT
## Good ones
GUT_mature <- read_csv(file.path(basepath, "csv_files", "GUT_predicted_mature_novel.csv"))
GUT_mature['Sample'] = 'GUT'
GUT_mature <- GUT_mature %>%
  filter(`significant randfold p-value` == "yes") %>%
  filter(`miRDeep2 score` >= 1) %>%
  filter(`mature read count` >= 10) %>%
  filter(`star read count` >= 1)
write.table(
  GUT_mature, 
  file.path(basepath, "csv_files", "GUT_predicted_mature_novel_filtered.csv"), 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)

## Hypothetical ones
GUT_mature_hypothetical <- read_csv(file.path(basepath, "csv_files", "GUT_predicted_mature_novel.csv"))
GUT_mature_hypothetical['Sample'] = 'GUT'
GUT_mature_hypothetical <- GUT_mature_hypothetical %>%
  filter(`significant randfold p-value` != "yes") %>%
  filter(`miRDeep2 score` < 1) %>%
  filter(`mature read count` < 10) %>%
  filter(`star read count` < 1)
write.table(
  GUT_mature_hypothetical, 
  file.path(basepath, "csv_files", "GUT_predicted_mature_novel_filtered_hypothetical.csv"), 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)


# HEMOLIMPH
## Good ones
HEM_mature <- read_csv(file.path(basepath, "csv_files", "HEM_predicted_mature_novel.csv"))
HEM_mature['Sample'] = 'HEM'
HEM_mature <- HEM_mature %>%
  filter(`significant randfold p-value` == "yes") %>%
  filter(`miRDeep2 score` >= 1) %>%
  filter(`mature read count` >= 10) %>%
  filter(`star read count` >= 1)
write.table(
  HEM_mature, 
  file.path(basepath, "csv_files", "HEM_predicted_mature_novel_filtered.csv"), 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)

## Hypothetical ones
HEM_mature_hypothetical <- read_csv(file.path(basepath, "csv_files", "HEM_predicted_mature_novel.csv"))
HEM_mature_hypothetical['Sample'] = 'HEM'
HEM_mature_hypothetical <- HEM_mature_hypothetical %>%
  filter(`significant randfold p-value` != "yes") %>%
  filter(`miRDeep2 score` < 1) %>%
  filter(`mature read count` < 10) %>%
  filter(`star read count` < 1)
write.table(
  HEM_mature_hypothetical, 
  file.path(basepath, "csv_files", "HEM_predicted_mature_novel_filtered_hypothetical.csv"), 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)


nrow(GLA_mature)
nrow(GUT_mature)
nrow(HEM_mature)