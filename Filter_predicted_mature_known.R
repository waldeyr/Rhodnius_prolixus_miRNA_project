# Este script (maduros conhecidos):
# 1) lê os resultados do mirDeep2 para cada amostra
# 2) filtra os resultados segundo alguns critérios usados em trabalhos anteriores
# 3) cria um diagrama de Venn
# 4) salva os resultados após filtro

library("readr")
library("dplyr")

# Takes this script dir as absolute path
basepath = dirname(rstudioapi::getSourceEditorContext()$path)
# SALIVARY GLAND
## GOOD ONES
GLA_mature <- read_csv(file.path(basepath, "csv_files", "GLA_predicted_mature_known.csv"), show_col_types = FALSE)
GLA_mature['Sample'] = 'GLA'
GLA_mature <- GLA_mature %>%
  filter(`significant randfold p-value` == "yes") %>%
  filter(`miRDeep2 score` >= 1) %>%
  filter(`mature read count` >= 10) %>%
  filter(`star read count` >= 1)
write.table(
  GLA_mature, 
  file.path(basepath, "csv_files", "GLA_predicted_mature_known_filtered.csv"), 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)

## HYPOTHETICAL ONES
GLA_mature_hypothetical <- read_csv(file.path(basepath, "csv_files", "GLA_predicted_mature_known.csv"), show_col_types = FALSE)
GLA_mature_hypothetical['Sample'] = 'GLA'
GLA_mature_hypothetical <- GLA_mature_hypothetical %>%
  filter(`significant randfold p-value` != "yes") %>%
  filter(`miRDeep2 score` < 1) %>%
  filter(`mature read count` < 10) %>%
  filter(`star read count` < 1)
write.table(
  GLA_mature_hypothetical, 
  file.path(basepath, "csv_files", "GLA_predicted_mature_known_filtered_hypothetical.csv"), 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)

# GUT
## GOOD ONES
GUT_mature <- read_csv(file.path(basepath, "csv_files", "GUT_predicted_mature_known.csv"), show_col_types = FALSE)
GUT_mature['Sample'] = 'GUT'
GUT_mature <- GUT_mature %>%
  filter(`significant randfold p-value` == "yes") %>%
  filter(`miRDeep2 score` >= 1) %>%
  filter(`mature read count` >= 10) %>%
  filter(`star read count` >= 1)
write.table(
  GUT_mature, 
  file.path(basepath, "csv_files", "GUT_predicted_mature_known_filtered.csv"), 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)
## HYPOTHETICAL ONES
GUT_mature_hypothetical <- read_csv(file.path(basepath, "csv_files", "GUT_predicted_mature_known.csv"), show_col_types = FALSE)
GUT_mature_hypothetical['Sample'] = 'GUT'
GUT_mature_hypothetical <- GUT_mature_hypothetical %>%
  filter(`significant randfold p-value` != "yes") %>%
  filter(`miRDeep2 score` < 1) %>%
  filter(`mature read count` < 10) %>%
  filter(`star read count` < 1 )
write.table(
  GUT_mature_hypothetical, 
  file.path(basepath, "csv_files", "GUT_predicted_mature_known_filtered_hypothetical.csv"), 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)

# HEMOLIMPH
## GOOD ONES
HEM_mature <- read_csv(file.path(basepath, "csv_files", "HEM_predicted_mature_known.csv"), show_col_types = FALSE)
HEM_mature['Sample'] = 'HEM'
HEM_mature <- HEM_mature %>%
  filter(`significant randfold p-value` == "yes") %>%
  filter(`miRDeep2 score` >= 1) %>%
  filter(`mature read count` >= 10) %>%
  filter(`star read count` >= 1)
write.table(
  HEM_mature, 
  file.path(basepath, "csv_files", "HEM_predicted_mature_known_filtered.csv"),
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)
## HYPOTHETICAL ONES
HEM_mature_hypothetical <- read_csv(file.path(basepath, "csv_files", "HEM_predicted_mature_known.csv"), show_col_types = FALSE)
HEM_mature_hypothetical['Sample'] = 'HEM'
HEM_mature_hypothetical <- HEM_mature_hypothetical %>%
  filter(`significant randfold p-value` != "yes") %>%
  filter(`miRDeep2 score` < 1) %>%
  filter(`mature read count` < 10) %>%
  filter(`star read count` < 1)
write.table(
  HEM_mature_hypothetical, 
  file.path(basepath, "csv_files", "HEM_predicted_mature_known_filtered_hypothetical.csv"), 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)

nrow(GLA_mature)
nrow(GLA_mature_hypothetical)
nrow(GUT_mature)
nrow(GUT_mature_hypothetical)
nrow(HEM_mature)
nrow(HEM_mature_hypothetical)
