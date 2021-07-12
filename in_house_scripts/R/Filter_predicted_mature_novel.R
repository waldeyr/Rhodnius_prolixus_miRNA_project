# Este script (maduros novos):
# 1) lê os resultados do mirDeep2 para cada amostra
# 2) filtra os resultados segundo alguns critérios usados em trabalhos anteriores
# 3) cria um diagrama de Venn
# 4) salva os resultados após filtro

library("readr")
library("dplyr")
library("VennDiagram")
library("RColorBrewer")

# GUT
## Good ones
RP1G_mature <- read_csv(file.path(getwd(), "files_pipeline", "RP1G_predicted_mature_novel.csv"))
RP1G_mature['Sample'] = 'RP1G'
RP1G_mature <- RP1G_mature %>%
  filter(`significant randfold p-value` == "yes") %>%
  filter(`miRDeep2 score` >= 1) %>%
  filter(`mature read count` >=10) %>%
  filter(`star read count` > 0)
write.table(
  RP1G_mature, 
  file.path(getwd(), "files_pipeline", "RP1G_predicted_mature_novel_filtered.tab"), 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)

## Hypothetical ones
RP1G_mature_hypothetical <- read_csv(file.path(getwd(), "files_pipeline", "RP1G_predicted_mature_novel.csv"))
RP1G_mature_hypothetical['Sample'] = 'RP1G'
RP1G_mature_hypothetical <- RP1G_mature_hypothetical %>%
  filter(`significant randfold p-value` != "yes") %>%
  filter(`miRDeep2 score` < 1) %>%
  filter(`mature read count` < 10) %>%
  filter(`star read count` <= 0)
write.table(
  RP1G_mature_hypothetical, 
  file.path(getwd(), "files_pipeline", "RP1G_predicted_mature_novel_filtered_hypothetical.tab"), 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)


# HEMOLIMPH
## Good ones
RP2H_mature <- read_csv(file.path(getwd(), "files_pipeline", "RP2H_predicted_mature_novel.csv"))
RP2H_mature['Sample'] = 'RP2H'
RP2H_mature <- RP2H_mature %>%
  filter(`significant randfold p-value` == "yes") %>%
  filter(`miRDeep2 score` >= 1) %>%
  filter(`mature read count` >=10) %>%
  filter(`star read count` > 0)
write.table(
  RP2H_mature, 
  file.path(getwd(), "files_pipeline", "RP2H_predicted_mature_novel_filtered.tab"), 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)

## Hypothetical ones
RP2H_mature_hypothetical <- read_csv(file.path(getwd(), "files_pipeline", "RP2H_predicted_mature_novel.csv"))
RP2H_mature_hypothetical['Sample'] = 'RP2H'
RP2H_mature_hypothetical <- RP2H_mature_hypothetical %>%
  filter(`significant randfold p-value` != "yes") %>%
  filter(`miRDeep2 score` < 1) %>%
  filter(`mature read count` < 10) %>%
  filter(`star read count` <= 0)
write.table(
  RP2H_mature_hypothetical, 
  file.path(getwd(), "files_pipeline", "RP2H_predicted_mature_novel_filtered_hypothetical.tab"), 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)


# SALIVARY GLAND
## Good ones
RPGland_mature <- read_csv(file.path(getwd(), "files_pipeline", "RP-Gland_predicted_mature_novel.csv"))
RPGland_mature['Sample'] = 'RPGland'
RPGland_mature <- RPGland_mature %>%
  filter(`significant randfold p-value` == "yes") %>%
  filter(`miRDeep2 score` >= 1) %>%
  filter(`mature read count` >=10) %>%
  filter(`star read count` > 0)
write.table(
  RPGland_mature, 
  file.path(getwd(), "files_pipeline", "RP-Gland_predicted_mature_novel_filtered.tab"), 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)

## Hypothetical ones
RPGland_mature_hypothetical <- read_csv(file.path(getwd(), "files_pipeline", "RP-Gland_predicted_mature_novel.csv"))
RPGland_mature_hypothetical['Sample'] = 'RP2H'
RPGland_mature_hypothetical <- RPGland_mature_hypothetical %>%
  filter(`significant randfold p-value` != "yes") %>%
  filter(`miRDeep2 score` < 1) %>%
  filter(`mature read count` < 10) %>%
  filter(`star read count` <= 0)
write.table(
  RPGland_mature_hypothetical, 
  file.path(getwd(), "files_pipeline", "RP-Gland_predicted_mature_novel_filtered_hypothetical.tab"), 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)

nrow(RPGland_mature)
nrow(RP1G_mature)
nrow(RP2H_mature)