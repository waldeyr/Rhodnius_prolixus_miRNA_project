# Este script (maduros conhecidos):
# 1) lê os resultados do mirDeep2 para cada amostra
# 2) filtra os resultados segundo alguns critérios usados em trabalhos anteriores
# 3) cria um diagrama de Venn
# 4) salva os resultados após filtro

library("readr")
library("dplyr")
library("VennDiagram")
library("RColorBrewer")


# SALIVARY GLAND
## GOOD ONES
RPGland_mature <- read_csv("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP-Gland_mature.csv")
RPGland_mature['Sample'] = 'RPGland'
RPGland_mature <- RPGland_mature %>%
  filter(`significant randfold p-value` == "yes") %>%
  filter(`miRDeep2 score` >= 1) %>%
  filter(`mature read count` >= 10) %>%
  filter(`star read count` > 0)
write.table(
  RPGland_mature, 
  "~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RPGland_mature_known_filtered.tab", 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)

## HYPOTHETICAL ONES
RPGland_mature_hypothetical <- read_csv("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP-Gland_mature.csv")
RPGland_mature_hypothetical['Sample'] = 'RPGland'
RPGland_mature_hypothetical <- RPGland_mature_hypothetical %>%
  filter(`significant randfold p-value` != "yes") %>%
  filter(`miRDeep2 score` < 1) %>%
  filter(`mature read count` < 10) %>%
  filter(`star read count` <= 0)
write.table(
  RPGland_mature_hypothetical, 
  "~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RPGland_mature_known_hypothetical.tab", 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)

# GUT
## GOOD ONES
RP1G_mature <- read_csv("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP1G_mature.csv")
RP1G_mature['Sample'] = 'RP1G'
RP1G_mature <- RP1G_mature %>%
  filter(`significant randfold p-value` == "yes") %>%
  filter(`miRDeep2 score` >= 1) %>%
  filter(`mature read count` >= 10) %>%
  filter(`star read count` > 0)
write.table(
  RP1G_mature, 
  "~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP1G_mature_known_filtered.tab", 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)
## HYPOTHETICAL ONES
RP1G_mature_hypothetical <- read_csv("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP1G_mature.csv")
RP1G_mature_hypothetical['Sample'] = 'RP1G'
RP1G_mature_hypothetical <- RP1G_mature_hypothetical %>%
  filter(`significant randfold p-value` != "yes") %>%
  filter(`miRDeep2 score` < 1) %>%
  filter(`mature read count` < 10) %>%
  filter(`star read count` <= 0)
write.table(
  RP1G_mature_hypothetical, 
  "~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP1G_mature_known_filtered_hypothetical.tab", 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)

# HEMOLIMPH
## GOOD ONES
RP2H_mature <- read_csv("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP2H_mature.csv")
RP2H_mature['Sample'] = 'RP2H'
RP2H_mature <- RP2H_mature %>%
  filter(`significant randfold p-value` == "yes") %>%
  filter(`miRDeep2 score` >= 1) %>%
  filter(`mature read count` >=10) %>%
  filter(`star read count` > 0)
write.table(
  RP2H_mature, 
  "~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP2H_mature_known_filtered.tab", 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)
## HYPOTHETICAL ONES
RP2H_mature_hypothetical <- read_csv("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP2H_mature.csv")
RP2H_mature_hypothetical['Sample'] = 'RP2H'
RP2H_mature_hypothetical <- RP2H_mature_hypothetical %>%
  filter(`significant randfold p-value` != "yes") %>%
  filter(`miRDeep2 score` < 1) %>%
  filter(`mature read count` < 10) %>%
  filter(`star read count` <= 0)
write.table(
  RP2H_mature_hypothetical, 
  "~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP2H_mature_known_filtered_hypothetical.tab", 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)

nrow(RPGland_mature)
nrow(RP1G_mature)
nrow(RP2H_mature)

# Chart

myCol <- brewer.pal(3, "Pastel2")
venn <- venn.diagram(
  x = list(RP1G_mature$`consensus mature sequence`, RP2H_mature$`consensus mature sequence`, RPGland_mature$`consensus mature sequence`),
  category.names = c("Gut" , "Hemolymph" , "Salivary gland"),
  filename = '~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/miRNA.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 960 , 
  width = 960 , 
  resolution = 300,
  compression = "lzw",
  
  #title
  main ="Mature miRNA annotated with miRBase",
  sub = "Unique and shared miRNA among the three tissues",
  main.cex = .9,
  sub.cex = .5,
  main.fontfamily = "sans",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)
