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
RP1G_mature <- read_csv("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP1G_mature_new.csv")
RP1G_mature['Sample'] = 'RP1G'
RP1G_mature <- RP1G_mature %>%
  filter(`significant randfold p-value` == "yes") %>%
  filter(`miRDeep2 score` >= 1) %>%
  filter(`mature read count` >=10) %>%
  filter(`star read count` > 0)
write.table(
  RP1G_mature, 
  "~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP1G_mature_novel_filtered.tab", 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)
## Hypothetical ones
RP1G_mature_hypothetical <- read_csv("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP1G_mature_new.csv")
RP1G_mature_hypothetical['Sample'] = 'RP1G'
RP1G_mature_hypothetical <- RP1G_mature_hypothetical %>%
  filter(`significant randfold p-value` != "yes") %>%
  filter(`miRDeep2 score` < 1) %>%
  filter(`mature read count` < 10) %>%
  filter(`star read count` <= 0)
write.table(
  RP1G_mature_hypothetical, 
  "~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP1G_mature_novel_filtered_hypothetical.tab", 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)




# HEMOLIMPH
## Good ones
RP2H_mature <- read_csv("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP2H_mature_new.csv")
RP2H_mature['Sample'] = 'RP2H'
RP2H_mature <- RP2H_mature %>%
  filter(`significant randfold p-value` == "yes") %>%
  filter(`miRDeep2 score` >= 1) %>%
  filter(`mature read count` >=10) %>%
  filter(`star read count` > 0)
write.table(
  RP2H_mature, 
  "~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP2H_mature_novel_filtered.tab", 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)
## Hypothetical ones
RP2H_mature_hypothetical <- read_csv("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP2H_mature_new.csv")
RP2H_mature_hypothetical['Sample'] = 'RP2H'
RP2H_mature_hypothetical <- RP2H_mature_hypothetical %>%
  filter(`significant randfold p-value` != "yes") %>%
  filter(`miRDeep2 score` < 1) %>%
  filter(`mature read count` < 10) %>%
  filter(`star read count` <= 0)
write.table(
  RP2H_mature_hypothetical, 
  "~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP2H_mature_novel_filtered_hypothetical.tab", 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)




# SALIVARY GLAND
## Good ones
RPGland_mature <- read_csv("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP-Gland_mature_new.csv")
RPGland_mature['Sample'] = 'RPGland'
RPGland_mature <- RPGland_mature %>%
  filter(`significant randfold p-value` == "yes") %>%
  filter(`miRDeep2 score` >= 1) %>%
  filter(`mature read count` >=10) %>%
  filter(`star read count` > 0)
write.table(
  RPGland_mature, 
  "~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RPGland_mature_novel_filtered.tab", 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)
## Hypothetical ones
RPGland_mature_hypothetical <- read_csv("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP-Gland_mature_new.csv")
RPGland_mature_hypothetical['Sample'] = 'RP2H'
RPGland_mature_hypothetical <- RPGland_mature_hypothetical %>%
  filter(`significant randfold p-value` != "yes") %>%
  filter(`miRDeep2 score` < 1) %>%
  filter(`mature read count` < 10) %>%
  filter(`star read count` <= 0)
write.table(
  RPGland_mature_hypothetical, 
  "~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RPGland_mature_novel_filtered_hypothetical.tab", 
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
  filename = '~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram//miRNA_novel.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 960 , 
  width = 960 , 
  resolution = 300,
  compression = "lzw",
  
  #title
  main ="Novel Mature miRNA",
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


## UNUSED DRAFTS

# df = rbind(RP1G_mature, RP2H_mature, RPGland_mature)
# 
# df2 <- within(df, { count_presence_of_mature_on_samples <- ave(Sample, `mature miRBase miRNA`, FUN=function(x) length(unique(x)))})
# #View(df2)
# 
# df3 <- df %>%
#   filter(`significant randfold p-value` == "yes") %>%
#   filter(`miRDeep2 score` >= 1) %>%
#   #filter(substr(`estimated probability that the miRNA candidate is a true positive`, 1, 2) > 51) %>%
#   filter(`mature read count`>=10)
# View(df3)