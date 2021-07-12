library("readr")
library("VennDiagram")
library("RColorBrewer")

# GUT
RP1G_predicted_mature_novel_filtered <- read_delim(file.path(getwd(), "files_pipeline", "RP1G_predicted_mature_novel_filtered.tab"), 
                                                   "\t", escape_double = FALSE, trim_ws = TRUE)

# HEMOLIMPH
RP2H_mature_novel_filtered <- read_delim(file.path(getwd(), "files_pipeline", "RP2H_predicted_mature_novel_filtered.tab"), 
                                         "\t", escape_double = FALSE, trim_ws = TRUE)


# SALIVARY GLAND
RP_Gland_mature_novel_filtered <- read_delim(file.path(getwd(), "files_pipeline", "RP-Gland_predicted_mature_novel_filtered.tab"), 
                                             "\t", escape_double = FALSE, trim_ws = TRUE)

nrow(RP1G_predicted_mature_novel_filtered)
nrow(RP2H_mature_novel_filtered)
nrow(RP_Gland_mature_novel_filtered)

# Chart
myCol <- brewer.pal(3, "Accent")
venn <- venn.diagram(
  x = list(RP1G_predicted_mature_novel_filtered$`consensus mature sequence`, RP2H_mature_novel_filtered$`consensus mature sequence`, RP_Gland_mature_novel_filtered$`consensus mature sequence`),
  category.names = c("Gut" , "Hemolymph" , "Salivary gland"),
  filename = file.path(getwd(), "files_pipeline", "Venn_predicted_filtered_miRNA_novel.png"),
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 960 , 
  width = 960 , 
  resolution = 300,
  compression = "lzw",
  
  #title
  main ="Novel Mature miRNA",
  sub = "Unique and shared novel mature miRNA among the three tissues",
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
  cat.cex = .6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

#Remove log files
system( paste("rm ", file.path(getwd(), "files_pipeline", "*.log")) )