library("readr")
library("VennDiagram")
library("RColorBrewer")
library("dplyr")
library("grDevices")
library(readxl)

basepath = dirname(rstudioapi::getSourceEditorContext()$path)

Known_miRNA <- read_excel(file.path(basepath, "csv_files", "mirnas_venn.xlsx"), sheet = "known_venn")

# View(Known_miRNA[Known_miRNA$Tissue == 'Gut',])

x <- list(
  Known_miRNA[Known_miRNA$Tissue == 'Gut',]$`miRNA`,
  Known_miRNA[Known_miRNA$Tissue == 'Hemolymph',]$`miRNA`,
  Known_miRNA[Known_miRNA$Tissue == 'Salivary gland',]$`miRNA`
)
names <- c("Gut", "Hemolymph", "Salivary gland")
# View(x)

# # Chart
myCol <- brewer.pal(3, "Accent")
venn <- venn.diagram(
  x,
  category.names = names,
  filename = NULL,
  # Output features
  imagetype="png" ,
  height = 960 , 
  width = 960 , 
  resolution = 300,
  compression = "lzw",
  #title
  #main ="Known Mature miRNA",
  #sub = "Unique and shared novel mature miRNA among the three tissues",
  main.cex = 2,
  sub.cex = .9,
  main.fontfamily = "sans",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .9,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = .9,
  cat.fontface = "bold",
  cat.default.pos = "text",
  # euler.d = 3,
  cat.pos = c(-27, 27, 167),
  cat.dist = c(.06, .06, -.04),
  cat.fontfamily = "sans",
  rotation = 1
)

pdf(file = file.path(basepath, "Plots", "Figure_03_Distribution profile of Rhodnius prolixus miRnomes.pdf"))
grid.draw(venn)
dev.off()

#Remove log files
#system( paste("rm ", file.path(basepath, "Plots", "*.log")) )