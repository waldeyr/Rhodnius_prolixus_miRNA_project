library(dplyr)
library(ggplot2)
library(readxl)

basepath = dirname(rstudioapi::getSourceEditorContext()$path)
# Reading data

################################################################################
# Taxonomy
################################################################################
mirbase_organisms <- read_delim(file.path(basepath, "txt_files", "mirbase_organisms.txt"), 
                                "\t", 
                                escape_double = FALSE, 
                                trim_ws = TRUE, 
                                show_col_types = FALSE)
colnames(mirbase_organisms) <- c("organism",	"division",	"Species",	"tree",	"NCBI-taxid")

Known_miRNA <- read_excel(file.path(basepath, "csv_files", "Known_miRNA.xlsx"), sheet = "Known_miRNA")
Known_miRNA$organism <- substr(Known_miRNA$`mature miRBase miRNA`, 1, 3)
Known_miRNA <- left_join(
  Known_miRNA,
  mirbase_organisms,
  by = c("organism")
)
# View(Known_miRNA)

result_basepath = dirname(rstudioapi::getSourceEditorContext()$path)

plot <- Known_miRNA %>%
 filter(!(`predicted mature seq. in accordance with miRBase mature seq.` %in% "NA")) %>%
 ggplot() +
 aes(x = `total read count`, y = Species, fill = Tissue) +
 geom_boxplot(alpha=.25) +
 geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .2, alpha=0.6)+
 scale_fill_viridis_d(option = "viridis", direction = 1) +
 scale_x_continuous(trans = "log10", breaks = c(100, 1000, 10000, 100000, 1000000)) +
  xlab(expression(paste("Abundance of ", italic("Rhodnius prolixus"), " miRNA reads")))+
 theme_light() +
 theme(
   legend.position = "none",
   text = element_text(size = 20),
   axis.text.y = element_text(face = "italic"),
   axis.text.x = element_text(size=10, angle=60, hjust=1)
   ) +
 facet_wrap(vars(Tissue))
plot
ggsave(file = file.path(result_basepath, "Plots", "Figure_02_The abundance of Rhodnius prolixus miRNA reads aligned against known
miRNAs database.pdf"), plot, width = 360, height = 360, units = "mm")

save.image(file='Environment.RData')
