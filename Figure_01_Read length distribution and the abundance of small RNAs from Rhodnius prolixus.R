# Packages
library(readr)
library(ggplot2)
library(gridExtra)
# Takes this script dir as absolute path
basepath = dirname(rstudioapi::getSourceEditorContext()$path)
# Reading data
GUT <- read_delim(file.path(basepath, "arf_files", "GUT.arf"), 
                   "\t", 
                   escape_double = FALSE, 
                   col_names = FALSE, 
                   trim_ws = TRUE)
GUT['Tissue'] = "Gut"
###################
HEM <- read_delim(file.path(basepath, "arf_files", "HEM.arf"),  
                   "\t", 
                   escape_double = FALSE, 
                   col_names = FALSE, 
                   trim_ws = TRUE)
HEM['Tissue'] =  "Hemolymph"
###################
GLA <- read_delim(file.path(basepath, "arf_files", "GLA.arf"), 
                   "\t", 
                   escape_double = FALSE, 
                   col_names = FALSE, 
                   trim_ws = TRUE)
GLA['Tissue'] =  "Salivary Gland"
###################
df = rbind(GUT[, c("Tissue","X2")], HEM[, c("Tissue","X2")], GLA[, c("Tissue","X2")])
###################
All_trimmed_reads_lenght_distribution <- read_delim(file.path(basepath, "csv_files", "All_trimmed_reads_lenght_distribution.csv"),
                                                    "\t", 
                                                    escape_double = FALSE, 
                                                    trim_ws = TRUE)

# Plot

# numeric scale
options(scipen=10000)

plot_mapped <- ggplot(df) +
  aes(x = X2, fill = Tissue) +
  facet_wrap(~Tissue, scales = 'free_x') +
  geom_histogram(bins=45) +
  theme(legend.position = "none",
        text = element_text(size = 20)
        )+
  scale_fill_viridis_d(option = "viridis") +
  labs(x = "Length", y = "Number of reads", title = "Length distribution of the mapped reads",  fill = "Tissue") +
  scale_x_continuous(breaks = seq(15, 50, by = 5)) +
  theme_light() +
  ylim(0L, 400000L)

plot_trimmed <- ggplot(All_trimmed_reads_lenght_distribution) +
 facet_wrap(~Tissue, scales = "free_x") +  
 aes(x = Length, y = Count, fill = Tissue) +
 geom_area(size = 1L) +
theme(legend.position = "none",
      text = element_text(size = 20)
      )+
 scale_fill_viridis_d(option = "viridis") +
 labs(y = "Number of reads", title = "Length distribution of the trimmed raw reads") +
  scale_x_continuous(breaks = seq(15, 50, by = 5)) +
  theme_light()

grid.arrange(plot_trimmed, plot_mapped, nrow = 2)
g <- arrangeGrob(plot_trimmed, plot_mapped, nrow = 2)
ggsave(file = file.path(basepath, "Plots", "Figure_01_Read length distribution and the abundance of small RNAs from Rhodnius prolixus.pdf"), g, width = 300, height = 300, units = "mm")

