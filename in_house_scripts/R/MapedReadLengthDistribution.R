library(readr)
library(ggplot2)
library(gridExtra)
RP1G <- read_delim("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/RP1G.arf", 
                   "\t", escape_double = FALSE, col_names = FALSE, 
                   trim_ws = TRUE)

RP1G['Tissue'] = "Gut"
RP2H <- read_delim("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/RP2H.arf", 
                   "\t", escape_double = FALSE, col_names = FALSE, 
                   trim_ws = TRUE)
RP2H['Tissue'] =  "Hemolymph"
RPGland <- read_delim("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/RPGland.arf", 
                   "\t", escape_double = FALSE, col_names = FALSE, 
                   trim_ws = TRUE)
RPGland['Tissue'] =  "Salivary gland"

df = rbind(RP1G[, c("Tissue","X2")], RP2H[, c("Tissue","X2")], RPGland[, c("Tissue","X2")])
options(scipen=10000)


plot_mapped <- ggplot(df) +
  aes(x = X2, fill = Tissue) +
  #geom_histogram(bins = 30L) +
  facet_wrap(~Tissue, scales = 'free_x') +
  geom_histogram(bins=45) +
  scale_fill_viridis_d(option = "viridis") +
  labs(x = "Length", y = "Number of reads", title = "Length distribution of the mapped reads", caption = "Tissue", fill = "Tissue") +
  scale_x_continuous(breaks = seq(15, 50, by = 5)) +
  theme_light() +
  theme(legend.position = "right") +
  ylim(0L, 400000L)

#plot



all_trimmed_reads_lenght_distribution <- read_delim("Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/all_trimmed_reads_lenght_distribution.csv", 
                                                    "\t", escape_double = FALSE, trim_ws = TRUE)
plot_trimmed <- ggplot(all_trimmed_reads_lenght_distribution) +
 facet_wrap(~Tissue, scales = "free_x") +  
 aes(x = Length, y = Count, fill = Tissue) +
 geom_area(size = 1L) +
 scale_fill_viridis_d(option = "viridis") +
 labs(y = "Number of reads", title = "Length distribution of the trimmed raw reads") +
  scale_x_continuous(breaks = seq(15, 50, by = 5)) +
  theme_light()
#plot2


grid.arrange(plot_trimmed, plot_mapped, nrow = 2)
