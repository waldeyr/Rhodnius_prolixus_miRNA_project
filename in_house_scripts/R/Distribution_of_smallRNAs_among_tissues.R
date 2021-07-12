library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(plotly)

#RP1G = Gut
RP1G_mapped_blast_results <- read_delim("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/RP1G_mapped_blast_results.tab", 
                                        "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(RP1G_mapped_blast_results) <- c("read", "rna_central", "identity", "alignment_length", "mismatches", "gap_opens", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score")
RP1G_mapped_blast_results <- RP1G_mapped_blast_results %>% distinct(read, .keep_all = TRUE)
#View(RP1G_mapped_blast_results)

#RP2H = Hemolymph
RP2H_mapped_blast_results <- read_delim("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/RP2H_mapped_blast_results.tab", 
                                        "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(RP2H_mapped_blast_results) <- c("read", "rna_central", "identity", "alignment_length", "mismatches", "gap_opens", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score")
RP2H_mapped_blast_results <- RP2H_mapped_blast_results %>% distinct(read, .keep_all = TRUE)
#head(RP2H_mapped_blast_results)

#RPGland = Salivary gland
RPGland_mapped_blast_results <- read_delim("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/RPGland_mapped_blast_results.tab", 
                                        "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(RPGland_mapped_blast_results) <- c("read", "rna_central", "identity", "alignment_length", "mismatches", "gap_opens", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score")
RPGland_mapped_blast_results <- RPGland_mapped_blast_results %>% distinct(read, .keep_all = TRUE)
#head(RPGland_mapped_blast_results)

#Taxonomy
rnacentral_id_taxonomy_annotation <- read_delim("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/RnaCentral/rnacentral_id_taxonomy_annotation.txt", 
                                                "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(rnacentral_id_taxonomy_annotation) <- c("rna_central", "taxonomy", "annotation")
#View(rnacentral_id_taxonomy_annotation)

RP1G_mapped_blast_results_rnacentral <- merge(RP1G_mapped_blast_results, rnacentral_id_taxonomy_annotation, by=c("rna_central","rna_central")) 
temp <- str_split_fixed(as.character(RP1G_mapped_blast_results_rnacentral$taxonomy)," ", 2) 
RP1G_mapped_blast_results_rnacentral['genus'] <- temp[,1]
RP1G_mapped_blast_results_rnacentral$rna_type <- ""


# GUT 
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)ribosomal", RP1G_mapped_blast_results_rnacentral$annotation), "rRNA", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)rRNA", RP1G_mapped_blast_results_rnacentral$annotation), "rRNA", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)tRNA", RP1G_mapped_blast_results_rnacentral$annotation), "tRNA", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)spliceosomal", RP1G_mapped_blast_results_rnacentral$annotation), "spliceosomal RNA", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)small nucleolar", RP1G_mapped_blast_results_rnacentral$annotation), "snoRNA", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)signal recognition", RP1G_mapped_blast_results_rnacentral$annotation), "srpRNA", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)antisense", RP1G_mapped_blast_results_rnacentral$annotation), "asRNA", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)Nuclear RNase P", RP1G_mapped_blast_results_rnacentral$annotation), "Nuclear RNase P", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)transfer", RP1G_mapped_blast_results_rnacentral$annotation), "tRNA", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)rNase MRP", RP1G_mapped_blast_results_rnacentral$annotation), "RNase MRP", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)sRNA", RP1G_mapped_blast_results_rnacentral$annotation), "Bacterial small RNA", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)RybB", RP1G_mapped_blast_results_rnacentral$annotation), "RybB RNA", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)Glm", RP1G_mapped_blast_results_rnacentral$annotation), "GlmY RNA", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)CsrB", RP1G_mapped_blast_results_rnacentral$annotation), "CsrB MRP", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)Hfq", RP1G_mapped_blast_results_rnacentral$annotation), "Hfq binding sRNA", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)TPP riboswitch", RP1G_mapped_blast_results_rnacentral$annotation), "TPP riboswitch", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)7sk", RP1G_mapped_blast_results_rnacentral$annotation), "7SK RNA", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)nhaA-I", RP1G_mapped_blast_results_rnacentral$annotation), "nhaA-I RNA", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)thermometer", RP1G_mapped_blast_results_rnacentral$annotation), "RNA thermometer", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)KDPG-aldolase", RP1G_mapped_blast_results_rnacentral$annotation), "KDPG-aldolase RNA motif", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)uncharacterized", RP1G_mapped_blast_results_rnacentral$annotation), "uncharacterized hypothetical RNA", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)GcvB", RP1G_mapped_blast_results_rnacentral$annotation), "GcvB RNA", RP1G_mapped_blast_results_rnacentral$rna_type)
RP1G_mapped_blast_results_rnacentral$rna_type <- ifelse( is.na(RP1G_mapped_blast_results_rnacentral$annotation), "uncharacterized hypothetical RNA", RP1G_mapped_blast_results_rnacentral$rna_type)

################################################################################
# Join Gut miRNA filtered results
RP1G_mature_novel_filtered <- read_delim("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP1G_mature_novel_filtered.tab", 
                                         "\t", escape_double = FALSE, comment = "#", 
                                         trim_ws = TRUE)
RP1G_mature_novel_filtered$rna_type <- "probable novel miRNA"
RP1G_mature_novel_filtered$class <- "miRNA"
RP1G_mature_novel_filtered_to_plot <- RP1G_mature_novel_filtered %>%
  select("rna_type", "class")
################################################################################
RP1G_mature_known_filtered <- read_delim("~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP1G_mature_known_filtered.tab", 
                                         "\t", escape_double = FALSE, comment = "#", 
                                         trim_ws = TRUE)
RP1G_mature_known_filtered$rna_type <- "probable known miRNA"
RP1G_mature_known_filtered$class <- "miRNA"
RP1G_mature_known_filtered_to_plot <- RP1G_mature_known_filtered %>%
  select("rna_type","class")
################################################################################
RP1G_mapped_blast_results_rnacentral['class'] <- "other small RNA"
RP1G_mapped_blast_results_rnacentral_to_plot <- RP1G_mapped_blast_results_rnacentral %>%
  select("rna_type", "class")
RP1G_to_plot <- bind_rows( RP1G_mapped_blast_results_rnacentral_to_plot,
                      bind_rows( 
                        RP1G_mature_novel_filtered_to_plot, 
                        RP1G_mature_known_filtered_to_plot)
  )

# HEMOLIMPH
RP2H_mapped_blast_results_rnacentral <- merge(RP2H_mapped_blast_results, rnacentral_id_taxonomy_annotation, by=c("rna_central","rna_central")) 
temp <- str_split_fixed(as.character(RP2H_mapped_blast_results_rnacentral$taxonomy)," ", 2) 
RP2H_mapped_blast_results_rnacentral['genus'] <- temp[,1]
RP2H_mapped_blast_results_rnacentral$rna_type <- ""
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)ribosomal", RP2H_mapped_blast_results_rnacentral$annotation), "rRNA", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)rRNA", RP2H_mapped_blast_results_rnacentral$annotation), "rRNA", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)tRNA", RP2H_mapped_blast_results_rnacentral$annotation), "tRNA", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)spliceosomal", RP2H_mapped_blast_results_rnacentral$annotation), "spliceosomal RNA", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)small nucleolar", RP2H_mapped_blast_results_rnacentral$annotation), "snoRNA", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)signal recognition", RP2H_mapped_blast_results_rnacentral$annotation), "srpRNA", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)antisense", RP2H_mapped_blast_results_rnacentral$annotation), "asRNA", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)Nuclear RNase P", RP2H_mapped_blast_results_rnacentral$annotation), "Nuclear RNase P", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)transfer", RP2H_mapped_blast_results_rnacentral$annotation), "tRNA", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)rNase MRP", RP2H_mapped_blast_results_rnacentral$annotation), "RNase MRP", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)sRNA", RP2H_mapped_blast_results_rnacentral$annotation), "Bacterial small RNA", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)RybB", RP2H_mapped_blast_results_rnacentral$annotation), "RybB RNA", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)Glm", RP2H_mapped_blast_results_rnacentral$annotation), "GlmY RNA", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)CsrB", RP2H_mapped_blast_results_rnacentral$annotation), "CsrB MRP", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)Hfq", RP2H_mapped_blast_results_rnacentral$annotation), "Hfq binding sRNA", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)TPP riboswitch", RP2H_mapped_blast_results_rnacentral$annotation), "TPP riboswitch", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)7sk", RP2H_mapped_blast_results_rnacentral$annotation), "7SK RNA", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)nhaA-I", RP2H_mapped_blast_results_rnacentral$annotation), "nhaA-I RNA", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)thermometer", RP2H_mapped_blast_results_rnacentral$annotation), "RNA thermometer", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)KDPG-aldolase", RP2H_mapped_blast_results_rnacentral$annotation), "KDPG-aldolase RNA motif", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)uncharacterized", RP2H_mapped_blast_results_rnacentral$annotation), "uncharacterized hypothetical RNA", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)GcvB", RP2H_mapped_blast_results_rnacentral$annotation), "GcvB RNA", RP2H_mapped_blast_results_rnacentral$rna_type)
RP2H_mapped_blast_results_rnacentral$rna_type <- ifelse( is.na(RP2H_mapped_blast_results_rnacentral$annotation), "uncharacterized hypothetical RNA", RP2H_mapped_blast_results_rnacentral$rna_type)

# Join Hemolimph miRNA filtered results
RP2H_mature_novel_filtered <- read_delim("Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP2H_mature_novel_filtered.tab", 
                                         "\t", escape_double = FALSE, comment = "#", 
                                         trim_ws = TRUE)
RP2H_mature_novel_filtered$rna_type <- "probable novel miRNA"

RP2H_mature_known_filtered <- read_delim("Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RP2H_mature_known_filtered.tab", 
                                         "\t", escape_double = FALSE, comment = "#", 
                                         trim_ws = TRUE)
RP2H_mature_known_filtered$rna_type <- "probable known miRNA"

RP2H_mapped_blast_results_rnacentral_to_plot <- RP2H_mapped_blast_results_rnacentral %>%
  select("rna_type")
RP2H_mature_known_filtered_to_plot <- RP2H_mature_known_filtered %>%
  select("rna_type")
RP2H_mature_novel_filtered_to_plot <- RP2H_mature_novel_filtered %>%
  select("rna_type")
RP2H_to_plot <- bind_rows(bind_rows(RP2H_mapped_blast_results_rnacentral_to_plot, RP2H_mature_novel_filtered_to_plot), RP2H_mature_known_filtered_to_plot)

# SALIVARY GLAND
RPGland_mapped_blast_results_rnacentral <- merge(RPGland_mapped_blast_results, rnacentral_id_taxonomy_annotation, by=c("rna_central","rna_central")) 
temp <- str_split_fixed(as.character(RPGland_mapped_blast_results_rnacentral$taxonomy)," ", 2) 
RPGland_mapped_blast_results_rnacentral['genus'] <- temp[,1]
RPGland_mapped_blast_results_rnacentral$rna_type <- ""
RPGland_mapped_blast_results_rnacentral$rna_type <- ""
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)ribosomal", RPGland_mapped_blast_results_rnacentral$annotation), "rRNA", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)rRNA", RPGland_mapped_blast_results_rnacentral$annotation), "rRNA", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)tRNA", RPGland_mapped_blast_results_rnacentral$annotation), "tRNA", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)spliceosomal", RPGland_mapped_blast_results_rnacentral$annotation), "spliceosomal RNA", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)small nucleolar", RPGland_mapped_blast_results_rnacentral$annotation), "snoRNA", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)signal recognition", RPGland_mapped_blast_results_rnacentral$annotation), "srpRNA", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)antisense", RPGland_mapped_blast_results_rnacentral$annotation), "asRNA", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)Nuclear RNase P", RPGland_mapped_blast_results_rnacentral$annotation), "Nuclear RNase P", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)transfer", RPGland_mapped_blast_results_rnacentral$annotation), "tRNA", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)rNase MRP", RPGland_mapped_blast_results_rnacentral$annotation), "RNase MRP", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)sRNA", RPGland_mapped_blast_results_rnacentral$annotation), "Bacterial small RNA", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)RybB", RPGland_mapped_blast_results_rnacentral$annotation), "RybB RNA", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)Glm", RPGland_mapped_blast_results_rnacentral$annotation), "GlmY RNA", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)CsrB", RPGland_mapped_blast_results_rnacentral$annotation), "CsrB MRP", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)Hfq", RPGland_mapped_blast_results_rnacentral$annotation), "Hfq binding sRNA", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)TPP riboswitch", RPGland_mapped_blast_results_rnacentral$annotation), "TPP riboswitch", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)7sk", RPGland_mapped_blast_results_rnacentral$annotation), "7SK RNA", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)nhaA-I", RPGland_mapped_blast_results_rnacentral$annotation), "nhaA-I RNA", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)thermometer", RPGland_mapped_blast_results_rnacentral$annotation), "RNA thermometer", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)KDPG-aldolase", RPGland_mapped_blast_results_rnacentral$annotation), "KDPG-aldolase RNA motif", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)uncharacterized", RPGland_mapped_blast_results_rnacentral$annotation), "uncharacterized hypothetical RNA", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( grepl("(?i)GcvB", RPGland_mapped_blast_results_rnacentral$annotation), "GcvB RNA", RPGland_mapped_blast_results_rnacentral$rna_type)
RPGland_mapped_blast_results_rnacentral$rna_type <- ifelse( is.na(RPGland_mapped_blast_results_rnacentral$annotation), "uncharacterized hypothetical RNA", RPGland_mapped_blast_results_rnacentral$rna_type)

# Join Salivary gland miRNA filtered results
RPGland_mature_novel_filtered <- read_delim("Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RPGland_mature_novel_filtered.tab", 
                                         "\t", escape_double = FALSE, comment = "#", 
                                         trim_ws = TRUE)
RPGland_mature_novel_filtered$rna_type <- "probable novel miRNA"

RPGland_mature_known_filtered <- read_delim("Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/VennDriagram/RPGland_mature_known_filtered.tab", 
                                         "\t", escape_double = FALSE, comment = "#", 
                                         trim_ws = TRUE)
RPGland_mature_known_filtered$rna_type <- "probable known miRNA"

RPGland_mapped_blast_results_rnacentral_to_plot <- RPGland_mapped_blast_results_rnacentral %>%
  select("rna_type")
RPGland_mature_known_filtered_to_plot <- RPGland_mature_known_filtered %>%
  select("rna_type")
RPGland_mature_novel_filtered_to_plot <- RPGland_mature_novel_filtered %>%
  select("rna_type")
RPGland_to_plot <- bind_rows(bind_rows(RPGland_mapped_blast_results_rnacentral_to_plot, RPGland_mature_novel_filtered_to_plot), RPGland_mature_known_filtered_to_plot)

rm(temp)

Gut_pie_data <- data.frame(RP1G_to_plot$rna_type)
colnames(Gut_pie_data) <- c("rna_type")
Gut_pie_data <- Gut_pie_data %>%
  group_by(rna_type) %>%
  count() %>%
  ungroup() %>%
  mutate(per=`n`/sum(`n`))
Gut_pie_data$label <- scales::percent(Gut_pie_data$per)
Gut_pie_data$tissue <- "Gut"
#Gut_pie_data

Hemolymph_pie_data <- data.frame(RP2H_to_plot$rna_type)
colnames(Hemolymph_pie_data) <- c("rna_type")
Hemolymph_pie_data <- Hemolymph_pie_data %>%
  group_by(rna_type) %>%
  count() %>%
  ungroup() %>%
  mutate(per=`n`/sum(`n`))
Hemolymph_pie_data$label <- scales::percent(Hemolymph_pie_data$per)
Hemolymph_pie_data$tissue <- "Hemolymph"
#Hemolymph_pie_data

Salivary_gland_pie_data <- data.frame(RP1G_to_plot$rna_type)
colnames(Salivary_gland_pie_data) <- c("rna_type")
Salivary_gland_pie_data <- Salivary_gland_pie_data %>%
  group_by(rna_type) %>%
  count() %>%
  ungroup() %>%
  mutate(per=`n`/sum(`n`))
  Salivary_gland_pie_data$label <- scales::percent(Salivary_gland_pie_data$per)
Salivary_gland_pie_data$tissue <- "Salivary gland"
#Salivary_gland_pie_data

All_pie <- rbind(rbind(Gut_pie_data, Hemolymph_pie_data),Salivary_gland_pie_data)

View(All_pie)

## set the levels in order we want
All_pie <- All_pie[order(factor(All_pie$tissue), factor(All_pie$per)),]


write.table(
  All_pie, 
  "~/Insync/GoogleDrive/projeto_paula/pipeline_mirDeep2/pipeline/inhouse_scripts/All_pie.tab", 
  append = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)


plot <- ggplot(All_pie, aes(x = rna_type, fill = rna_type, group = rna_type, weight = per)) +
 geom_bar() +
 scale_y_continuous(limits = c(0,1)) +
 scale_fill_viridis_d(option = "viridis") +
 labs(x = "small RNA", y = "Percentage among the tissues") +
 coord_flip() +
 geom_text(aes(y=per, label=label), position=position_dodge(width=0.9), vjust=0.5, hjust=-.2, size=3) +
 theme(legend.position = "none") +
 facet_wrap(vars(tissue), scales = "free", ncol=1, strip.position = "left")

plot