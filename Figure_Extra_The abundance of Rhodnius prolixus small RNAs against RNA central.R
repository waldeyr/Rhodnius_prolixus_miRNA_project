# Packages
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(plotly)
require(scales) # to access break formatting functions
# Takes this script dir as absolute path
basepath = dirname(rstudioapi::getSourceEditorContext()$path)
# Reading data


# Taxonomy
rnacentral_id_taxonomy_annotation <- read_delim(file.path(basepath, "txt_files", "taxonomy_annotation.txt"), 
                                                "\t", 
                                                escape_double = FALSE, 
                                                trim_ws = TRUE, 
                                                show_col_types = FALSE)
colnames(rnacentral_id_taxonomy_annotation) <- c("rna_central", "taxonomy", "annotation")
# View(rnacentral_id_taxonomy_annotation)

# GUT 
GUT_small_rna <- read_delim(file.path(basepath, "csv_files", "GUT_other_small_rna_blast_results.tab"), 
                            "\t", 
                            col_names = FALSE, 
                            escape_double = FALSE, 
                            trim_ws = TRUE,
                            show_col_types = FALSE)
colnames(GUT_small_rna) <- c("read", "rna_central", "identity", "alignment_length", "mismatches", "gap_opens", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score")
GUT_small_rna <- GUT_small_rna %>% distinct(read, .keep_all = TRUE)
GUT_small_rna <- GUT_small_rna %>%
  separate(read, c("read", "total_read_count"), "_x") # Insert read counts by spliting read column
GUT_small_rna$"total_read_count" <- sapply(GUT_small_rna$"total_read_count", as.numeric)
GUT_small_rna_rnacentral <- merge(GUT_small_rna, rnacentral_id_taxonomy_annotation, by=c("rna_central","rna_central")) 
temp <- str_split_fixed(as.character(GUT_small_rna_rnacentral$taxonomy)," ", 2) 
GUT_small_rna_rnacentral['genus'] <- temp[,1]
GUT_small_rna_rnacentral$rna_type <- ""
GUT_small_rna_rnacentral$rna_type <- ""
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)ribosomal", GUT_small_rna_rnacentral$annotation), "rRNA", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)rRNA", GUT_small_rna_rnacentral$annotation), "rRNA", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)tRNA", GUT_small_rna_rnacentral$annotation), "tRNA", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)spliceosomal", GUT_small_rna_rnacentral$annotation), "spliceosomal RNA", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)small nucleolar", GUT_small_rna_rnacentral$annotation), "snoRNA", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)signal recognition", GUT_small_rna_rnacentral$annotation), "srpRNA", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)antisense", GUT_small_rna_rnacentral$annotation), "asRNA", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)Nuclear RNase P", GUT_small_rna_rnacentral$annotation), "Nuclear RNase P", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)transfer", GUT_small_rna_rnacentral$annotation), "tRNA", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)rNase MRP", GUT_small_rna_rnacentral$annotation), "RNase MRP", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)sRNA", GUT_small_rna_rnacentral$annotation), "Bacterial small RNA", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)RybB", GUT_small_rna_rnacentral$annotation), "RybB RNA", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)Glm", GUT_small_rna_rnacentral$annotation), "GlmY RNA", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)CsrB", GUT_small_rna_rnacentral$annotation), "CsrB MRP", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)Hfq", GUT_small_rna_rnacentral$annotation), "Hfq binding sRNA", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)TPP riboswitch", GUT_small_rna_rnacentral$annotation), "TPP riboswitch", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)7sk", GUT_small_rna_rnacentral$annotation), "7SK RNA", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)nhaA-I", GUT_small_rna_rnacentral$annotation), "nhaA-I RNA", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)thermometer", GUT_small_rna_rnacentral$annotation), "RNA thermometer", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)KDPG-aldolase", GUT_small_rna_rnacentral$annotation), "KDPG-aldolase RNA motif", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)uncharacterized", GUT_small_rna_rnacentral$annotation), "uncharacterized hypothetical RNA", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)GcvB", GUT_small_rna_rnacentral$annotation), "GcvB RNA", GUT_small_rna_rnacentral$rna_type)
GUT_small_rna_rnacentral$rna_type <- ifelse( is.na(GUT_small_rna_rnacentral$annotation), "uncharacterized hypothetical RNA", GUT_small_rna_rnacentral$rna_type)

## Join GUT miRNA results (all predicted novel/known without filter)
GUT_predicted_mature_novel <- read_delim(file.path(basepath, "csv_files", "GUT_predicted_mature_novel.csv"), 
                                         delim = ",", 
                                         escape_double = FALSE, 
                                         comment = "#", 
                                         trim_ws = TRUE,
                                         show_col_types = FALSE)
names(GUT_predicted_mature_novel)[5] <- "total_read_count"
GUT_predicted_mature_novel$rna_type <- "R. prolixus novel miRNA"
GUT_predicted_mature_novel$class <- "miRNA"
GUT_predicted_mature_novel_to_plot <- GUT_predicted_mature_novel %>%
  select("rna_type", "class", "total_read_count")
GUT_predicted_mature_novel_to_plot$"total_read_count" <- sapply(GUT_predicted_mature_novel_to_plot$"total_read_count", as.numeric)

GUT_predicted_mature_known <- read_delim(file.path(basepath, "csv_files", "GUT_predicted_mature_known.csv"), 
                                         delim = ",", 
                                         escape_double = FALSE, 
                                         comment = "#", 
                                         trim_ws = TRUE,
                                         show_col_types = FALSE)
names(GUT_predicted_mature_known)[5] <- "total_read_count"
colnames(GUT_predicted_mature_known)[which(names(GUT_predicted_mature_known) == "total read count")] <- "total_read_count"
GUT_predicted_mature_known$rna_type <- "R. prolixus conserved miRNA"
GUT_predicted_mature_known$class <- "miRNA"
GUT_predicted_mature_known_to_plot <- GUT_predicted_mature_known %>%
  select("rna_type", "class", "total_read_count")
GUT_predicted_mature_known_to_plot$"total_read_count" <- sapply(GUT_predicted_mature_known_to_plot$"total_read_count", as.numeric)

# Bind miRNA + other smallRNA
GUT_small_rna_rnacentral['class'] <- "Other small RNA"
GUT_small_rna_rnacentral_to_plot <- GUT_small_rna_rnacentral %>%
  select("rna_type", "class", "total_read_count", "taxonomy", "annotation")
GUT_to_plot <- bind_rows( GUT_small_rna_rnacentral_to_plot,
                          bind_rows( 
                            GUT_predicted_mature_novel_to_plot, 
                            GUT_predicted_mature_known_to_plot)
)
# View(GUT_to_plot)

################################################################################
# HEMOLIMPH
################################################################################
HEM_small_rna <- read_delim(file.path(basepath, "csv_files", "HEM_other_small_rna_blast_results.tab"), 
                            "\t", 
                            col_names = FALSE, 
                            escape_double = FALSE, 
                            trim_ws = TRUE, 
                            show_col_types = FALSE)
colnames(HEM_small_rna) <- c("read", "rna_central", "identity", "alignment_length", "mismatches", "gap_opens", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score")
HEM_small_rna <- HEM_small_rna %>% distinct(read, .keep_all = TRUE)
HEM_small_rna <- HEM_small_rna %>%
  separate(read, c("read", "total_read_count"), "_x") # Insert read counts by spliting read column
HEM_small_rna$"total_read_count" <- sapply(HEM_small_rna$"total_read_count", as.numeric)
HEM_small_rna_rnacentral <- merge(HEM_small_rna, rnacentral_id_taxonomy_annotation, by=c("rna_central","rna_central")) 
temp <- str_split_fixed(as.character(HEM_small_rna_rnacentral$taxonomy)," ", 2) 
HEM_small_rna_rnacentral['genus'] <- temp[,1]
HEM_small_rna_rnacentral$rna_type <- ""
HEM_small_rna_rnacentral$rna_type <- ""
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)ribosomal", HEM_small_rna_rnacentral$annotation), "rRNA", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)rRNA", HEM_small_rna_rnacentral$annotation), "rRNA", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)tRNA", HEM_small_rna_rnacentral$annotation), "tRNA", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)spliceosomal", HEM_small_rna_rnacentral$annotation), "spliceosomal RNA", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)small nucleolar", HEM_small_rna_rnacentral$annotation), "snoRNA", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)signal recognition", HEM_small_rna_rnacentral$annotation), "srpRNA", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)antisense", HEM_small_rna_rnacentral$annotation), "asRNA", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)Nuclear RNase P", HEM_small_rna_rnacentral$annotation), "Nuclear RNase P", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)transfer", HEM_small_rna_rnacentral$annotation), "tRNA", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)rNase MRP", HEM_small_rna_rnacentral$annotation), "RNase MRP", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)sRNA", HEM_small_rna_rnacentral$annotation), "Bacterial small RNA", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)RybB", HEM_small_rna_rnacentral$annotation), "RybB RNA", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)Glm", HEM_small_rna_rnacentral$annotation), "GlmY RNA", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)CsrB", HEM_small_rna_rnacentral$annotation), "CsrB MRP", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)Hfq", HEM_small_rna_rnacentral$annotation), "Hfq binding sRNA", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)TPP riboswitch", HEM_small_rna_rnacentral$annotation), "TPP riboswitch", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)7sk", HEM_small_rna_rnacentral$annotation), "7SK RNA", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)nhaA-I", HEM_small_rna_rnacentral$annotation), "nhaA-I RNA", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)thermometer", HEM_small_rna_rnacentral$annotation), "RNA thermometer", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)KDPG-aldolase", HEM_small_rna_rnacentral$annotation), "KDPG-aldolase RNA motif", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)uncharacterized", HEM_small_rna_rnacentral$annotation), "uncharacterized hypothetical RNA", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)GcvB", HEM_small_rna_rnacentral$annotation), "GcvB RNA", HEM_small_rna_rnacentral$rna_type)
HEM_small_rna_rnacentral$rna_type <- ifelse( is.na(HEM_small_rna_rnacentral$annotation), "uncharacterized hypothetical RNA", HEM_small_rna_rnacentral$rna_type)

## Join HEM miRNA results (all predicted novel/known without filter)
HEM_predicted_mature_novel <- read_delim(file.path(basepath, "csv_files", "HEM_predicted_mature_novel.csv"), 
                                         delim = ",", 
                                         escape_double = FALSE, 
                                         comment = "#", 
                                         trim_ws = TRUE, 
                                         show_col_types = FALSE)
names(HEM_predicted_mature_novel)[5] <- "total_read_count"
HEM_predicted_mature_novel$rna_type <- "R. prolixus novel miRNA"
HEM_predicted_mature_novel$class <- "miRNA"
HEM_predicted_mature_novel_to_plot <- HEM_predicted_mature_novel %>%
  select("rna_type", "class", "total_read_count")
HEM_predicted_mature_novel_to_plot$"total_read_count" <- sapply(HEM_predicted_mature_novel_to_plot$"total_read_count", as.numeric)

HEM_predicted_mature_known <- read_delim(file.path(basepath, "csv_files", "HEM_predicted_mature_known.csv"), 
                                         delim = ",",
                                         escape_double = FALSE, 
                                         comment = "#", trim_ws = TRUE,
                                         show_col_types = FALSE)
names(HEM_predicted_mature_known)[5] <- "total_read_count"
colnames(HEM_predicted_mature_known)[which(names(HEM_predicted_mature_known) == "total read count")] <- "total_read_count"
HEM_predicted_mature_known$rna_type <- "R. prolixus conserved miRNA"
HEM_predicted_mature_known$class <- "miRNA"
HEM_predicted_mature_known_to_plot <- HEM_predicted_mature_known %>%
  select("rna_type", "class", "total_read_count")
HEM_predicted_mature_known_to_plot$"total_read_count" <- sapply(HEM_predicted_mature_known_to_plot$"total_read_count", as.numeric)

# Bind miRNA + other smallRNA
HEM_small_rna_rnacentral['class'] <- "Other small RNA"
HEM_small_rna_rnacentral_to_plot <- HEM_small_rna_rnacentral %>%
  select("rna_type", "class", "total_read_count", "taxonomy", "annotation")
HEM_to_plot <- bind_rows( HEM_small_rna_rnacentral_to_plot,
                          bind_rows( 
                            HEM_predicted_mature_novel_to_plot, 
                            HEM_predicted_mature_known_to_plot)
)
# View(HEM_to_plot)


################################################################################
# SALIVARY GLAND
################################################################################
GLA_small_rna <- read_delim(file.path(basepath, "csv_files", "GLA_other_small_rna_blast_results.tab"), 
                            "\t", 
                            col_names = FALSE, 
                            escape_double = FALSE, 
                            trim_ws = TRUE,
                            show_col_types = FALSE)
colnames(GLA_small_rna) <- c("read", "rna_central", "identity", "alignment_length", "mismatches", "gap_opens", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score")
GLA_small_rna <- GLA_small_rna %>% distinct(read, .keep_all = TRUE)
GLA_small_rna <- GLA_small_rna %>%
  separate(read, c("read", "total_read_count"), "_x") # Insert read counts by spliting read column
GLA_small_rna$"total_read_count" <- sapply(GLA_small_rna$"total_read_count", as.numeric)
GLA_small_rna_rnacentral <- merge(GLA_small_rna, rnacentral_id_taxonomy_annotation, by=c("rna_central","rna_central")) 
temp <- str_split_fixed(as.character(GLA_small_rna_rnacentral$taxonomy)," ", 2) 
GLA_small_rna_rnacentral['genus'] <- temp[,1]
GLA_small_rna_rnacentral$rna_type <- ""
GLA_small_rna_rnacentral$rna_type <- ""
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)ribosomal", GLA_small_rna_rnacentral$annotation), "rRNA", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)rRNA", GLA_small_rna_rnacentral$annotation), "rRNA", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)tRNA", GLA_small_rna_rnacentral$annotation), "tRNA", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)spliceosomal", GLA_small_rna_rnacentral$annotation), "spliceosomal RNA", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)small nucleolar", GLA_small_rna_rnacentral$annotation), "snoRNA", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)signal recognition", GLA_small_rna_rnacentral$annotation), "srpRNA", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)antisense", GLA_small_rna_rnacentral$annotation), "asRNA", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)Nuclear RNase P", GLA_small_rna_rnacentral$annotation), "Nuclear RNase P", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)transfer", GLA_small_rna_rnacentral$annotation), "tRNA", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)rNase MRP", GLA_small_rna_rnacentral$annotation), "RNase MRP", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)sRNA", GLA_small_rna_rnacentral$annotation), "Bacterial small RNA", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)RybB", GLA_small_rna_rnacentral$annotation), "RybB RNA", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)Glm", GLA_small_rna_rnacentral$annotation), "GlmY RNA", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)CsrB", GLA_small_rna_rnacentral$annotation), "CsrB MRP", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)Hfq", GLA_small_rna_rnacentral$annotation), "Hfq binding sRNA", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)TPP riboswitch", GLA_small_rna_rnacentral$annotation), "TPP riboswitch", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)7sk", GLA_small_rna_rnacentral$annotation), "7SK RNA", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)nhaA-I", GLA_small_rna_rnacentral$annotation), "nhaA-I RNA", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)thermometer", GLA_small_rna_rnacentral$annotation), "RNA thermometer", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)KDPG-aldolase", GLA_small_rna_rnacentral$annotation), "KDPG-aldolase RNA motif", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)uncharacterized", GLA_small_rna_rnacentral$annotation), "uncharacterized hypothetical RNA", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( grepl("(?i)GcvB", GLA_small_rna_rnacentral$annotation), "GcvB RNA", GLA_small_rna_rnacentral$rna_type)
GLA_small_rna_rnacentral$rna_type <- ifelse( is.na(GLA_small_rna_rnacentral$annotation), "uncharacterized hypothetical RNA", GLA_small_rna_rnacentral$rna_type)

## Join GLA miRNA results (all predicted novel/known without filter)
GLA_predicted_mature_novel <- read_delim(file.path(basepath, "csv_files", "GLA_predicted_mature_novel.csv"), 
                                         delim = ",", 
                                         escape_double = FALSE, 
                                         comment = "#", 
                                         trim_ws = TRUE,
                                         show_col_types = FALSE)
names(GLA_predicted_mature_novel)[5] <- "total_read_count"
GLA_predicted_mature_novel$rna_type <- "R. prolixus novel miRNA"
GLA_predicted_mature_novel$class <- "miRNA"
GLA_predicted_mature_novel_to_plot <- GLA_predicted_mature_novel %>%
  select("rna_type", "class", "total_read_count")
GLA_predicted_mature_novel_to_plot$"total_read_count" <- sapply(GLA_predicted_mature_novel_to_plot$"total_read_count", as.numeric)

GLA_predicted_mature_known <- read_delim(file.path(basepath, "csv_files", "GLA_predicted_mature_known.csv"), 
                                         delim = ",", 
                                         escape_double = FALSE, 
                                         comment = "#", 
                                         trim_ws = TRUE,
                                         show_col_types = FALSE)
names(GLA_predicted_mature_known)[5] <- "total_read_count"
colnames(GLA_predicted_mature_known)[which(names(GLA_predicted_mature_known) == "total read count")] <- "total_read_count"
GLA_predicted_mature_known$rna_type <- "R. prolixus conserved miRNA"
GLA_predicted_mature_known$class <- "miRNA"
GLA_predicted_mature_known_to_plot <- GLA_predicted_mature_known %>%
  select("rna_type", "class", "total_read_count")
GLA_predicted_mature_known_to_plot$"total_read_count" <- sapply(GLA_predicted_mature_known_to_plot$"total_read_count", as.numeric)

# Bind miRNA + other smallRNA
GLA_small_rna_rnacentral['class'] <- "Other small RNA"
GLA_small_rna_rnacentral_to_plot <- GLA_small_rna_rnacentral %>%
  select("rna_type", "class", "total_read_count", "taxonomy", "annotation")
GLA_to_plot <- bind_rows( GLA_small_rna_rnacentral_to_plot,
                          bind_rows( 
                            GLA_predicted_mature_novel_to_plot, 
                            GLA_predicted_mature_known_to_plot)
)
# View(GLA_to_plot)


## Remove temporary dataframe
rm(temp)
################################################################################
# Building THE PLOT
################################################################################
Gut_pie_data <- data.frame(GUT_to_plot$taxonomy, GUT_to_plot$total_read_count)
colnames(Gut_pie_data) <- c("Species", 'total_read_count')
Gut_pie_data <- within(Gut_pie_data, {
  f <- total_read_count < 100 # threshold 100 reads
  Species[f] <- 'Other Species'
})
Gut_pie_data <- Gut_pie_data[,-grep("f",colnames(Gut_pie_data))]
Gut_pie_data$Species <- Gut_pie_data$Species %>% replace_na('Not identified')
Gut_pie_data <- Gut_pie_data %>%
  group_by(Species) %>%
  summarize(abundance = sum(as.numeric(total_read_count, na.rm = TRUE)))
Gut_pie_data <- Gut_pie_data %>%
  mutate(per=paste0(round(abundance/sum(abundance)*100,2),"%"))
Gut_pie_data <- Gut_pie_data %>%
  mutate(label_abundance=paste0(abundance, " (",round(abundance/sum(abundance)*100,2),"%)"))
Gut_pie_data$tissue <- "Gut"
#View(Gut_pie_data)
########################################
Hem_pie_data <- data.frame(HEM_to_plot$taxonomy, HEM_to_plot$total_read_count)
colnames(Hem_pie_data) <- c("Species", 'total_read_count')
Hem_pie_data <- within(Hem_pie_data, {
  f <- total_read_count < 100 # threshold 100 reads
  Species[f] <- 'Other Species'
})
Hem_pie_data <- Hem_pie_data[,-grep("f",colnames(Hem_pie_data))]
Hem_pie_data$Species <- Hem_pie_data$Species %>% replace_na('Not identified')
Hem_pie_data <- Hem_pie_data %>%
  group_by(Species) %>%
  summarize(abundance = sum(as.numeric(total_read_count, na.rm = TRUE)))
Hem_pie_data <- Hem_pie_data %>%
  mutate(per=paste0(round(abundance/sum(abundance)*100,2),"%"))
Hem_pie_data <- Hem_pie_data %>%
  mutate(label_abundance=paste0(abundance, " (",round(abundance/sum(abundance)*100,2),"%)"))
Hem_pie_data$tissue <- "Hemolymph"
#View(Hem_pie_data)
#######################################
Gla_pie_data <- data.frame(GLA_to_plot$taxonomy, GLA_to_plot$total_read_count)
colnames(Gla_pie_data) <- c("Species", 'total_read_count')
Gla_pie_data <- within(Gla_pie_data, {
  f <- total_read_count < 100 # threshold 100 reads
  Species[f] <- 'Other Species'
})
Gla_pie_data <- Gla_pie_data[,-grep("f",colnames(Gla_pie_data))]
Gla_pie_data$Species <- Gla_pie_data$Species %>% replace_na('Not identified')
Gla_pie_data <- Gla_pie_data %>%
  group_by(Species) %>%
  summarize(abundance = sum(as.numeric(total_read_count, na.rm = TRUE)))
Gla_pie_data <- Gla_pie_data %>%
  mutate(per=paste0(round(abundance/sum(abundance)*100,2),"%"))
Gla_pie_data <- Gla_pie_data %>%
  mutate(label_abundance=paste0(abundance, " (",round(abundance/sum(abundance)*100,2),"%)"))
Gla_pie_data$tissue <- "Salivary Gland"
#View(Gla_pie_data)

All_pie <- rbind(rbind(Gla_pie_data, Gut_pie_data), Hem_pie_data )
#View(All_pie)

write.table(
  All_pie,
  file.path(basepath, "csv_files/All_pie_abundance_species.tab"),
  append = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE)


################################################################################
# THE PLOT:  Diversidade de Small RNAs
################################################################################
plot <- ggplot(
  All_pie, 
  aes( x = reorder(Species, -abundance), 
       fill   = Species, 
       group  = Species, 
       weight = abundance
  )
) +
  geom_bar(alpha=.75) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x  ),
                labels = trans_format("log10", math_format(10^.x  ))) +
  scale_fill_discrete() +
  labs(
    x = "Species", 
    y = "Mapped reads abundance across Species (%)"
  ) +
  #coord_flip() +
  geom_text(
    aes(
      y     = abundance,
      label = label_abundance),
    color   = "black",
    vjust   = -.4,
    # hjust   = -.15,
    size    = 2.5) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "italic"),
    axis.text.x = element_text(angle = 60, vjust = .96, hjust=1),
    # axis.title.x=element_blank(),
    # axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  ) +
  facet_wrap(
    vars(tissue),
    scales = "free",
    ncol = 1,
    strip.position = "left"
  )
plot
ggsave(file = file.path(basepath, "Plots", "Figure_Extra_The abundance of Rhodnius prolixus small RNAs against RNA central.pdf"), plot, width = 280, height = 360, units = "mm")
