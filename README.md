# *Rhodnius prolixus* miRNA project


![Rhodnius prolixus miRNA project](https://github.com/waldeyr/Rhodnius_prolixus_miRNA_project/blob/main/Diagrams/README.png)

This repository is related to the manuscript **"Insights Into the miRNome Profiles from the Salivary Glands, Gut and Hemolymph of the Chagas Disease Vector *Rhodnius prolixus*"**.
It contains scripts and some bioinformatics comands used to generate and anlyse data and plots.

## Abstract
microRNAs (miRNAs) are small non-coding RNAs that post-transcriptionally regulate gene expression by targeting mRNAs for degradation or translational repression. The growing interest in miRNAs over the last years has led to their description in numerous organisms, although there are no data on miRNAs in the Triatominae subfamily, the vectors of Trypanosoma cruzi, the causative agent of Chagas disease, which disables millions of people mainly in Latin America. Here, we provide the first insight into Rhodnius prolixus miRNome profiles (gut, hemolymph, and salivary glands tissues) using high-throughput sequencing. We identified, 52 mature miRNAs previously reported in Ecdysozoa taxa, including 39 ubiquitously expressed in the three tissues. In addition, 112, 73 and 78 novel miRNAs could be predicted in the gut, hemolymph, and salivary glands, respectively, from which 15 are shared among the three samples. The comprehensive resource of miRNAs from R. prolixus is available at NCBI (https://www.ncbi.nlm.nih.gov). To provide functional insight into miRNAs expression data, in silico analysis of the target genes of the top eight expressed miRNAs from salivary glands suggests R. prolixus may modulate the host’s gene expression at the post-transcriptional level. This study presents the first experimental evidence of miRNAs in a Triatominae species, disclosing these important regulatory molecules of vector biology. 


### Creating the environment and installing the tools

* [Conda environment](https://www.anaconda.com/products/individual)
* [seqkit](https://bioinf.shenwei.me/seqkit/) (version 0.16.1)
* [mirDeep2](https://www.mdc-berlin.de/content/mirdeep2-documentation) (version 0.1.3)
* [Biopython](https://biopython.org) (version 1.78)

```Bash
conda create -n bioinfo python=3.6
```

```Bash
conda activate bioinfo
```

```Bash
conda install -c bioconda seqkit
```

```Bash
conda install mirdeep2
```

```Bash
conda update mirdeep2
```

```Bash
conda install -c conda-forge biopython
```


## Trimming fasta headers to avoid head length problems using mirDeep2

```Bash
seqkit seq -i RprolixusV48.fa > out.fasta
```

## Parsing trimmed fastq files to fit them as fasta inputs for mirDeep2 using the in-house script *mirDepp2_fit_fastq_to_fasta.py*

```Bash
python3 mirDepp2_fit_fastq_to_fasta.py trimmed/RP1G.fastq RP1G_input_for_mirDeep2.fasta RP1
```

```Bash
python3 mirDepp2_fit_fastq_to_fasta.py trimmed/RP2H.fastq RP2H_input_for_mirDeep2.fasta RP2
```

```Bash
python3 mirDepp2_fit_fastq_to_fasta.py trimmed/RPGland.fastq RPGland_input_for_mirDeep2.fasta RPG
```

## Getting and parsing miRBase *Ecdysozoa* sequences to fit them as fasta inputs for mirDeep2 using the in-house script *adjustMirBaseEcdysozoaFastaHeaderFormirDeep2.py*

```Bash
python3 adjustMirBaseEcdysozoaFastaHeaderFormirDeep2.py
```

## Create reference indexes for the *R. prolixus* genome

```Bash
bowtie-build --large-index --bmax 16777216 --dcv 256 --threads 4 reference/RprolixusV48.fa reference/RprolixusV48
```

## Map reads to the *R. prolixus* genome (version 48) using mirDeep2 combined with a set of *Ecdysozoa* filtered miRNA downloaded from miRBase

### RPgland

```Bash
mirdeep2/bin/mapper.pl RPGland_input_for_mirDeep2.fasta -c -p reference/RprolixusV48 -t RPGland.arf -o 4 -i -j -n -q -v
```

```Bash
mirdeep2/bin/miRDeep2.pl RPGland_input_for_mirDeep2.fasta reference/RprolixusV48.fa RPGland.arf mirBase_mature_Ecdysozoa_for_meerDeep2.fasta none none
```

### RP1G
```Bash
mirdeep2/bin/mapper.pl RP1G_input_for_mirDeep2.fasta -c -p reference/RprolixusV48 -t RP1G.arf -o 4 -i -j -n -q -v
```

```Bash
mirdeep2/bin/miRDeep2.pl RP1G_input_for_mirDeep2.fasta reference/RprolixusV48.fa RP1G.arf mirBase_mature_Ecdysozoa_for_meerDeep2.fasta none none
```

### RP2H

```Bash
mirdeep2/bin/mapper.pl RP2H_input_for_mirDeep2.fasta -c -p reference/RprolixusV48 -t RP2H.arf -o 4 -i -j -n -q -v
```

```Bash
mirdeep2/bin/miRDeep2.pl RP2H_input_for_mirDeep2.fasta reference/RprolixusV48.fa RP2H.arf mirBase_mature_Ecdysozoa_for_meerDeep2.fasta none none
```

## Parsing the mirDeep2 resulting ARF files to Fasta using the in-house script *arf_to_fasta.py*

```Bash
python3 arf_to_fasta.py

```

## Getting and parsing to a tab separated file a miRBase *Ecdysozoa* Species names on a file using the in-house script *getMirBaseEcdysozoaSpecies.py*

```Bash
python3 getMirBaseEcdysozoaSpecies.py

```

## Identifying other small RNA (not miRNA) using Blast alignments of mapped reads against a set of small RNAs filtered for *Triatominae* downloaded from RNACentral

* Applied filter to download the sequences from RNA central: tax_string:"triatominae" and so_rna_type_name:"NcRNA"

```Bash
blastn -task blastn-short -db ../RnaCentral/tax_stringtriatominae_AND_so_rna_type_nameNcRNA  -query RPGland_mapped.fasta  -max_target_seqs 5 -max_hsps 1 -evalue 1e-2 -perc_identity 80 -num_threads 4 -outfmt 6 -out RPGland_mapped_blast_results.tab
```

```Bash
blastn -task blastn-short -db ../RnaCentral/tax_stringtriatominae_AND_so_rna_type_nameNcRNA  -query RP1G_mapped.fasta  -max_target_seqs 5 -max_hsps 1 -evalue 1e-2 -perc_identity 80 -num_threads 4 -outfmt 6 -out RP1G_mapped_blast_results.tab
```

```Bash
blastn -task blastn-short -db ../RnaCentral/tax_stringtriatominae_AND_so_rna_type_nameNcRNA  -query RP2H_mapped.fasta  -max_target_seqs 5 -max_hsps 1 -evalue 1e-2 -perc_identity 80 -num_threads 4 -outfmt 6 -out RP2H_mapped_blast_results.tab
```


## Filtering the valid results for miRNA (known and novel)

This step consists in run the R scripts `Filter_predicted_mature_known.R` and `Filter_predicted_mature_novel.R``. The criteria for filtering are the same described in the manuscript.

## To run these commands it is necessary to download the large files:

[Arf files](https://drive.google.com/drive/folders/1Ji4K8ozLo0RAAz9HTPWxXBO8aRr5IFKB?usp=share_link)
[CSV files](https://drive.google.com/drive/folders/1bxiD8OJxTuaKQGhhE7JkTPazWcRga-g2?usp=share_link)
[TXT files](https://drive.google.com/drive/folders/1YzTkRm1ZksuFG0UtFFgoFwNd_faVY5nb?usp=share_link)

The raw reads can be downloaded from NCBI from the links presented in the manuscript.