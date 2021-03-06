# Rhodnius prolixus miRNA project

![Rhodnius prolixus miRNA project](https://github.com/waldeyr/Rhodnius_prolixus_miRNA_project/raw/main/Diagrams/Pipeline.png)

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

## Parsing miRBase *Ecdysozoa* sequences to fit them as fasta inputs for mirDeep2 using the in-house script *adjustMirBaseEcdysozoaFastaHeaderFormirDeep2.py*

```Bash
python3 adjustMirBaseEcdysozoaFastaHeaderFormirDeep2.py
```

## Create reference indexes

```Bash
bowtie-build --large-index --bmax 16777216 --dcv 256 --threads 4 reference/RprolixusV48.fa reference/RprolixusV48
```

## Map reads to the Rhodnius prolixus genome (version 48) using mirDeep2 combined with a set of *Ecdysozoa* filtered miRNA downloaded from miRBase

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

## Parsing the mirDeep2 resulting ARF files to Fasta using an in-house script

```Bash
python3 arf_to_fasta.py

```

## Get a tab separated file with miRBase *Ecdysozoa* Species names using the in-house script *getMirBaseEcdysozoaSpecies.py*

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

