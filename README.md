# Rhodnius prolixus miRNA project

## Conda environment (https://www.anaconda.com/products/individual)

### Creating the environment and installing the tools

* seqkit (version 0.16.1)
* mirDeep2 (version 0.1.3)

`conda create -n bioinfo python=3.6`

`conda activate bioinfo`

`conda install -c bioconda seqkit`

`conda install mirdeep2`

`conda update mirdeep2`

`conda install -c conda-forge biopython`


## Trimming fasta headers to avoid head length problems using mirDeep2
`seqkit seq -i RprolixusV48.fa > out.fasta`

## Parsing trimmed fastq files to fit them as fasta inputs for mirDeep2 using the in-house script *mirDepp2_fit_fastq_to_fasta.py*

`python3 mirDepp2_fit_fastq_to_fasta.py trimmed/RP1G.fastq RP1G_input_for_mirDeep2.fasta RP1`

`python3 mirDepp2_fit_fastq_to_fasta.py trimmed/RP2H.fastq RP2H_input_for_mirDeep2.fasta RP2`

`python3 mirDepp2_fit_fastq_to_fasta.py trimmed/RPGland.fastq RPGland_input_for_mirDeep2.fasta RPG`

## Create reference indexes

`bowtie-build --large-index --bmax 16777216 --dcv 256 --threads 4 reference/RprolixusV48.fa reference/RprolixusV48`

## Map reads to the Rhodnius prolixus genome (version 48) using mirDeep2 combined with a set of *Ecdysozoa* filtered miRNA downloaded from miRBase

### RPgland

`mirdeep2/bin/mapper.pl RPGland_input_for_mirDeep2.fasta -c -p reference/RprolixusV48 -t RPGland.arf -o 4 -i -j -n -q -v`

`mirdeep2/bin/miRDeep2.pl RPGland_input_for_mirDeep2.fasta reference/RprolixusV48.fa RPGland.arf mirBase_mature_Ecdysozoa_for_meerDeep2.fasta none none`

### RP1G
`mirdeep2/bin/mapper.pl RP1G_input_for_mirDeep2.fasta -c -p reference/RprolixusV48 -t RP1G.arf -o 4 -i -j -n -q -v`

`mirdeep2/bin/miRDeep2.pl RP1G_input_for_mirDeep2.fasta reference/RprolixusV48.fa RP1G.arf mirBase_mature_Ecdysozoa_for_meerDeep2.fasta none none`

### RP2H

`mirdeep2/bin/mapper.pl RP2H_input_for_mirDeep2.fasta -c -p reference/RprolixusV48 -t RP2H.arf -o 4 -i -j -n -q -v`

`mirdeep2/bin/miRDeep2.pl RP2H_input_for_mirDeep2.fasta reference/RprolixusV48.fa RP2H.arf mirBase_mature_Ecdysozoa_for_meerDeep2.fasta none none`


### How to view the results?

Look at the .html files inside the individual experiments to them in an interactive browser
The same results are provided in a .csv file that can be opened in software like Microsoft/Excel or OpenOffice/Calc

* Sequences from RNA central: tax_string:"triatominae" and so_rna_type_name:"NcRNA"

`blastn -task blastn-short -db ../RnaCentral/tax_stringtriatominae_AND_so_rna_type_nameNcRNA  -query RPGland_mapped.fasta  -max_target_seqs 5 -max_hsps 1 -evalue 1e-2 -perc_identity 80 -num_threads 4 -outfmt 6 -out RPGland_mapped_blast_results.tab`

`blastn -task blastn-short -db ../RnaCentral/tax_stringtriatominae_AND_so_rna_type_nameNcRNA  -query RP1G_mapped.fasta  -max_target_seqs 5 -max_hsps 1 -evalue 1e-2 -perc_identity 80 -num_threads 4 -outfmt 6 -out RP1G_mapped_blast_results.tab`

`blastn -task blastn-short -db ../RnaCentral/tax_stringtriatominae_AND_so_rna_type_nameNcRNA  -query RP2H_mapped.fasta  -max_target_seqs 5 -max_hsps 1 -evalue 1e-2 -perc_identity 80 -num_threads 4 -outfmt 6 -out RP2H_mapped_blast_results.tab`

