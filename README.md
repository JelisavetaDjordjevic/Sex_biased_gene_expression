# Sex_biased_gene_expression

Repository for the collected scripts used in study: 
https://doi.org/10.1101/2021.01.23.427895 

*DATA*

Contains 
1. Timema californicum read counts- **timema_counts.csv**
2. dN/dS estimates- **all_stages_dnds.txt**
3. Drosophila melanogaster read counts- **droso_counts.csv**
4. Accession numbers of Timema californicum reads reads- **reads_accession_numbers.txt**

*SCRIPTS* 

*Quality control of the reads (removing adapters, trimming),
mapping reads to the genome and counting the reads per gene*

1. **QC_mapping_counting.sh**

*To analyze Timema data*

Get sex-biased genes and classify genes into categories
2. At the hatchling stage- 
   **hatclings_edger.R** 
   
3. At the juvenile stage-
   **juveniles_edger.R**
   
4. At the adult stage- 
   **adults_edger.R**

To make a heatmap of sex-biased genes
5. **heatmap_timema.R**

To plot dNdS values for different sex- bias gene categories and to statistically compare between gene categories groups

6.  **DnDs.R**

To calculate stage specificity of gene expression (Tau index) for sex-bias categories, and statistics to compare tau values between sex-bias gene categories

7. **tau_calculation.R**

*To analyze Drosophila data*

Get sex-biased genes and classify genes into categories
1. At the larval stage- 
   **droso_larva_edger.R** 
   
2. At the juvenile stage
   **droso_pupa_edger.R**
   
3. At the adult stage- 
   **droso_adult_edger.R**

To get venn diagram plot of three developemental stages, and to statistically test the intersection

4. **droso_venn_diagram.R**
