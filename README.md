# Sex_biased_gene_expression

Repository for the collected scripts used in study: 
https://doi.org/10.1101/2021.01.23.427895 

*DATA*

Contains 
1. Timema californicum read counts- **timema_counts.csv**
2. dN/dS estimates- **all_stages_dnds.txt**
3. Drosophila melanogaster read counts- **droso_counts.csv**

*SCRIPTS* 

*To analyze Timema data*

Get sex-biased genes and classify genes into categories
1. At the hatchling stage- 
   **hatclings_edger.R** 
   
2. At the juvenile stage
   **juveniles_edger.R**
   
3. At the adult stage- 
   **adults_edger.R**

To plot dNdS values for different sex- bias gene categories and to statistically compare between gene categories groups

4.  **DnDs.R**

To calculate stage specificity of gene expression (Tau index) for sex-bias categories, and statistics to compare tau values between sex-bias gene categories

5. **tau_calculation.R**

*To analyze Drosophila data*

1. At the larval stage- 
   **droso_larva_edger.R** 
   
2. At the juvenile stage
   **droso_pupa_edger.R**
   
3. At the adult stage- 
   **droso_adult_edger.R**
