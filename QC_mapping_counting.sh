#Download reference genome and annotation (GFF)
#accession numbers for Timema californicum reads in the table https://github.com/JelisavetaDjordjevic/Sex_biased_gene_expression/tree/main/Data/reads_accession_numbers.txt

#Genome and annotation for T. californicum available at:  https://doi.org/10.5281/zenodo.5636226 



###Quality control of the reads and trimming

## Remove the adapters
#### Cutadapt
module load UHTS/Quality_control/cutadapt/2.3

Trim_dir="adapt_trimmed"
mkdir $Trim_dir

for i in <path/to/reads>/*R1.fq.gz ; do 
R1_file=`echo $i`
	R2_file=`echo $i | sed 's/R1.fq.gz/R2.fq.gz/'`
	TR1a=`echo $R1_file | sed 's/.*\///' | sed 's/.fq.gz/_Trimmed.fq.gz/' `
	TR2a=`echo $R2_file | sed 's/.*\///'| sed 's/.fq.gz/_Trimmed.fq.gz/' `
	TR1=`echo $Trim_dir"/"$TR1a`
	TR2=`echo $Trim_dir"/"$TR2a`
	echo $R1_file
	echo $R2_file
	echo $TR1
	echo $TR2

cutadapt \
			-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
			-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
			-o $TR1 -p $TR2 \
			$R1_file $R2_file 

done

### trimmomatic ### 
module load UHTS/Analysis/trimmomatic/0.36

qTrim_dir="Adapt_trimmed_qtr"
mkdir $qTrim_dir

for i in $Trim_dir/*R1_Trimmed.fq.gz; do
	R1_file=`echo $i`
	R2_file=`echo $i | sed 's/R1_Trimmed.fq.gz/R2_Trimmed.fq.gz/'`
	qTR1a=`echo $R1_file | sed 's/.*\///' | sed 's/Trimmed.fq.gz/ADQT.fq/' `
	qTR2a=`echo $R2_file | sed 's/.*\///' | sed 's/Trimmed.fq.gz/ADQT.fq/' `
	qTR1=`echo $qTrim_dir"/"$qTR1a`
	qTR2=`echo $qTrim_dir"/"$qTR2a`
	
	qTR1a_UP=`echo $R1_file | sed 's/.*\///' | sed 's/Trimmed.fq.gz/ADQT_UP.fq/' `
	qTR2a_UP=`echo $R2_file | sed 's/.*\///' | sed 's/Trimmed.fq.gz/ADQT_UP.fq/' `
	qTR1_UP=`echo $qTrim_dir"/"$qTR1a_UP`
	qTR2_UP=`echo $qTrim_dir"/"$qTR2a_UP`
	
	echo $R1_file
	echo $R2_file
	echo $qTR1
	echo $qTR2
	echo $qTR1_UP
	echo $qTR2_UP
	trimmomatic PE -phred33 $R1_file $R2_file $qTR1 $qTR1_UP $qTR2 $qTR2_UP LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:80

done
### Remove unpaired files in Adapt_trimmed_qtr ( rm Adapt_trimmed_qtr/*UP.fq) ####

###Build index for star

module load UHTS/Aligner/STAR/2.6.0c
mkdir indexes

STAR --runThreadN 12 \
     --runMode genomeGenerate \
     --genomeDir </path/to/indexes> \
     --genomeFastaFiles </path/to/genome/2_Tcm_b3v08.fa> \
     --sjdbGTFfile </path/to/annotation/2_Tcm_b3v08.max.func.ipr.gff> \
     --sjdbGTFtagExonParentTranscript Parent \
     --sjdbOverhang 99 \
     --genomeChrBinNbits 15  --limitGenomeGenerateRAM 79000000000

###map with star
module load UHTS/Aligner/STAR/2.6.0c

for i in <path/to/trimmed_reads>/*R1_ADQT.fq.gz; do
R1_file=`echo $i`
R2_file=`echo $i | sed 's/R1_ADQT.fq.gz/R2_ADQT.fq.gz/'`

STAR --genomeDir </path/to/indexes> \
     --readFilesIn $R1_file $R2_file  --readFilesCommand zcat --outSAMtype BAM Unsorted --runThreadN 12 \
     --outFileNamePrefix ${R1_file%%_R1_ADQT.fq.gz}

done

##to check mapping and summarize results use multiqc (it use any log files)
module load UHTS/Analysis/MultiQC/1.3
multiqc .

## Sort bam files 

module load UHTS/Analysis/samtools/1.8

for i in <path/to/mapped/>*out.bam ; do
        file=`echo $i`
        samtools sort -n $file > `echo $file | sed 's/.out.bam/.out.sorted.bam/'`
done

### Use "Maker_gff_to_HTseq_gff.py" to change gff file for HTseq to count
#download "Maker_gff_to_HTseq_gff.py" from : https://github.com/DarrenJParker/Useful_scripts/blob/master/Maker_gff_to_HTseq_gff.py

python3 Maker_gff_to_HTseq_gff.py -i </path/to/annotation.gff> \
	-o 2_Tcm_b3v08_forHTSeq

### Counting gene units (including introns) with HTseq 

module load UHTS/Analysis/HTSeq/0.9.1

for i in </path/to/sorted_bam>/*sorted.bam ; do
file=`echo $i`

htseq-count -f bam -r name -s reverse -t gene -i ID -m union --nonunique none $i \
</path/to/fixedgff/2_Tcm_b3v08_forHTSeq.gff> > `echo $file | sed 's/Aligned.out.sorted.bam/.counts/'`

done

### Convert HTseq counts for edgeR (using HTSeq_to_edgeR.py )
#download HTSeq_to_edgeR.py from: https://github.com/DarrenJParker/Useful_scripts/blob/master/HTSeq_to_edgeR.py

python3 HTSeq_to_edgeR.py -i /path/to/counts -o for_edgeR

