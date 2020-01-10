#CAST and C57 allele-specific expression
Here we provide tools to classify molecules to their allele of origin for the CAST X C57 F1 mouse cells.

First, sequencing data should be processed using zUMIs from fastq files to aligned bam files and UMI count tables.
Of note, the genome positions with strain-specific variation should be masked with N to avoid a mapping bias towards the reference allele.
[SNPsplit](https://github.com/FelixKrueger/SNPsplit) can be used to generate the N-masked genome fasta file.  


 `zUMIs-master.sh -y mouse_cross.yaml`

Based on the zUMIs output, you can run the allele-specific expression script.
It requires only the config file used for zUMIs and a VCF file of CAST specific SNPs.
In this repository, we provide the VCF file used for the publication analyses. This file contains CAST/EiJ strain specific SNPs, obtained from the
mouse genome project dbSNP version 142 and filtered for variants clearly observed in existing CAST/EiJ x C57/Bl6J F1 data.

 `Rscript get_variant_overlap_CAST.R --help`
 `Rscript get_variant_overlap_CAST.R --yaml mouse_cross.yaml --vcf CAST.SNPs.validated.vcf.gz`

For users with a working zUMIs installation, the script does not require additional dependencies.
The output contains files for both directly assigned molecules and total UMI counts broken down by the observed gene-wise allele-fractions.
