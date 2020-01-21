# Smart-seq3

This repository contains the scripts and pipelines used to process and analyse Smart-seq3 libraries, as described in Hagemann-Jensen et al. 2020. 

We here provide the code to perform the following steps, that are expanded upon in the dedicated sub-folders.

## 1) Processing of Smart-seq3 data with zUMI. 
We show how fastq files are efficiently processed to BAM files in a manner that simultaneously distinguishes 5' from internal reads, and error-corrects both cell barcodes and molecular barcodes.

## 2) Scripts to reconstruct RNA molecules based on the zUMI prepared BAM files.
Using python script stitcher.py we in silico reconstruct RNA molecules based on the read pair alignments in the zUMI generated BAM files. Note that for RNA reconstruction, paired-end sequencing data is required. This step results in a new BAM file where each entry is a reconstructed molecule.

## 4) Scripts to assign reconstructed RNA molecules to allelic origins.
Stand-alone scripts that assigns ...

## 4) Scripts to assign reconstructed RNA molecules to transcript isoforms.
Using a couple of python scripts, we assign each RNA molecule to a set of compatible isoforms (including unique assignments). The resulting assignments are reported in tab-delimited text files.

## 5) Notebooks.
Here we post notebooks that show the analysis workflows for selected analyses from Hagemann-Jensen et al. as R or Python Jupyter notebooks.

