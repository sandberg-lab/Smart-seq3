# Smart-seq3

This repository contains the scripts and pipelines used to process and analyse Smart-seq3 libraries, as described in Hagemann-Jensen et al. 2020. 

We here provide the code to perform the following steps, that are expanded upon in the dedicated sub-folders.

### 1) Processing of Smart-seq3 data with zUMIs. 
We show how fastq files are efficiently processed to BAM files in a manner that simultaneously distinguishes 5' from internal reads, and error-corrects both cell barcodes and molecular barcodes using [zUMIs](https://github.com/sdparekh/zUMIs).

First, you should obtain raw fastq files *without demultiplexing*, as the data will be processed in a pooled fashion. When running the [bcl2fastq](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) conversion, be sure to keep index read fastq files.

Example for a dual-index, 150 bp PE run: 
`bcl2fastq --use-bases-mask Y150N,I8,I8,Y150N --no-lane-splitting --create-fastq-for-index-reads -R /mnt/storage1/NextSeqNAS/191011_NB502120_0154_AHVG7JBGXB`

Next, prepare your config file in [YAML format for zUMIs](https://github.com/sdparekh/zUMIs/wiki/Usage#setup-using-the-yaml-config-file). The UMI sequence needs to be correctly extracted from 5' reads in Smart-seq3. These will always be the first Illumina read and are recognized by our unique 11bp tag sequence. Thus, you need to set the following settings:

```
file1:
    name: /mnt/storage2/temp_workdir/Undetermined_S0_L003_R1_001.fastq.gz
    base_definition:
      - cDNA(23-150)
      - UMI(12-19)
    find_pattern: ATTGCGCAATG
```

You can find an [example YAML file here](https://github.com/sandberg-lab/Smart-seq3/blob/master/allele_level_expression/mouse_cross.yaml).

Note that we advise caution when using STARs 2-pass mapping mode, as we have observed some spurious novel splice junctions being used that may distort molecule reconstructions.

### 2) Scripts to reconstruct RNA molecules based on the zUMIs prepared BAM files.
Using our python script [*stitcher.py*](https://github.com/sandberg-lab/Smart-seq3/tree/master/stitcher) we in silico reconstruct RNA molecules based on the read pair alignments in the zUMIs generated BAM files. Note that for RNA reconstruction, paired-end sequencing data is required. This step results in a new BAM file where each entry is a reconstructed molecule.

https://github.com/sandberg-lab/Smart-seq3/tree/master/stitcher

### 3) Scripts to assign reconstructed RNA molecules to allelic origins.
We provide a stand-alone Rscript that assigns molecules to their allele of origin.

https://github.com/sandberg-lab/Smart-seq3/tree/master/allele_level_expression

### 4) Scripts to assign reconstructed RNA molecules to transcript isoforms.
Using a [couple of python scripts](https://github.com/sandberg-lab/Smart-seq3/tree/master/ss3iso), we assign each RNA molecule to a set of compatible isoforms (including unique assignments). The resulting assignments are reported in tab-delimited text files.

https://github.com/sandberg-lab/Smart-seq3/tree/master/ss3iso

### 5) Notebooks.
Here we post notebooks that show the analysis workflows for selected analyses from Hagemann-Jensen et al. as R or Python Jupyter notebooks.

