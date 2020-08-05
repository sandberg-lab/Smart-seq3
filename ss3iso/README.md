# Welcome to ss3iso 

<img src="https://github.com/sandberg-lab/Smart-seq3/blob/master/ss3iso/isoform_reconstruction.png" alt="UMI reads" width="600"/>

ss3iso is a Python pipeline developed for isoform reconstruction of UMI-linking fragments from Smart-seq3. For detailed information, please read our paper [Single-cell RNA counting at allele- and isoform-resolution using Smart-seq3](https://www.biorxiv.org/content/10.1101/817924v1).

ss3iso uses [zUMIs](https://github.com/sdparekh/zUMIs) output BAM tagged with corrected cell and UMI barcodes as input. The pipeline requires GTF annotations (Ensembl, RefSeq or Gencode) and needs to be specified by **gtf_source** in configuration file.

## Dependencies

Make sure the following softwares and Python packages are installed before running ss3iso.

```
Python3
tabix
bedtools (v2.26.0)
samtools

optparse (python module)
glob (python module)
configparser (python module)
re (python module)
pybedtools (python module)
subprocess (python module)
pysam (python module)
pandas (python module)
collections (python module)
numpy (python module)
multiprocessing (python module)
functools (python module)
```

## Installation

Checkout ss3iso repository to your prefered folder on a computing server using following command. No futher installation is needed. 

``` git clone https://github.com/sandberg-lab/Smart-seq3/ss3iso.git ```

## Usage

Execute ss3iso pipeline using the following command line.
```
python ss3_isoform.py -i [path/to/inputBAM] -c [path/to/configuration file] -e [experiment] -o [path/to/output directory] -p [number of processes] -s [species] -P -Q
```

Options:
```
-i, --inputBAM: input ZUMIs BAM path
-c, --config: the required pipeline configuration file
-e, --experiment: the name of the experiment/study
-o, --outputDir: the output directory
-p, --process: the number of processes for parallel computing (default: 8)
-s, --species: the species under study (default: hg38)
-P, --Preprocess: run preprocessing on input BAM
-Q, --Quantification: run isoform reconstruction and quantification
```

## Changelog

## Getting help
If you have any questions and suggestions on our pipeline, feel free to contact us by email (ping.chen@ki.se, rickard.sandberg@ki.se).

