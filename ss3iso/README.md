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
-i, --inputBAM: input ZUMIs BAM path. Note: Use '*filtered.tagged.Aligned.out.bam.ex.featureCounts.UBfix.sort.bam' generated by zUMIs. Every read should have a UB:Z tag.
-c, --config: the required pipeline configuration file
-e, --experiment: the name of the experiment/study
-o, --outputDir: the output directory
-p, --process: the number of processes for parallel computing (default: 8)
-s, --species: the species under study (default: hg38)
-P, --Preprocess: run preprocessing on input BAM
-Q, --Quantification: run isoform reconstruction and quantification
```

Example contents in the input BAM:
```
NB502120:154:HVG7JBGXB:2:21104:11500:9869       163     1       14409   3       85M65S  =       14692   410     GCTCAGTTCTTTATTGATTGGTGTGCCGTTTTCTCTGGAAGCCTCTTAAGAACACAGTGGCGCAGGCTGGGTGGAGCCGTCCCCCATGAAGTACAGGCAGACAAGTCCCCGCCCCAGCTGTGTGGCCTCAAGCCAGCCTTCCACTCCTTG  AAAAAEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEAEEAEEEEE/EEAEEEEEEEAEAAAAE<AEEEAEEA/EE<EAE/AEE/AEAAEEEE6//  NH:i:2  HI:i:1  AS:i:206        nM:i:2  BX:Z:AAGCCGTTTGAACGCT   BC:Z:AAGCCGTTTGAACGCT   UX:Z:TAATCTCT   XS:Z:Unassigned_Ambiguity       UB:Z:TAATCTCT
NB502120:154:HVG7JBGXB:1:13112:25990:9712       163     1       14414   255     150M    =       14749   602     GTTCTTTATTGATTGGTGTGCCGTTTTCTCTGGAAGCCTCTTAAGAACACTGTGGCGCAGGCTGGGTGGAGCCGTCCCCCCATGGAGCACAGGCAGACAGAAGTCCCCGCCCCAGCTGTGTGGCCTCAAGCCAGCCTTCCGCTCCTTGAA  AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAE/EEEEE<AAAAEEEEEE/EEAEEEEAEEEEEEEE<EEEEAA<EEEEEE/EAEEE/EAEEAAE<AA//6A/EA<<6/6<6  NH:i:1  HI:i:1  AS:i:275        nM:i:1  BX:Z:AAGCCGTTGAGGTTAG   BC:Z:AAGCCGTTGAGGTTAG   UX:Z:GCCAAGGG   XS:Z:Assigned3  XN:i:1  XT:Z:ENSG00000227232    UB:Z:GCCAAGGG
NB502120:154:HVG7JBGXB:4:23406:11214:8076       163     1       14414   3       98M2D51M1S      =       14692   405     GTTCTTTATTGATTGGTGTGCCGTTTTCTCTGGAAGCCTCTTAAGAACACAGTGGCGCAGGCTGGGTGGAGCCGTCCCCCCATGAAGTACAGGCAGACAAGTCCCCGCCCCAGCTGTGTGGCCTCAAGCCAGCCTTCCGCTCCTTGAAGC  6AAAAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEE<E/EEEEEEEE6E/EEAEEEEEEAEEAA<<AAA/AEE<AEAEAEEEEEE<E//E//EAAEEEE//66<AEAAEEEEE6AE<A/6A6/AA<A</<6A//  NH:i:2  HI:i:1  AS:i:262        nM:i:3  BX:Z:AAGCCGTTGAGGTTAG   BC:Z:AAGCCGTTGAGGTTAG   UX:Z:GTGACTCT   XS:Z:Assigned1  XN:i:1  XT:Z:ENSG00000227232    UB:Z:GTGACTCT
NB502120:154:HVG7JBGXB:1:12205:16478:4384       163     1       14414   3       80M70S  =       14692   405     GTTCTTTATTGATTGGTGTGCCGTTTTCTCTGGAAGCCTCTTAAGAACACAGTGGCGCAGGCTGGGTGGAGCCGTCCCCCATGAAGTACAGGCAGACAAGTCCCCGCCCCAGCTGTGTGGCCTCAAGCCAGCCTTCCACTCCTTGAAGCT  AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/EEEEEEEEEEEAEEEEEAEEEEEEE/EEEEEE6AAEAEEA/EEEEEEEEEEEEEEAEEEAAEEEA<E<EEEA<AAEEEEA  NH:i:2  HI:i:1  AS:i:201        nM:i:2  BX:Z:AAGCCGTTTGAACGCT   BC:Z:AAGCCGTTTGAACGCT   UX:Z:TAATCTCT   XS:Z:Assigned1  XN:i:1  XT:Z:ENSG00000227232    UB:Z:TAATCTCT
NB502120:154:HVG7JBGXB:1:21202:15138:15924      163     1       14414   255     42M5I102M1S     =       14683   405     GTTCTTTATTGATTGGTGTGCCATTTTCTCTGGAAGCCTCTTTAGAGAAGAACACAGTGGCGCAGGCTGGGTGGAGCCGTCCCCCCATGGAGCACAGGCAGACAGAAGTCCCCCCCCCAGCTGTGTGGCCTCAGGCCAGCCTTCCGCTCC  AAAAAAEEEEEEEEEAEEEEEEEEEEEAE//EEEEEEEEEE<EAEEEEEE/EEEEEEEEEEEEEEEEEEEE6EAAEEEEEEEEEEEEEAE<AEAE/EEEE/<EEEEEEEEAEE/<<<</EA6/6<////AAAEE/<EE/E</6AAAAAAA  NH:i:1  HI:i:1  AS:i:216        nM:i:9  BX:Z:GCATGTCTGAACCTGT   BC:Z:GCATGTCTGAACCTGT   UX:Z:GCATCTGG   XS:Z:Assigned3  XN:i:1  XT:Z:ENSG00000227232    UB:Z:GCATCTGG
NB502120:154:HVG7JBGXB:4:11511:5522:12140       163     1       14414   3       80M70S  =       14668   382     GTTCTTTATTGATTGGTGTGCCGTTTTCTCTGGAAGCCTCTTAAGAACACAGTGGCGCAGGCTGGGTGGAGCCGTCCCCCATGAAGTACAGGCAGACAAGTCCCCGCCCCAGCTGTGTGGCCTCAAGCCAGCCTTCCACTCCTTGAAGCT  AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE6E6AAAEA<EEAAE<E<EEEEEEEEEEAEEE<E<EAAA<A<AAEEEA<  NH:i:2  HI:i:1  AS:i:204        nM:i:1  BX:Z:TACCGTCTTAGCAAGC   BC:Z:TACCGTCTTAGCAAGC   UX:Z:AACGGTGT   XS:Z:Assigned1  XN:i:1  XT:Z:ENSG00000227232    UB:Z:AACGGTGT
NB502120:154:HVG7JBGXB:4:21509:2326:6303        163     1       14414   255     150M    =       14692   405     GTTCTTTATTGATTGGTGTGCCGTTTTCTCTGGAAGCCTCTTAAGAACACAGTGGCGCAGGCTGGGTGGAGCCGTCCCCCCATGGAGCACAGGCAGACAGAAGTCCCCGCCCCAGCTGTGTGGCCTCAAGCCAGCCTTCCGCTCCTTGAA  AAAAAEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAEEAEEEAEEEEEEAEEEEEAEEEEEAEE<E<AAAEE/AEA/E<EEEEE<EEEEEEAEEEEAEEAEEEE<<6EE  NH:i:1  HI:i:1  AS:i:275        nM:i:0  BX:Z:TTCCGTTCCCTCTTCA   BC:Z:TTCCGTTCCCTCTTCA   UX:Z:TGCATCTC   XS:Z:Assigned1  XN:i:1  XT:Z:ENSG00000227232    UB:Z:TGCATCTC
NB502120:154:HVG7JBGXB:4:23502:2511:18043       163     1       14414   255     150M    =       14692   405     GTTCTATATTGATTGGTGTGCCGTTTTCTCTGGAAGCCTCTTAAGAACACAGTGGCGCAGGCTGGGTGGAGCCGACCCCCCATGGAGCACAGGCAGACACAAGTCCCCGCCCCAGCTGTGTGGCCTCAAGCCAGCCTACCGCTCCTTGAA  AAAAA/EEEE/EE//6AE/EEEAEEEE/EAEAEAEEEEEE/EEE6EEE/AE//AAA/AEE/E6EAEE6EEEAEE/<<EEAEA/A<<///A/<<EAEEEA/EA//A/6<<A/<E/<AE<EEEA/EE//AA/E<//EA//<A<E<<A/6<AA  NH:i:1  HI:i:1  AS:i:265        nM:i:5  BX:Z:TTCCGTTCCCTCTTCA   BC:Z:TTCCGTTCCCTCTTCA   UX:Z:TGCATCTC   XS:Z:Assigned1  XN:i:1  XT:Z:ENSG00000227232    UB:Z:TGCATCTC
NB502120:154:HVG7JBGXB:1:11308:5428:8170        163     1       14419   3       150M    =       14757   604     TTATTGATTGGTGTGCCGTTTTCTCTGGAAGCCTCTTAAGAACACAGTGGCGCAGGCTGGGTGGAGCCGTCCCCCCATGGAGCACAGGCAGACAGAAGTCCCCGCCCCAGCTGTGTGGCCTCAAGCCAGCCTACCGCTCCATGAAGCAGG  A/AAAE6///AAE//E/E//E/EAEEE/E66EEEE<EEEE//EAEE6EE/EEEA/AEE/AE//E/EEEEA/EEEEE/AEE/EEEE//EEEE/E<EE/EEEEEEAEEEEE/E</AA<EAEE<EEA/E/AE<E/AE6<A/AE/////<<///  NH:i:2  HI:i:1  AS:i:268        nM:i:4  BX:Z:CGCAAGAACGGTTGTT   BC:Z:CGCAAGAACGGTTGTT   UX:Z:GCTGGGCG   XS:Z:Assigned3  XN:i:1  XT:Z:ENSG00000227232    UB:Z:GCTGGGCG
NB502120:154:HVG7JBGXB:1:21206:19630:11590      163     1       14433   3       79M2D70M1S      =       14692   386     GCCGTTTTCTCTGGAAGCCTCTTAAGAACACAGTGGCGCAGGCTGGGTGGAGCCGTCCCCCCATGAAGTACAGGCAGACAAGTCCCCGCCCCAGCAGTGTGGCCTCAAGCCAGCCTTCCGCTCCTTGAAGCTGGTCTCCACACAGTGCTA  AAAAAEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEA<EEEEEEA<EEEEE/EEEAEAAEE/EEEEA/E/EEEEEAEA6A/<AEE/<AA/AEE<E<EE<6AAE/AE<<EEE<AAAEAA//<A/A  NH:i:2  HI:i:1  AS:i:260        nM:i:4  BX:Z:TTGGAACCCGGTTGTT   BC:Z:TTGGAACCCGGTTGTT   UX:Z:GCAGATTC   XS:Z:Assigned1  XN:i:1  XT:Z:ENSG00000227232    UB:Z:GCAGATTC
NB502120:154:HVG7JBGXB:1:21202:9431:15912       163     1       14453   255     149M1S  =       14692   365     CTTAAGAACACTGTGGCGCAGGCTGGGTGGAGCCGTCCCCCCATGGAGCACAGGCAGACAGAAGTCCCCGCCCCAGCTGTGTGGCCTCAAGCCAGCCTTCCGCTCCTTGAAGCTGGTCTCCACACAGTGCTGGTTCCGTCACCCCCTCCC  AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEAEEEEEEEEEAEAEEEE</EEEEEEAAAEEEEE/EEAEA<EAEEEEAEEEEEE/A/EE<E/E<AEEE/A</AAEEE/A<<<6  NH:i:1  HI:i:1  AS:i:269        nM:i:2  BX:Z:AAGCCGTTCCACATAG   BC:Z:AAGCCGTTCCACATAG   UX:Z:GCATCCTG   XS:Z:Assigned1  XN:i:1  XT:Z:ENSG00000227232    UB:Z:GCATCCTG
NB502120:154:HVG7JBGXB:3:21405:1669:3957        163     1       14453   255     149M1S  =       14692   365     CTTAAGAACACTGTGGCGCAGGCTGGGTGGAGCCGTCCCCCCATGGAGCACAGGCAGACAGAAGTCCCCGCCCCAGCTGTGTGGCCTCAAGCCAGCCTTCCGCTCCTTGAAGCTGGTCTCCACACAGTGCTGGTTCCGTCACCCCCTCCC  AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEE/EEEEEAEE6EEEEEEEEEEEAEEEEEEEEE/EE//EE/EEAE/A/<AEE/EEEA<EA/E/EE<<EAE/EEAEEE<A<AAEEE//<E<<A</<EEE6/A<EEE<AAA  NH:i:1  HI:i:1  AS:i:269        nM:i:2  BX:Z:AAGCCGTTCCACATAG   BC:Z:AAGCCGTTCCACATAG   UX:Z:GCATCCTG   XS:Z:Assigned1  XN:i:1  XT:Z:ENSG00000227232    UB:Z:GCATCCTG
NB502120:154:HVG7JBGXB:1:22107:18366:16253      163     1       14453   3       59M2D90M1S      =       14692   366     CTTAAGAACACAGTGGCGCAGGCTGGGTGGAGCCGTCCCCCCATGAAGTACAGGCAGACAAGTCCCCGCCCCAGCTGTGTGGCCTCAAGCCAGCCTTCCGCTCCTTGAAGCTGGTCTCCACACAGTGCTAGTTCCATCACCCCCTCCCAG  AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEE/EEEEAEEEEEEEEEAEEEEEEEEEEEEEAE6EE6EEE<<EEAEEEE/EAEEEEAEE<EEEA/E/A/<EAEEEEEAEEE<<66<<<<6<AAA/  NH:i:2  HI:i:1  AS:i:258        nM:i:5  BX:Z:TCCAAGTCGAACCTGT   BC:Z:TCCAAGTCGAACCTGT   UX:Z:GGCTCTCG   XS:Z:Assigned1  XN:i:1  XT:Z:ENSG00000227232    UB:Z:GGCTCTCG
NB502120:154:HVG7JBGXB:3:21601:21748:12734      163     1       14453   3       149M1S  =       14728   403     CTTAAGAACACTGTGGCGCAGGCTGGGTGGAGCCGTCCCCCCATGGAGCACAGGCAGACAGAAGTCCCCGCCCCAGCTGTGTGGCCTCAAGCCAGCCTTCCGCTCCTTGAAGCTGGTCTCCACACAGTGCTGGTTCCGTCACCCCCTCCC  AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEAEAEEEEEEEEEEEEEEEEE<E/EEEEEEEEEA<<EAEEEEAEEEEEEEEEE/AEE6EEEEEEEEEEEEE/E/AEEEAEEEEEAAEAEE<<E<A/AEAEEA<A<A/E<A  NH:i:2  HI:i:1  AS:i:273        nM:i:1  BX:Z:TCCAAGTCGAACCTGT   BC:Z:TCCAAGTCGAACCTGT   UX:Z:ACTTGGGT   XS:Z:Assigned1  XN:i:1  XT:Z:ENSG00000227232    UB:Z:ACTTGGGT
NB502120:154:HVG7JBGXB:2:13107:5473:12860       163     1       14455   255     150M    =       14692   364     TAAGAACACTGTGGCGCAGGCTGGGTGGAGCCGTCCCCCCATGGAGCACAGGCAGACAGAAGTCCCCGCCCCAGCTGTGTGGCCTCAAGCCAGCCTTCCGCTCCTTGAAGCTGGTCTCCACACAGTGCTGGTTCCGTCACCCCCTCCCAA  A6AAAEEEEEEEEEAEEAEEEAEEEAEE6E6AEEEEEEEEEEEEEEEEEAEEEEE/EEEEEEEAA</6A<A/AEAAEEE6EEEEEEE6EEEAEEEA<AEEEEEE/<EAEE<AEAAEAAAAE<A/AA<AA/EEE/EEEE/<6<<<</<<//  NH:i:1  HI:i:1  AS:i:273        nM:i:1  BX:Z:CACCTAACCAGATTCG   BC:Z:CACCTAACCAGATTCG   UX:Z:GTATCGAC   XS:Z:Assigned1  XN:i:1  XT:Z:ENSG00000227232    UB:Z:GTATCGAC
NB502120:154:HVG7JBGXB:1:23309:15624:7432       163     1       14455   255     150M    =       14692   362     TAAGAACACTGTGGCGCAGGCTGGGTGGAGCCGTCCCCCCATGGAGCACAGGCAGACAGAAGTCCCCGCCCCAGCTGTGTGGCCTCAAGCCAGCCTTCCGCTCCATGAAGCTGGTCTCCACACAGTGCTGGTTCCGTCCCCCCCTCCCAA  A<A66EEEEEEEE6A6/6E6EEAE/6EE/EEEEEA<EEE/E/EEE6EEE/EEEAEEAE/AEE/AAE//EE<E//<EE/A<AAEEEEAEEEA6/EEE6/<EA/E//A<EE/EEAEEE</<EE/EA/AEEAEE<EE//A//A6/6//AEA/E  NH:i:1  HI:i:1  AS:i:265        nM:i:4  BX:Z:TTAGGCCACCATCCAA   BC:Z:TTAGGCCACCATCCAA   UX:Z:GAGTTGAG   XS:Z:Assigned1  XN:i:1  XT:Z:ENSG00000227232    UB:Z:GAGTTGAG
```

## Getting help
If you have any questions and suggestions on our pipeline, feel free to contact us by email (ping.chen@ki.se, rickard.sandberg@ki.se).

