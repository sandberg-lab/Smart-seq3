# stitcher.py

_stitcher.py_ reconstructs molecules from Smart-seq3 data processed with zUMIs => 2.6.0 https://github.com/sdparekh/zUMIs and outputs a .sam file which can be used for further analysis.

## System Requirements for stitcher.py

_stitcher.py_  is a python3 script with dependencies:
```
pandas
numpy
pysam
joblib
pygtrie
pyinterval
```
No further installation is needed.

## Usage

stitcher.py [-h] [--i input] [--o output] [--g gtf] [--t threads]
                   [--cells cells] [--contig contig] [-v]

**arguments:**
```
  -h, --help       show this help message and exit
  --i input        Input .bam file
  --o output       Output .sam file
  --g gtf          gtf file with gene information
  --t threads      Number of threads
  --cells cells    List of cell barcodes to stitch molecules (text file, one cell barcode per line).
  --contig contig  Restrict stitching to contig
  -v, --version    show program's version number and exit
```


## Output 

_stitcher.py_ writes its results to a .sam file as the reads are being processed. Some of the fields have a slightly different interpretation than usual. The algorithm does not handle insertions at the moment, and remove those before constructing the molecule. The query name is in the format "cell:gene:umi". The D character in the CIGAR string indicates missing coverage. The MAPQ is always 255. 

The .sam file also contain many additional custom tags:

```
NR : Number of reads used to stitch.
ER : Number of reads covering an exon.
IR : Number of reads covering an intron.
BC : Cell barcode (same as zUMIs).
XT : Gene barcode (same as ZUMIs < 2.6.0).
UB : UMI (same as zUMIs).
EL : Locations of read ends (strand dependent).
NC : If there is a conflict in the reconstruction, the number of conflicting intervals.
IL : If there is a conflict in the reconstruction, the intervals where there is a conflict.
```

## Example 

```
python3 stitcher.py --i smartseq3_file.bam --o smartseq3_molecules.sam --g Mus_musculus.GRCm38.91.chr.clean.gtf --t 10 --contig chr1 --cells cells.txt
```
