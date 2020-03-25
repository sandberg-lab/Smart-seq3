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
portion
rpy2
```
No further installation is needed.

## Usage

stitcher.py [-h] [--i input] [--o output] [--g gtf] [--counts counts] 
            [--t threads] [--cells cells] [--contig contig] [v]

**arguments:**
```
  -h, --help       show this help message and exit
  --i input        Input .bam file
  --o output       Output .sam file
  --g gtf          gtf file with gene information
  --counts counts  zUMIs .rds file with read counts 
  --t threads      Number of threads
  --cells cells    List of cell barcodes to stitch molecules (text file, one cell barcode per line).
  --contig contig  Restrict stitching to contig
  -v, --version    show program's version number and exit
```
## Input

_stitcher.py_ takes .bam file processed with zUMIs => 2.6.0, together with a gtf file and the zUMIs .rds expression file, to reconstruct molecules for the genes in the gtf file.

As optional parameter, the user can specify the number of threads used for parallelization, the cells to process, and restrict the reconstruction to a given contig.

## Output 

_stitcher.py_ writes its results to a .sam file as the reads are being processed. Some of the fields have a slightly different interpretation than usual. The algorithm does not handle insertions at the moment, and remove those before stitching the molecule. The query name is in the format "cell:gene:umi". The D character in the CIGAR string indicates missing coverage. The MAPQ is always 255. 

In some cases the UMI reads contain conflicting information regarding the splicing of the molecule. This could either be due to differences in mapping or due to UMI collisions. In those cases, _stitcher.py_ prefer the unspliced option and writes two additional tags which contain the number of intervals which have conflicting information and the intervals themselves. 

The .sam file contain many additional custom tags:

```
NR : Number of reads used to stitch.
ER : Number of reads covering an exon.
IR : Number of reads covering an intron.
BC : Cell barcode (same as zUMIs).
XT : Gene barcode (same as ZUMIs < 2.6.0).
UB : UMI (same as zUMIs).
EL : Locations of read ends (strand dependent).
NC : Conflict in the reconstruction, the number of conflicting bases.
IL : Conflict in the reconstruction, the intervals where there is a conflict. Written as start1,end1,start2,end2,...
```

## Example 

```
python3 stitcher.py --i smartseq3_file.bam --o smartseq3_molecules.sam --counts zUMIs_output/expression/Smartseq3.dgecounts.rds --g Mus_musculus.GRCm38.91.chr.clean.gtf --t 10 --contig chr1 --cells cells.txt
```
