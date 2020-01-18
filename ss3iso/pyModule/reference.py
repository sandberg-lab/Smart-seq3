#!/usr/bin/python
# Developer: Ping Chen
# Contact: ping.chen@ki.se
# Date: 2020-01-10
# Version: 0.1.3

import re
import os
import subprocess
import pysam
import pandas as pd
from collections import defaultdict
import numpy as np
import multiprocessing as mp
from functools import partial
import pybedtools
import glob
import warnings
from .informative_reads import *

def _ref_transcript_struc(df, total_n_ex, gene_id):
    
    coordinates = df.groupby(by="blockCount").apply(lambda x: ';'.join(list(x.apply(lambda y: '%s-%s' %(y[1],y[2]), axis=1).values)))
    ex_idx = list(coordinates.index)
    ex_idx.sort()
    coordinates = coordinates[ex_idx]
    junc = ['%s^%s' %(ex_idx[ii], ex_idx[ii+1]) for ii in range(len(ex_idx)-1)]
    out = [gene_id, ','.join(map(str,ex_idx)), ','.join(list(coordinates.values)), ','.join(junc), str(total_n_ex)]
    
    return pd.Series(out)

def _build_gene_ref(indir, outdir, gene_info, gene):

    print(gene)
    os.system('echo "%s" >> %s/_log' %(gene, outdir))

    rcds = '\t'.join(map(str,gene_info.loc[gene].values))
    
    obj = geneObj(None, None, indir)
    obj.get_exon_coordinates(rcds.strip())
    obj.outdir = outdir
    
    curr_df = sm.query('gene=="%s"' %(obj.gene)).iloc[:,[0,1,2,3,5]]
    if curr_df.shape[0]==0: return None 
    curr_df.to_csv('%s/_%s' %(outdir, obj.gene), sep="\t", index=False, header=False)
    
    gene_bed = pybedtools.BedTool('%s/_%s' %(outdir, obj.gene))
    tmp = gene_bed.intersect(obj.ex_bed, wa=True, wb=True)
    intersect = tmp.to_dataframe()
    res = intersect.groupby(by="score").apply(_ref_transcript_struc, obj.ex_bed.to_dataframe()['score'].max(), obj.gene)
    
    return res

   
def build_reference(conf, indir):
    
    outdir = '%s/.reference' %(indir)
    if not os.path.exists(outdir): os.makedirs(outdir)
    
    gene_info = pd.read_table('%s/gene.gff' %(indir), header=None, index_col=None, sep='\t')
    gene_info.index = gene_info.apply(lambda x: x[8].split(';')[0].split('=')[1], axis=1)
    
    if conf['annotation']['gtf_source'] == 'ensembl':
        t_idx = 3
    else:
        t_idx = 2
        
    exon_ref = pd.read_table('%s/exon.gff' %(indir), header=None, index_col=None, sep="\t")
    genes = pd.DataFrame([val.split(';')[1].split('=')[1] for val in exon_ref.iloc[:,8].values])
    trans = pd.DataFrame([val.split(';')[t_idx].split('=')[1] for val in exon_ref.iloc[:,8].values])
    annot = pd.concat([genes,trans], axis=1)
    
    global sm
    sm = exon_ref.iloc[:,[0,3,4,6]]
    sm = pd.concat([sm, annot], axis=1)
    sm.columns = ['chrom','start','end','strand','gene','transcript']
    sm['start'] = sm['start'] - 1
    
    keptfiles = glob.glob('%s/keptReads/*/*.csv' %(indir))
    genes = ['_'.join(val.split('/')[-1].split('_')[:-2]) for val in keptfiles]
    chrom_dict = {'_'.join(val.split('/')[-1].split('_')[:-2]): val.split('/')[-2] for val in keptfiles}
    os.system('> %s/_log' %(outdir))
    
    pool = mp.Pool(processes=int(conf['expression']['nproc']))
    func = partial(_build_gene_ref, indir, outdir, gene_info)
    results = pool.map(func, genes, chunksize=1)
    
    filtered = [item for item in results if item is not None]
    
    out_df = pd.concat(filtered, axis=0)
    out_df.index.name = 'Transcript'
    out_df.reset_index(inplace=True)
    out_df.columns = ['Transcript','Gene','Exon_Index','Exon_Loc','Junction','Total_n_exons']
    
    out_chr = pd.DataFrame([chrom_dict[gene] for gene in out_df['Gene'].values], columns=['chrom'])
    ref_iso = pd.concat([out_df, out_chr], axis=1)
    
    p = subprocess.Popen('rm -rf %s/../.tmp/*' %(outdir), shell=True)
    (output, err) = p.communicate()  
    
    p = subprocess.Popen('rm -rf %s/../.tempDir' %(outdir), shell=True)
    (output, err) = p.communicate()  
    
    os.system('mkdir %s/../.tempDir' %(outdir))
    
    p = subprocess.Popen('rm -rf %s' %(outdir), shell=True)
    (output, err) = p.communicate()  
    
    return ref_iso