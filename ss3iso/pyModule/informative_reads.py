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

def get_exons(chrom, start, end, strand, gid, outdir):
    
    p = subprocess.Popen(['tabix', '%s/exon.sorted.gff.gz' %(outdir), '%s:%s-%s' %(chrom, start, end)], stdout=subprocess.PIPE)
    rcds = p.communicate()[0].decode("utf-8")
    
    rcds = [item for item in rcds.strip().split('\n') if re.search(gid, item)]
    
    tmpfile = '%s/.tempDir/_temp_%s.bed' %(outdir, gid)
    outF = open(tmpfile, "w")
    for rcd in rcds:
        items = rcd.split('\t')
        outF.write('\t'.join([items[0], str(int(items[3])-1), items[4], '.', '.', items[6]]))
        outF.write('\n')
    outF.close()
    
    os.system('sort -k1,1 -k2,2n %s | bedtools merge -s -d -1 -i - > %s/.tempDir/_temp_%s.merged.bed' %(tmpfile, outdir, gid))
    
    ordered_exons = pd.read_table('%s/.tempDir/_temp_%s.merged.bed' %(outdir, gid), header=None, index_col=None, sep="\t")
    
    return ordered_exons

class geneObj(object):
    
    def __init__(self, in_bam_uniq, in_bam_multi, outdir):
        
        self.strand_flags = {'+': [99, 147], '-': [83, 163]}
            
        self.gene = None
        self.exons = None
        self.ex_bed = None
        self.in_bam_uniq = in_bam_uniq
        self.in_bam_multi = in_bam_multi
        self.outdir = outdir
        self.chrom = None
        self.start = None
        self.end = None
        self.strand = None
        self.uniq_aligned_reads = None
        self.multi_aligned_reads = None
        self.uniq_r_bclist = None
        
    def get_exon_coordinates(self, gene):
        
        fds = gene.split('\t')
        gene_id = fds[-1].split(';')[0].split('=')[1]
        self.gene = gene_id
        
        exons = get_exons(fds[0], fds[3], fds[4], fds[6], gene_id, self.outdir)
        
        exon_idx = pd.DataFrame(list(range(1,exons.shape[0]+1)))
        exons = pd.concat([exons, exon_idx], axis=1)
        exons.to_csv('%s/.tempDir/_%s' %(self.outdir, self.gene), index=False, header=False, sep="\t")
        self.exons = exons
        self.ex_bed = pybedtools.BedTool('%s/.tempDir/_%s' %(self.outdir, self.gene))
        
        self.chrom = str(exons.iloc[0,0])
        self.start = np.min([exons.iloc[0,1], exons.iloc[-1,2]])
        self.end = np.max([exons.iloc[0,1], exons.iloc[-1,2]])
        self.strand = exons.iloc[0,3]
        
        return
    
    def get_aligned_reads(self, n_read_limit, passed_cells):
        
        samfile = pysam.AlignmentFile(self.in_bam_uniq, "rc")
        try:
            r_iterator = samfile.fetch(self.chrom, int(self.start), int(self.end))
        except:
            return None
        
        nreads = len([r_idx for r_idx, x in enumerate(r_iterator) if x.flag in self.strand_flags[self.strand]])
        if nreads > n_read_limit: return self.gene
              
        r_iterator = samfile.fetch(self.chrom, int(self.start), int(self.end))
        read_dict = {r_idx: _make_dict(x, self.chrom, self.strand, self.gene, r_idx) for r_idx, x in enumerate(r_iterator) if x.flag in self.strand_flags[self.strand] and list(filter(regx1.match, x.to_dict()['tags']))[0].replace('BC:Z:','') in passed_cells}
        samfile.close()
        
        df = [read_dict[r_idx]['r_blocks'] for r_idx in read_dict.keys()]
        
        if len(df) == 0: return None
        pd.concat(df, axis=0).to_csv('%s/.tempDir/_%s_reads_blocks.bed' %(self.outdir, self.gene), index=False, sep="\t", header=False)
        read_bed = pybedtools.BedTool('%s/.tempDir/_%s_reads_blocks.bed' %(self.outdir, self.gene))
        
        tmp = self.ex_bed.intersect(read_bed, wa=True, wb=True)
        if os.stat(tmp.fn).st_size == 0:
            return None
        
        intersect_all = tmp.to_dataframe()
        read_idx_list = list(set(intersect_all.iloc[:,9].values))
    
        ex_coord = ','.join(self.exons.apply(lambda x: '%s-%s' %(x[1],x[2]), axis=1).values)
        
        aligned_reads = [_make_list_aligned_reads2(r_idx, read_dict, intersect_all, ex_coord) for r_idx in read_idx_list]
        
        colnames = ['name', 'flag', 'ref_name', 'ref_pos', 'map_quality', 'cigar',
                    'next_ref_name', 'next_ref_pos', 'length', 'seq', 'qual', 'tags',
                    'read_mapped_position', 'geneid', 'Exon_Index', 'Category', 'BC', 'UB', 'exon_coordinates']
        self.uniq_aligned_reads = pd.DataFrame(aligned_reads, columns=colnames).drop_duplicates()
        self.uniq_r_bclist = list(set(self.uniq_aligned_reads.apply(lambda x: '%s+%s' %(x['BC'], x['UB']), axis=1).values))
        self.uniq_aligned_reads.insert(19, 'MapFlag', 'unique')
        
        return None
    
    def get_aligned_reads_from_multi(self, passed_cells):
        
        samfile = pysam.AlignmentFile(self.in_bam_multi, "rc")
        try:
            r_iterator = samfile.fetch(self.chrom, int(self.start), int(self.end))
        except:
            return None
        
        read_dict = {r_idx: _make_dict2(x, self.chrom, self.strand, self.gene, self.uniq_r_bclist, r_idx) for r_idx, x in enumerate(r_iterator) if x.flag in self.strand_flags[self.strand] and list(filter(regx1.match, x.to_dict()['tags']))[0].replace('BC:Z:','') in passed_cells}
        df = [read_dict[r_idx]['r_blocks'] for r_idx in read_dict.keys() if read_dict[r_idx] is not None]
        
        if len(df) == 0: return None
        pd.concat(df, axis=0).to_csv('%s/.tempDir/_%s_reads_blocks.bed' %(self.outdir, self.gene), index=False, sep="\t", header=False)
        read_bed = pybedtools.BedTool('%s/.tempDir/_%s_reads_blocks.bed' %(self.outdir, self.gene))
       
        tmp = self.ex_bed.intersect(read_bed, wa=True, wb=True)
        if os.stat(tmp.fn).st_size == 0:
            return None
        
        intersect_all = tmp.to_dataframe()
        read_idx_list = list(set(intersect_all.iloc[:,9].values))
        
        ex_coord = ','.join(self.exons.apply(lambda x: '%s-%s' %(x[1],x[2]), axis=1).values)
        aligned_reads = [_make_list_aligned_reads2(r_idx, read_dict, intersect_all, ex_coord) for r_idx in read_idx_list]
        
        colnames = ['name', 'flag', 'ref_name', 'ref_pos', 'map_quality', 'cigar',
                    'next_ref_name', 'next_ref_pos', 'length', 'seq', 'qual', 'tags',
                    'read_mapped_position', 'geneid', 'Exon_Index', 'Category', 'BC', 'UB', 'exon_coordinates']
        self.multi_aligned_reads = pd.DataFrame(aligned_reads, columns=colnames).drop_duplicates()
        self.multi_aligned_reads.insert(19, 'MapFlag', 'multi')
        
        return None

def _initialize_make_list_aligned():
    
    global my_read_dict
    global my_intersect_all
    global my_ex_coord
    
regx1 = re.compile("BC:Z:")
regx2= re.compile("UB:Z:")
def get_aligned_reads_mp(obj, nproc, passed_cells):
    
    global my_read_dict
    global my_intersect_all
    global my_ex_coord
    
    my_intersect_all = None
    my_read_dict = None
    my_ex_coord = None
    
    samfile = pysam.AlignmentFile(obj.in_bam_uniq, "rc")
    try:
        r_iterator = samfile.fetch(obj.chrom, int(obj.start), int(obj.end))
    except:
        return obj
    
    rcds = np.array([[r_idx, x.to_dict(), x.get_blocks()] for r_idx, x in enumerate(r_iterator) if x.flag in obj.strand_flags[obj.strand] and list(filter(regx1.match, x.to_dict()['tags']))[0].replace('BC:Z:','') in passed_cells])
    pool = mp.Pool(processes=nproc)
    func = partial(_make_dict_mp, obj.chrom, obj.strand, obj.gene)
    read_dict_list = pool.map(func, rcds, chunksize=1)
    pool.close()
    
    my_read_dict = {}
    tmp = [my_read_dict.update(elemt) for elemt in read_dict_list]
    df = [my_read_dict[r_idx]['r_blocks'] for r_idx in my_read_dict.keys()]
    samfile.close()
        
    if len(df) == 0: return obj
    pd.concat(df, axis=0).to_csv('%s/.tempDir/_%s_reads_blocks.bed' %(obj.outdir, obj.gene), index=False, sep="\t", header=False) 
    read_bed = pybedtools.BedTool('%s/.tempDir/_%s_reads_blocks.bed' %(obj.outdir, obj.gene))
        
    tmp = obj.ex_bed.intersect(read_bed, wa=True, wb=True)
    if os.stat(tmp.fn).st_size == 0:
            return obj
    
    my_intersect_all = tmp.to_dataframe()
    read_idx_list = list(set(my_intersect_all.iloc[:,9].values))
    
    my_ex_coord = ','.join(obj.exons.apply(lambda x: '%s-%s' %(x[1],x[2]), axis=1).values)
    
    pool = mp.Pool(processes=nproc, initializer=_initialize_make_list_aligned)
    aligned_reads = pool.map(_make_list_aligned_reads_mp, read_idx_list, chunksize=1)
    pool.close()
    
    colnames = ['name', 'flag', 'ref_name', 'ref_pos', 'map_quality', 'cigar',
                'next_ref_name', 'next_ref_pos', 'length', 'seq', 'qual', 'tags',
                'read_mapped_position', 'geneid', 'Exon_Index', 'Category', 'BC', 'UB', 'exon_coordinates']
    obj.uniq_aligned_reads = pd.DataFrame(aligned_reads, columns=colnames).drop_duplicates()
    obj.uniq_r_bclist = list(set(obj.uniq_aligned_reads.apply(lambda x: '%s+%s' %(x['BC'], x['UB']), axis=1).values))
    obj.uniq_aligned_reads.insert(19, 'MapFlag', 'unique')
        
    return obj

def get_aligned_reads_from_multi_mp(obj, nproc, passed_cells):
    
    global my_read_dict
    global my_intersect_all
    global my_ex_coord
    global my_uniq_r_bclist
    
    my_intersect_all = None
    my_read_dict = None
    my_ex_coord = None
    my_uniq_r_bclist = obj.uniq_r_bclist.copy()
        
    samfile = pysam.AlignmentFile(obj.in_bam_multi, "rc")
    try:
        r_iterator = samfile.fetch(obj.chrom, int(obj.start), int(obj.end))
    except:
        return obj
        
    rcds = np.array([[r_idx, x.to_dict(), x.get_blocks()] for r_idx, x in enumerate(r_iterator) if x.flag in obj.strand_flags[obj.strand] and list(filter(regx1.match, x.to_dict()['tags']))[0].replace('BC:Z:','') in passed_cells])
    pool = mp.Pool(processes=nproc)
    func = partial(_make_dict2_mp, obj.chrom, obj.strand, obj.gene)
    
    read_dict_list = pool.map(func, rcds, chunksize=1)
    pool.close()
    
    my_read_dict = {}
    tmp = [my_read_dict.update(elemt) for elemt in read_dict_list if elemt is not None]  # fast!!
    df = [my_read_dict[r_idx]['r_blocks'] for r_idx in my_read_dict.keys()]
    samfile.close()
        
    if len(df) == 0: return obj
    pd.concat(df, axis=0).to_csv('%s/.tempDir/_%s_multi_reads_blocks.bed' %(obj.outdir, obj.gene), index=False, sep="\t", header=False)
    read_bed = pybedtools.BedTool('%s/.tempDir/_%s_multi_reads_blocks.bed' %(obj.outdir, obj.gene))
       
    tmp = obj.ex_bed.intersect(read_bed, wa=True, wb=True)
    if os.stat(tmp.fn).st_size == 0:
            return obj
        
    my_intersect_all = tmp.to_dataframe()
    read_idx_list = list(set(my_intersect_all.iloc[:,9].values))
        
    my_ex_coord = ','.join(obj.exons.apply(lambda x: '%s-%s' %(x[1],x[2]), axis=1).values)
    
    pool = mp.Pool(processes=nproc, initializer=_initialize_make_list_aligned)
    aligned_reads = pool.map(_make_list_aligned_reads_mp, read_idx_list, chunksize=1)
    pool.close()
    
    colnames = ['name', 'flag', 'ref_name', 'ref_pos', 'map_quality', 'cigar',
                'next_ref_name', 'next_ref_pos', 'length', 'seq', 'qual', 'tags',
                'read_mapped_position', 'geneid', 'Exon_Index', 'Category', 'BC', 'UB', 'exon_coordinates']
    
    obj.multi_aligned_reads = pd.DataFrame(aligned_reads, columns=colnames).drop_duplicates()
    obj.multi_aligned_reads.insert(19, 'MapFlag', 'multi') 
    
    return obj

def gtf2exon(gtf, outdir, include_spikein=False):
    
    filename = '%s/exon.gff' %(outdir)
    outF = open(filename, "w")
    
    if include_spikein:
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                fds = line.strip().split('\t')
                if fds[0] in ['diySpike']:
                    if fds[2] == 'exon':
                        annot = fds[8].replace(' "', '=').replace('"; ',';').replace('";','').replace('; ',';').replace(' ','=').split(';')
                        annot = 'loc=%s:%s-%s:%s;%s' %(fds[0],fds[3],fds[4],fds[6],';'.join([annot[i] for i in [0,2,1]]))
                        outF.write('%s\t%s\n' %('\t'.join(fds[:8]), annot))
                else:
                    if fds[2] == 'exon':
                        annot = fds[8].replace(' "', '=').replace('"; ',';').replace('";','').replace('; ',';').replace(' ','=')
                        annot = 'loc=%s:%s-%s:%s;%s' %(fds[0],fds[3],fds[4],fds[6],annot)
                        outF.write('%s\t%s\n' %('\t'.join(fds[:8]), annot))
        outF.close()
    else:
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                fds = line.strip().split('\t')
                if fds[2] == 'exon':
                    annot = fds[8].replace(' "', '=').replace('"; ',';').replace('";','').replace('; ',';').replace(' ','=')
                    annot = 'loc=%s:%s-%s:%s;%s' %(fds[0],fds[3],fds[4],fds[6],annot)
                    outF.write('%s\t%s\n' %('\t'.join(fds[:8]), annot))
        outF.close()
    
    os.system('sort -k1,1 -k4,4n %s | bgzip > %s/exon.sorted.gff.gz' %(filename, outdir))
    os.system('tabix -p gff %s/exon.sorted.gff.gz' %(outdir))
    os.system('zless %s/exon.sorted.gff.gz | bedtools merge -i - -s -d -1 -c 1 -o count > %s/exon_merged.bed' %(outdir, outdir)) 
    
    return

def _make_dict(x, chrom, strand, gene, r_idx):
   
    print(r_idx)
    curr = x.to_dict()
    r_blocks  = pd.DataFrame(x.get_blocks(), columns=['start','end'])
    r_blocks.insert(0,'chr', chrom)
    r_blocks.insert(3,'strand', strand)
    r_blocks.insert(4,'rid', r_idx)
      
    curr['r_blocks'] = r_blocks
    curr['read_mapped_position'] = ','.join(r_blocks.apply(lambda x: '%s-%s' %(x[1],x[2]), axis=1).values)   # start: 0-based; end: 1-based
    curr['geneid'] = gene
    
    return curr

def _make_dict_mp(chrom, strand, gene, rcd):
   
    r_idx, curr_dict, block = rcd
   
    print(r_idx)
    
    r_blocks  = pd.DataFrame(block, columns=['start','end'])
    r_blocks.insert(0,'chr', chrom)
    r_blocks.insert(3,'strand', strand)
    r_blocks.insert(4,'rid', r_idx)
      
    curr_dict['r_blocks'] = r_blocks
    curr_dict['read_mapped_position'] = ','.join(r_blocks.apply(lambda x: '%s-%s' %(x[1],x[2]), axis=1).values)   # start: 0-based; end: 1-based
    curr_dict['geneid'] = gene
    
    return {r_idx: curr_dict}

def _make_dict2(x, chrom, strand, gene, uniq_r_bclist, r_idx):
    
    regx1 = re.compile("BC:Z:")
    regx2= re.compile("UB:Z:")
   
    print(r_idx)
    curr = x.to_dict()
    
    bc = list(filter(regx1.match, curr['tags']))[0].replace('BC:Z:','')
    ub = list(filter(regx2.match, curr['tags']))[0].replace('UB:Z:','')
    if '%s+%s' %(bc, ub) not in uniq_r_bclist: return None
    
    r_blocks  = pd.DataFrame(x.get_blocks(), columns=['start','end'])
    r_blocks.insert(0,'chr', chrom)
    r_blocks.insert(3,'strand', strand)
    r_blocks.insert(4,'rid', r_idx)
      
    curr['r_blocks'] = r_blocks
    curr['read_mapped_position'] = ','.join(r_blocks.apply(lambda x: '%s-%s' %(x[1],x[2]), axis=1).values)   # start: 0-based; end: 1-based
    curr['geneid'] = gene
    
    return curr

def _make_list_aligned_reads2(r_idx, read_dict, intersect_all, ex_coord):
    
    regx1 = re.compile("BC:Z:")
    regx2= re.compile("UB:Z:")
    
    exon_category_dict = {1: 'exon'}
    
    print(r_idx)
    curr = read_dict[r_idx].copy()
    
    intersect = sorted(set(intersect_all.loc[intersect_all['blockCount']==r_idx,'score'].values))
    n_intersect = len(intersect)

    curr['Exon_Index'] = ','.join(map(str,intersect))
    curr['Category'] = exon_category_dict.get(n_intersect, 'junction')
    
    bc = list(filter(regx1.match, curr['tags']))[0]
    ub = list(filter(regx2.match, curr['tags']))[0]
    
    curr['BC'] = bc.replace('BC:Z:','')
    curr['UB'] = ub.replace('UB:Z:','')
    curr['tags'] = ';'.join(curr['tags'])
    curr['exon_coordinates'] = ex_coord
    del curr['r_blocks']
    
    return list(pd.Series(curr).values)

def _make_list_aligned_reads_mp(r_idx):
    
    regx1 = re.compile("BC:Z:")
    regx2= re.compile("UB:Z:")
    
    exon_category_dict = {1: 'exon'}
    print(r_idx)
    curr = my_read_dict[r_idx].copy()
    
    intersect = sorted(set(my_intersect_all.loc[my_intersect_all['blockCount']==r_idx,'score'].values))
    n_intersect = len(intersect)
    
    curr['Exon_Index'] = ','.join(map(str,intersect))
    curr['Category'] = exon_category_dict.get(n_intersect, 'junction')
    
    bc = list(filter(regx1.match, curr['tags']))[0]
    ub = list(filter(regx2.match, curr['tags']))[0]
    
    curr['BC'] = bc.replace('BC:Z:','')
    curr['UB'] = ub.replace('UB:Z:','')
    curr['tags'] = ';'.join(curr['tags'])
    curr['exon_coordinates'] = my_ex_coord
    del curr['r_blocks']
    
    return list(pd.Series(curr).values)

def gtf2gene(gtf, outdir, field):
    
    filename = '%s/gene.gff' %(outdir)
    outF = open(filename, "w")
    
    if field == 'gene':   
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line): continue
                fds = line.strip().split('\t')
                if fds[2] == 'gene':
                    annot = fds[8].replace(' "', '=').replace('"; ',';').replace('";','').replace('; ',';').replace(' ','=')
                    outF.write('%s\t%s\n' %('\t'.join(fds[:8]), annot))
    else:
        gene_dict = defaultdict(dict)
        with open(gtf, 'r') as f:
            for line in f:
                fds = line.strip().split('\t')
                if fds[2] == 'transcript':
                    annot = fds[8].replace(' "', '=').replace('"; ',';').replace('";','').replace('; ',';').replace(' ','=')
                    curr_gene = annot.split(';')[0].split('=')[1]
                    if curr_gene not in gene_dict.keys():
                        gene_dict[curr_gene]['start'] = []
                        gene_dict[curr_gene]['end'] = []
                        gene_dict[curr_gene]['transcript'] = []
                    gene_dict[curr_gene]['chrom'] = fds[0]
                    gene_dict[curr_gene]['start'].append(int(fds[3]))
                    gene_dict[curr_gene]['end'].append(int(fds[4]))
                    gene_dict[curr_gene]['strand'] = fds[6]
                    gene_dict[curr_gene]['transcript'].append(annot.split(';')[1].split('=')[1])
                    
        for gene in gene_dict.keys():
            outF.write('%s\trefseq\tgene\t%s\t%s\t.\t%s\t.\tgene_id=%s;transcript_id=%s;gene_name=%s\n' %(gene_dict[gene]['chrom'], np.min(gene_dict[gene]['start']), np.max(gene_dict[gene]['end']), gene_dict[gene]['strand'], gene, ','.join(gene_dict[gene]['transcript']), gene))
    
    outF.close()
    return

def _fetch_exonic_reads(outdir, in_bam):
    
    exBed = '%s/exon_merged.bed' %(outdir)
    in_bam_filename = 'ex_%s' %(os.path.basename(in_bam))
    os.system('bedtools intersect -abam %s -b %s -wa -u > %s/%s' %(in_bam, exBed, outdir, in_bam_filename))
    os.system('samtools index %s/%s' %(outdir, in_bam_filename))
    
    return

def _get_reads(in_bam_uniq, in_bam_multi, outdir, chrom, nproc, n_read_limit, passed_cells, mRds, gene):
    
    gobj = geneObj(in_bam_uniq, in_bam_multi, outdir)
    gobj.get_exon_coordinates(gene)
    
    if not os.path.exists('%s/keptReads/%s/%s_aligned_reads.csv' %(outdir, chrom, gobj.gene)): 
        os.system('echo "Start gene %s..." >> %s/keptReads/%s/_log' %(gobj.gene, outdir, chrom))
    else:
        os.system('echo "Gene %s exists in output directory...Skip..." >> %s/keptReads/%s/_log' %(gobj.gene, outdir, chrom))
        return
    
    report_gene = None
    
    if nproc < 2:
        report_gene = gobj.get_aligned_reads(n_read_limit, passed_cells)
    else:
        gobj = get_aligned_reads_mp(gobj, nproc, passed_cells)
        
    if report_gene is not None: return report_gene
    if gobj.uniq_aligned_reads is None: return None
    
    if mRds:
        if nproc < 2:
            gobj.get_aligned_reads_from_multi(passed_cells)
        else:
            gobj = get_aligned_reads_from_multi_mp(gobj, nproc, passed_cells)
    
        if gobj.multi_aligned_reads is not None:
            aligned = pd.concat([gobj.uniq_aligned_reads, gobj.multi_aligned_reads], axis=0)
        else:
            aligned = gobj.uniq_aligned_reads
    else:
        aligned = gobj.uniq_aligned_reads
       
    os.system('echo "%s has %s aligned reads..." >> %s/keptReads/%s/_log' %(gobj.gene, aligned.shape[0], outdir, chrom))
    
    p = subprocess.Popen('rm %s/.tempDir/_%s*' %(outdir, gobj.gene), shell=True)
    (output, err) = p.communicate()  
    
    if aligned.shape[0] > 0:
        aligned.to_csv('%s/keptReads/%s/%s_aligned_reads.csv' %(outdir, chrom, gobj.gene), sep="\t", index=False, header=False)
    return None


def fetch_gene_reads(in_bam_uniq, in_bam_multi, conf, species, outdir, spikein=False):
    
    gtf = conf['annotation']['%s_%s_gtf' %(species,conf['annotation']['gtf_source'])]
    if conf['annotation']['gtf_source'] == 'refseq':
        field = 'transcript'
    else:
        field = 'gene'
    gtf2exon(gtf, outdir, spikein)
    gtf2gene(gtf, outdir, field)
    
    pool = mp.Pool(2)
    func = partial(_fetch_exonic_reads, outdir)
    pool.map(func, [in_bam_uniq, in_bam_multi], chunksize=1)
    in_bam_uniq = '%s/ex_%s' %(outdir, os.path.basename(in_bam_uniq))
    in_bam_multi = '%s/ex_%s' %(outdir, os.path.basename(in_bam_multi))
    
    cells_to_use = list(pd.read_table(conf['annotation']['zumi_keptbarcode'], header=None, index_col=None, sep=",").iloc[:,0].values)
    
    cmd = 'cut -f1 %s/gene.gff | sort | uniq' %(outdir)
    p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    rcds = p.communicate()[0].decode("utf-8")
    chrom_list = [item for item in rcds.strip().split('\n')]
    
    if not os.path.exists('%s/.tempDir' %outdir): os.makedirs('%s/.tempDir' %outdir)
    if not os.path.exists('%s/keptReads' %outdir): os.makedirs('%s/keptReads' %outdir)
    
    for chrom in chrom_list:
        
        print('...for genes on %s' %(chrom))
        if not os.path.exists('%s/keptReads/%s' %(outdir,chrom)): os.makedirs('%s/keptReads/%s' %(outdir,chrom))
        os.system('> %s/keptReads/%s/_log' %(outdir,chrom))
        
        os.system('echo "*** genes on %s ***" >> %s/keptReads/%s/_log' %(chrom, outdir, chrom))
        cmd = 'grep ^%s[[:space:]] %s/gene.gff' %(chrom, outdir)
        p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        rcds = p.communicate()[0].decode("utf-8")
        
        genes = [item for item in rcds.strip().split('\n')]
        
        pool = mp.Pool(processes=int(conf['expression']['nproc']))
        func = partial(_get_reads, in_bam_uniq, in_bam_multi, outdir, chrom, 1, int(conf['expression']['n_read_limit']), cells_to_use, False) 
        report_genes = pool.map(func, genes, chunksize=1)
        pool.close()
            
        report_genes = list(filter(None, report_genes))
        for gname in report_genes:
            print(gname)
            gene = [gg for gg in genes if re.search(gname, gg)][0]
            results = _get_reads(in_bam_uniq, in_bam_multi, outdir, chrom, int(conf['expression']['nproc']), int(conf['expression']['n_read_limit']), cells_to_use, False, gene) 
        
        p = subprocess.Popen('rm -rf %s/.tmp' %(outdir), shell=True)
        (output, err) = p.communicate()
        
        if not os.path.exists('%s/.tmp' %outdir): os.makedirs('%s/.tmp' %outdir)
        
    return