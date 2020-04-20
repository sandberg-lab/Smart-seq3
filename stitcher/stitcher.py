#!/usr/bin/env python
# V 1.0
# anton jm larsson anton.larsson@ki.se
import argparse
import pysam
import pandas as pd
import numpy as np
import pygtrie
import portion as P
import itertools
import sys
import time
import os
from joblib import delayed,Parallel
from multiprocessing import Process, JoinableQueue
__version__ = '1.0'
nucleotides = ['A', 'T', 'C', 'G']
nuc_dict = {'A':0, 'T':1, 'C':2, 'G':3, 'N': 4}
np.seterr(divide='ignore')

import rpy2.robjects as robjects
from rpy2.robjects.packages import STAP

def make_ll_array(e):
    y = np.array([e[0]/3,e[0]/3,e[0]/3,e[0]/3])
    if e[1] != 4:
        y[e[1]] = 1-e[0]
    return np.log10(y)


def intervals_extract(iterable): 
    iterable = sorted(set(iterable)) 
    for key, group in itertools.groupby(enumerate(iterable), 
    lambda t: t[1] - t[0]): 
        group = list(group) 
        yield [group[0][1], group[-1][1]] 

def interval(t):
    return P.from_data([(True,i[0],i[1], True) for i in t])

def get_time_formatted(time):
    day = time // (24 * 3600)
    time = time % (24 * 3600)
    hour = time // 3600
    time %= 3600
    minutes = time // 60
    time %= 60
    seconds = time
    s = ''.join(['{} day{}, '.format(day, 's'*(1 != day))*(0 != day), 
                 '{} hour{}, '.format(hour,'s'*(1 != hour))*(0 != hour), 
                 '{} minute{}, '.format(minutes,'s'*(1 != minutes))*(0 != minutes), 
                 '{:.2f} second{}, '.format(seconds,'s'*(1 != seconds))*(0 != seconds)])
    s = s[:-2]
    s = s + '.'
    return s

def get_insertions_locs(cigtuples):
    insertion_locs = []
    l = 0
    for c in cigtuples:
        if c[0] == 0:
            l += c[1]
        elif c[0] == 1:
            for i in range(c[1]):
                insertion_locs.append(l)
                l += 1
    return insertion_locs

def get_skipped_tuples(cigtuples, ref_positions):
    skipped_locs = []
    l = 0
    for c in cigtuples:
        if c[0] == 0:
            l += c[1]
        elif c[0] == 3:
            skipped_locs.append((ref_positions[l-1]+1, ref_positions[l]-1))
    return skipped_locs

def stitch_reads(read_d, mol_dict=None, cell = None, gene = None, umi = None):
    master_read = {}
    seq_df = None
    qual_df = None
    nreads = len(read_d)
    reverse_read1 = []
    read_ends = [0]*nreads
    read_starts = [0]*nreads
    exonic_list = [0]*nreads
    intronic_list = [0]*nreads
    seq_series_list = []
    qual_series_list = []
    for i,read in enumerate(read_d):
        if read.has_tag('GE'):
            exonic = True
        else:
            exonic = False
        if read.has_tag('GI'):
            intronic = True
        else:
            intronic = False
        p_x = list(10**(-np.float_(np.array(read.query_alignment_qualities))/10))
        seq = [nuc_dict[c] for c in read.query_alignment_sequence]
        cigtuples = read.cigartuples
        insertion_locs = get_insertions_locs(cigtuples)
        try:
            for loc in insertion_locs:
                    del seq[loc]
                    del p_x[loc]
        except IndexError:
            return (False, read.query_name + '\n')
            continue
        ref_positions = read.get_reference_positions()
        skipped_intervals = get_skipped_tuples(cigtuples, ref_positions)

        if mol_dict is None:
            if read.is_read1:
                reverse_read1.append(read.is_reverse)
            else:
                read_starts[i] = read.reference_start
                read_ends[i] = read.reference_end
        else:
            if mol_dict[cell][gene][umi].is_reverse:
                if not read.is_reverse:
                    read_starts[i] = read.reference_start
            else:
                if read.is_reverse:
                    read_ends[i] = read.reference_end
        exonic_list[i] = exonic
        intronic_list[i] = intronic
        seq_series = pd.Series(seq, index=ref_positions)
        seq_series.name = read.query_name
        seq_series_list.append(seq_series)
        qual_series = pd.Series(p_x, index=ref_positions)
        qual_series.name = read.query_name
        qual_series_list.append(qual_series)
        if len(master_read) == 0:
            master_read['skipped_intervals'] = skipped_intervals
        else:
            master_read['skipped_intervals'].extend(skipped_intervals)
    master_read['SN'] = read.reference_name
    qual_df = pd.DataFrame(qual_series_list).fillna(3).transpose()
    seq_df = pd.DataFrame(seq_series_list).fillna(4).astype(int).transpose()
    merged_df = pd.DataFrame(np.rec.fromarrays((qual_df.values, seq_df.values)).tolist(),
                      columns=qual_df.columns,
                      index=qual_df.index)

    qual_probs = 10**(merged_df.applymap(make_ll_array).apply(lambda x: np.concatenate(list(x))).sum(axis=1)).values.reshape(qual_df.shape[0], 4)

    normed_probs = qual_probs/qual_probs.sum(axis=1)[:, np.newaxis]
    prob_max = np.max(normed_probs, axis=1)
    master_read['seq'] = ''.join([nucleotides[x] if p > 0.3 else 'N' for p, x in zip(prob_max, np.argmax(normed_probs, axis=1))])
    master_read['phred'] = np.nan_to_num(np.rint(-10*np.log10(1-prob_max+1e-13)))
    if mol_dict is None:
        v, c = np.unique(reverse_read1, return_counts=True)
        m = c.argmax()
        master_read['is_reverse'] = v[m]
    else:
        master_read['is_reverse'] = mol_dict[cell][gene][umi].is_reverse
    if master_read['is_reverse']:
        master_read['ends'] = list(set(read_starts)-{0})
    else:
        master_read['ends'] = list(set(read_ends)-{0})
    master_read['ref_intervals'] = interval(intervals_extract(seq_df.index.values))
    master_read['skipped_intervals'] = interval(list(set(master_read['skipped_intervals'])))
    master_read['del_intervals'] =  ~(master_read['ref_intervals'] | master_read['skipped_intervals'])
    master_read['NR'] = nreads
    master_read['IR'] = np.sum(intronic_list)
    master_read['ER'] = np.sum(exonic_list)
    master_read['cell'] = cell
    master_read['gene'] = gene
    master_read['umi'] = umi
    return (True, convert_to_sam(master_read))


def assemble_reads(bamfile,gene_to_stitch, cell_set, q):
    readtrie = pygtrie.StringTrie()
    bam = pysam.AlignmentFile(bamfile, 'rb')
    gene_of_interest = gene_to_stitch['gene_id']
    for read in bam.fetch(gene_to_stitch['seqid'], gene_to_stitch['start'], gene_to_stitch['end']):
        cell = read.get_tag('BC')
        if cell_set is not None:
            if cell not in cell_set:
                continue
        umi = read.get_tag('UB')
        if umi == '':
            continue
        if read.has_tag('GE'):
            gene_exon = read.get_tag('GE')
        else:
            gene_exon = 'Unassigned'
        if read.has_tag('GI'):
            gene_intron = read.get_tag('GI')
        else:
            gene_intron = 'Unassigned'

        # if it maps to the intron or exon of a gene
        if gene_intron != 'Unassigned' or gene_exon != 'Unassigned':
            # if it is a junction read
            if gene_intron == gene_exon:
                gene = gene_intron
            # if it's an only intronic read
            elif gene_intron != 'Unassigned' and gene_exon == 'Unassigned':
                gene = gene_intron
            # if it's an only exonic read
            elif gene_exon != 'Unassigned' and gene_intron == 'Unassigned':
                gene = gene_exon
            # if the exon and intron gene tag contradict each other
            else:
                continue
        else:
            continue
        if read.is_paired and not read.is_unmapped and not read.mate_is_unmapped and gene == gene_of_interest \
        and read.is_proper_pair:
            node = '{}/{}/{}'.format(cell,gene,umi)
            if readtrie.has_node(node):
                readtrie[node].append(read)
            else:
                readtrie[node] = [read]
    mol_list = []
    mol_append = mol_list.append
    for node, mol in readtrie.iteritems():
        info = node.split('/')
        read_names = [r.query_name for r in mol]
        if 2*len(set(read_names)) == len(mol):
            mol_append(stitch_reads(mol, None, info[0], info[1], info[2]))
        else:
            mol_append((False, '{} does not have all reads within the annotated gene\n'.format(node)))
    del readtrie
    q.put(True, mol_list)
    return gene_of_interest


def make_POS_and_CIGAR(stitched_m):
    CIGAR = ''
    conflict = False
    interval_list = []
    ref_and_skip_intersect = stitched_m['ref_intervals'] & stitched_m['skipped_intervals']
    nreads_conflict = 0
    if not ref_and_skip_intersect.empty:
        conflict = True
        nreads_conflict = len(list(P.iterate(ref_and_skip_intersect, step=1))) 
        stitched_m['skipped_intervals'] = stitched_m['skipped_intervals'] - ref_and_skip_intersect
        interval_list = [i for t in P.to_data(ref_and_skip_intersect) for i in t[1:-1]]
    ref_tuples = [(i[1] if i[0] else i[1]+1, i[2] if i[3] else i[2]-1) for i in P.to_data(stitched_m['ref_intervals'])]
    if stitched_m['skipped_intervals'].empty:
        skipped_tuples = []
    else:
        skipped_tuples = [(i[1] if i[0] else i[1]+1, i[2] if i[3] else i[2]-1) for i in P.to_data(stitched_m['skipped_intervals'])]
    if stitched_m['del_intervals'].empty:
        del_tuples = []
    else:
        del_tuples = [(i[1] if i[0] else i[1]+1, i[2] if i[3] else i[2]-1) for i in P.to_data(stitched_m['del_intervals'])[1:-1]]
    POS = ref_tuples[0][0] + 1
    tuple_dict = {'M': ref_tuples, 'N': skipped_tuples, 'D': del_tuples}
    l = []
    while sum(len(t) for t in tuple_dict.values()) > 0:
        pos_dict = {k:v[0][0] for k,v in tuple_dict.items() if len(v) > 0}
        c = min(pos_dict, key=pos_dict.get)
        n_bases = np.int_(tuple_dict[c[0]][0][1]-tuple_dict[c[0]][0][0])+1
        CIGAR += '{}{}'.format(n_bases,c[0])
        if c[0] == 'M':
            l.append(n_bases)
        del tuple_dict[c[0]][0]
    return POS, CIGAR, conflict, nreads_conflict, interval_list

def convert_to_sam(stitched_m):
    sam_dict = {}
    POS, CIGAR, conflict, nreads_conflict, interval_list = make_POS_and_CIGAR(stitched_m)
    sam_dict['QNAME'] = '{}:{}:{}'.format(stitched_m['cell'],stitched_m['gene'],stitched_m['umi'])
    sam_dict['FLAG'] = str(16*stitched_m['is_reverse'])
    sam_dict['RNAME'] = stitched_m['SN']
    sam_dict['POS'] = str(POS)
    sam_dict['MAPQ'] = str(255)
    sam_dict['CIGAR'] = CIGAR
    sam_dict['RNEXT'] = '*'
    sam_dict['PNEXT'] = str(0)
    sam_dict['TLEN'] = str(0)
    sam_dict['SEQ'] = stitched_m['seq']
    sam_dict['QUAL'] = "".join([chr(int(p)) for p in np.clip(stitched_m['phred'],0,126-33)+33])
    sam_dict['NR'] = 'NR:i:{}'.format(stitched_m['NR'])
    sam_dict['ER'] = 'ER:i:{}'.format(stitched_m['ER'])
    sam_dict['IR'] = 'IR:i:{}'.format(stitched_m['IR'])
    sam_dict['BC'] = 'BC:Z:{}'.format(stitched_m['cell'])
    sam_dict['XT'] = 'XT:Z:{}'.format(stitched_m['gene'])
    sam_dict['UB'] = 'UB:Z:{}'.format(stitched_m['umi'])
    sam_dict['EL'] = 'EL:B:I,{}'.format(','.join([str(e) for e in stitched_m['ends']]))
    if conflict:
        sam_dict['NC'] = 'NC:i:{}'.format(nreads_conflict)
        sam_dict['IL'] = 'IL:B:I,{}'.format(','.join([str(e) for e in interval_list]))
    return '\t'.join(list(sam_dict.values())) + '\n'

def yield_reads(read_dict):
    for cell in read_dict:
        for gene in read_dict[cell]:
            #print('\t', gene)
            for umi in read_dict[cell][gene]:
                #print('\t\t', umi)
                yield read_dict[cell][gene][umi], None, cell, gene, umi


def create_write_function(filename, bamfile, version):
    bam = pysam.AlignmentFile(bamfile, 'rb')
    header = bam.header['SQ']
    def write_sam_file(q):
        error_file = open('{}_error.log'.format(os.path.splitext(filename)[0]), 'w')
        with open(filename, 'w') as samfile:
            # write header
            samfile.write('@HD\tVN:1.0\tSO:unknown\n')
            for SQ in header:
                samfile.write('@SQ\tSN:{}\tLN:{}\n'.format(SQ['SN'],SQ['LN']))
            samfile.write('@PG\tID:stitcher\tVN:{}\n'.format(version))
            while True:
                good, mol_list = q.get()
                if good is None: break
                if good:
                    for success, mol in mol_list:
                        if success:
                            samfile.write(mol)
                        else:
                            error_file.write(mol)
                q.task_done()
            samfile.truncate()
            q.task_done()
        error_file.close()
        return None
    return write_sam_file

def extract(d, keys):
    return dict((k, d[k]) for k in keys if k in d)

def construct_stitched_molecules(infile, outfile,gtffile,counts, cells, contig, threads, version):


    if cells is not None:
        cell_set = set([line.rstrip() for line in open(cells)])
    else:
        cell_set = None

    gene_list = []
    with open(gtffile, 'r') as f:
        for line in f:
            l = line.split('\t')
            if len(l) < 8:
                continue
            if l[2] == 'gene':
                if contig is not None:
                    if l[0] == contig:
                        gene_list.append({'gene_id': l[8].split(' ')[1].replace('"', '').strip(';'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4])})
                    else:
                        continue
                else:
                    gene_list.append({'gene_id': l[8].split(' ')[1].replace('"', '').strip(';'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4])})
    gene_df = pd.DataFrame(gene_list)
    gene_df.index = gene_df['gene_id']
    if counts is not None:
        mfunc = 'to_df <- function(dobj){return(as.data.frame(as.matrix(dobj)))}'
        rsparse2pandas = STAP(mfunc, "to_df")
        readRDS = robjects.r['readRDS']
        zumis_data = readRDS(counts)
        zd = dict(zip(zumis_data.names, list(zumis_data)))

        for key in ['readcount']:#, 'readcount']:
            zd[key] = dict(zip(zd[key].names, list(zd[key])))
            zd[key]['inex'] = dict(zip(zd[key]['inex'].names, list(zd[key]['inex'])))
            zd[key]['inex']['all'] = rsparse2pandas.to_df(zd[key]['inex']['all'])

        total_counts = zd['readcount']['inex']['all'].sum(axis=1).sort_values(ascending=False)
        total_counts.name = 'total_counts'
        
        gene_df = gene_df.join(total_counts).fillna(0).sort_values('total_counts', ascending=False)
        
    params = Parallel(n_jobs=threads, verbose = 3, backend='multiprocessing')(delayed(assemble_reads)(infile, gene, cell_set) for g,gene in gene_df.iterrows())


    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Stitch together molecules from reads sharing the same UMI')
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-o','--output', metavar='output', type=str, help='Output .sam file')
    parser.add_argument('-g','--gtf', metavar='gtf', type=str, help='gtf file with gene information')
    parser.add_argument('--counts', metavar='counts', type=str, help='zUMIs .rds file with read counts')
    parser.add_argument('-t', '--threads', metavar='threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--cells', default=None, metavar='cells', type=str, help='List of cell barcodes to stitch molecules')
    parser.add_argument('--contig', default=None, metavar='contig', type=str, help='Restrict stitching to contig')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    args = parser.parse_args()
    infile = args.input
    if infile is None:
        raise Exception('No input file provided.')
    outfile = args.output
    if outfile is None:
        raise Exception('No output file provided.')
    gtffile = args.gtf  
    if gtffile is None:
        raise Exception('No gtf file provided.')
    counts = args.counts
    threads = int(args.threads)
    cells = args.cells
    contig = args.contig
    m = Manager()
    q = m.JoinableQueue()
    p = Process(target=create_write_function(filename=outfile, bamfile=infile, version=__version__), args=(q,))
    p.start()
    
    print('Stitching reads for {}'.format(infile))
    
    start = time.time()
    construct_stitched_molecules(infile, outfile, gtffile,counts, cells, contig, threads,q, __version__)
    q.put((None,None))
    p.join()
    end = time.time()
    
    print('Finished writing stitched molecules from {} to {}, took {}'.format(infile, outfile, get_time_formatted(end-start)))
