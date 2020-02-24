# V 0.4
# anton jm larsson anton.larsson@ki.se
import argparse
import pysam
import pandas as pd
import numpy as np
from interval import interval
import sys
import time
from joblib import delayed,Parallel

nucleotides = ['A', 'T', 'C', 'G']
def make_ll_array(e):
    p=e[0]
    i=e[1]
    y = np.repeat(p/3, 4)
    if i != 4:
        y[i] = 1-p
    return np.log10(y)

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

def get_skipped_intervals(cigtuples, ref_positions):
    skipped_locs = interval()
    l = 0
    for c in cigtuples:
        if c[0] == 0:
            l += c[1]
        elif c[0] == 3:
            skipped_locs = skipped_locs | interval[ref_positions[l-1]+1, ref_positions[l]-1]
    return skipped_locs

def get_ref_intervals(ref_positions):
    ref_intervals = interval[ref_positions[0], ref_positions[1]] 
    for i in range(len(ref_positions)-2):
        if (ref_positions[i+2] - ref_positions[i+1]) == 1:
            ref_intervals = ref_intervals | interval[ref_positions[i+1],ref_positions[i+2]]
    return ref_intervals

def get_del_intervals(ref_skip_union):
    prev_interval = None
    del_intervals = interval()
    for x in ref_skip_union:
        if prev_interval is not None:
            if prev_interval[1] + 1 == x[0]:
                prev_interval = x
                continue
            else:
                del_intervals = del_intervals | interval[prev_interval[1]+1, x[0]-1]
            prev_interval = x
        else:
            prev_interval = x
    return del_intervals

def get_read_info(read):
    if read.has_tag('GE'):
        exonic = True
    else:
        exonic = False
    if read.has_tag('GI'):
        intronic = True
    else:
        intronic = False
    p_x = list(10**(-np.float_(np.array(read.query_alignment_qualities))/10))
    seq = [char for char in read.query_alignment_sequence]
    cigtuples = read.cigartuples
    insertion_locs = get_insertions_locs(cigtuples)
    for loc in insertion_locs:
            del seq[loc]
            del p_x[loc]
    ref_positions = read.get_reference_positions()
    try:
        skipped_intervals = get_skipped_intervals(cigtuples, ref_positions)
    except IndexError:
        print(read)
    new_read = {'name': read.query_name, 'p_x': p_x, 'seq': seq, 'ref_positions': ref_positions,'skipped_intervals':skipped_intervals,
               'intronic': intronic, 'exonic': exonic,
                'is_reverse': read.is_reverse, 'reference_name': read.reference_name, 'reference_start': read.reference_start, 
               'reference_end': read.reference_end, 'is_read1': read.is_read1}
    return new_read

def stitch_reads(read_d, mol_dict=None, cell = None, gene = None, umi = None):
    master_read = {}
    seq_df = None
    qual_df = None
    reverse_read1 = []
    read_ends = []
    read_starts = []
    exonic_list = []
    intronic_list = []
    for read in read_d:

        if mol_dict is None:
            if read['is_read1']:
                reverse_read1.append(read['is_reverse'])
            else:
                read_starts.append(read['reference_start'])
                read_ends.append(read['reference_end'])
        else:
            if mol_dict[cell][gene][umi]['is_reverse']:
                if not read['is_reverse']:
                    read_starts.append(read['reference_start'])
            else:
                if read['is_reverse']:
                    read_ends.append(read['reference_end'])
        exonic_list.append(read['exonic'])
        intronic_list.append(read['intronic'])
        seq_series = pd.Series(read['seq'], index=read['ref_positions'])
        seq_series.name = read['name']
        qual_series = pd.Series(read['p_x'], index=read['ref_positions'])
        qual_series.name = read['name']
        if seq_df is None:
            seq_df = pd.DataFrame(seq_series)
            qual_df = pd.DataFrame(qual_series)
        else:
            seq_df = seq_df.join(seq_series,how='outer', rsuffix = '_right')
            qual_df = qual_df.join(qual_series,how='outer', rsuffix = '_right')
        if len(master_read) == 0:
            master_read['skipped_intervals'] = read['skipped_intervals']
        else:
            master_read['skipped_intervals'] = master_read['skipped_intervals'] | read['skipped_intervals']
    master_read['SN'] = read['reference_name']
    qual_df = qual_df.replace(np.nan, 3)
    seq_df = seq_df.replace(['A','T', 'C', 'G', np.nan, 'N'],[0,1,2,3,4,4])
    merged_df = pd.DataFrame(np.rec.fromarrays((qual_df.values, seq_df.values)).tolist(), 
                      columns=qual_df.columns,
                      index=qual_df.index)
    qual_probs = 10**(merged_df.applymap(make_ll_array).apply(lambda x: np.concatenate(list(x))).sum(axis=1)).values.reshape(qual_df.shape[0], 4)
    normed_probs = qual_probs/qual_probs.sum(axis=1)[:, np.newaxis]
    prob_max = np.max(normed_probs, axis=1)
    master_read['seq'] = ''.join([nucleotides[x] if p > 0.3 else 'N' for p, x in zip(prob_max, np.argmax(normed_probs, axis=1))])
    master_read['phred'] = np.rint(-10*np.log10(1-prob_max+1e-13))
    if mol_dict is None:
        v, c = np.unique(reverse_read1, return_counts=True)
        m = c.argmax()
        master_read['is_reverse'] = v[m]
    else:
        master_read['is_reverse'] = mol_dict[cell][gene][umi].is_reverse
    if master_read['is_reverse']:
        master_read['ends'] = list(set(read_starts))
    else:
        master_read['ends'] = list(set(read_ends))
    master_read['ref_intervals'] = get_ref_intervals(seq_df.index)
    ref_skip_union = (master_read['ref_intervals'] | master_read['skipped_intervals'])
    master_read['del_intervals'] =  get_del_intervals(ref_skip_union)
    master_read['NR'] = len(read_d)
    master_read['intronic_reads'] = sum(intronic_list)
    master_read['exonic_reads'] = sum(exonic_list)
    master_read['cell'] = cell
    master_read['gene'] = gene
    master_read['umi'] = umi
    return master_read

def make_read_dict(bamfile,contig, read_dict = {}, cell_list = None):
    for read in bamfile.fetch(contig):
        if read.is_paired and not read.is_unmapped and not read.mate_is_unmapped:
            cell = read.get_tag('BC')
            if cell_list is not None:
                 if cell not in cell_list:
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
                elif gene_intron != 'Unassigned' and gene_intron != gene_exon:
                    gene = gene_intron
                # if it's an only exonic read
                elif gene_exon != 'Unassigned' and gene_intron != gene_exon:
                    gene = gene_exon
                # if the exon and intron gene tag contradict each other
                else:
                    print('contradiction')
                    continue
            else:
                continue

            umi = read.get_tag('UB')
            
            if cell in read_dict.keys():
                if gene in read_dict[cell].keys():
                    if umi in read_dict[cell][gene].keys():
                        read_dict[cell][gene][umi].append(get_read_info(read))
                    else:
                        
                        read_dict[cell][gene][umi] = [get_read_info(read)]
                else:
                    read_dict[cell][gene] = {}
                    read_dict[cell][gene][umi] = [get_read_info(read)]
            else:
                read_dict[cell] = {}
                read_dict[cell][gene] = {}
                read_dict[cell][gene][umi] = [get_read_info(read)]
            
    return read_dict


def make_POS_and_CIGAR(stitched_m):
    CIGAR = ''
    ref_tuples = [(int(i.inf),int(i.sup)) for i in stitched_m['ref_intervals']]
    skipped_tuples = [(int(i.inf),int(i.sup)) for i in stitched_m['skipped_intervals']]
    del_tuples = [(int(i.inf),int(i.sup)) for i in stitched_m['del_intervals']]
    POS = ref_tuples[0][0] + 1
    tuple_dict = {'M': ref_tuples, 'N': skipped_tuples, 'D': del_tuples}
    conflict = len(stitched_m['ref_intervals'] & stitched_m['skipped_intervals']) > 0
    if conflict:
        return POS, '*', conflict
    l = []
    while sum(len(t) for t in tuple_dict.values()) > 0:
        pos_dict = {k:v[0][0] for k,v in tuple_dict.items() if len(v) > 0}
        c = min(pos_dict, key=pos_dict.get)
        n_bases = np.int_(tuple_dict[c[0]][0][1]-tuple_dict[c[0]][0][0])+1
        CIGAR += '{}{}'.format(n_bases,c[0])
        if c[0] == 'M':
            l.append(n_bases)
        del tuple_dict[c[0]][0]
    conflict = sum(l) != len(stitched_m['seq'])
    if conflict:
        CIGAR = '*'
    return POS, CIGAR, sum(l) != len(stitched_m['seq'])

def convert_to_sam(stitched_m):
    sam_dict = {}
    POS, CIGAR, conflict = make_POS_and_CIGAR(stitched_m)
    sam_dict['QNAME'] = '{}:{}:{}'.format(stitched_m['cell'],stitched_m['gene'],stitched_m['umi'])
    sam_dict['FLAG'] = str(16*stitched_m['is_reverse']+4*conflict)
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
    return '\t'.join(list(sam_dict.values())) + '\n'

def yield_reads(read_dict):
    for cell in read_dict:
        for gene in read_dict[cell]:
            #print('\t', gene)
            for umi in read_dict[cell][gene]:
                #print('\t\t', umi)
                yield read_dict[cell][gene][umi], None, cell, gene, umi

def write_sam_file(stitched_mols, filename, bamfile):
    with open(filename, 'w') as samfile:
        # write header
        samfile.write('@HD\tVN:1.0\tSO:unknown\n')
        for SQ in bamfile.header['SQ']:
            samfile.write('@SQ\tSN:{}\tLN:{}\n'.format(SQ['SN'],SQ['LN']))
        samfile.write('@PG\tID:stitcher\tVN:0.1\n')
        for mol in stitched_mols:
            samfile.write(convert_to_sam(mol))
        samfile.truncate()
    return None

def extract(d, keys):
    return dict((k, d[k]) for k in keys if k in d)

def construct_stitched_molecules(infile, outfile, cells, contig, threads):
    print('Gathering reads for {}'.format(infile))
    start = time.time()
    bamfile = pysam.AlignmentFile(infile, 'rb')
    if cells is not None:
        cell_list = [line.rstrip() for line in open(cells)]
    else:
        cell_list = None
    read_dict = make_read_dict(bamfile, contig, read_dict={}, cell_list = cell_list)
    end = time.time()
    print('Finished gathering reads for {}, took {}'.format(infile, get_time_formatted(end-start)))

    print('Stitching reads into molecules for {}'.format(infile))
    start = time.time()
    stitched_mols = Parallel(n_jobs=threads, verbose = 3)(delayed(stitch_reads)(*d) for d in yield_reads(read_dict))
    end = time.time()
    print('Finished stitching reads into molecules for {}, took {}'.format(infile, get_time_formatted(end-start)))

    print('Writing stitched molecules from {} to {}'.format(infile, outfile))
    start = time.time()
    write_sam_file(stitched_mols, outfile, bamfile)
    end = time.time()
    print('Finished writing stitched molecules from {} to {}, took {}'.format(infile, outfile, get_time_formatted(end-start)))
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Stitch together molecules from reads sharing the same UMI')
    parser.add_argument('--i',metavar='input', type=str, nargs=1, help='Input .bam file')
    parser.add_argument('--o', metavar='output', type=str, nargs=1, help='Output .sam file')
    parser.add_argument('--t', metavar='threads', type=int, nargs=1, default=1, help='Number of threads')
    parser.add_argument('--cells', default=None, metavar='cells', type=str, nargs=1, help='List of cell barcodes to stitch molecules')
    parser.add_argument('--contig', default=None, metavar='contig', type=str, nargs=1, help='Restrict stitching to contig')
    args = parser.parse_args()
    infile = args.i[0]
    outfile = args.o[0]
    threads = args.t
    if args.cells is None:
        cells = args.cells
    else:
        cells = args.cells[0]
    if args.contig is None:
        contig = args.contig
    else:
        contig = args.contig[0]
    construct_stitched_molecules(infile, outfile, cells, contig, threads)
