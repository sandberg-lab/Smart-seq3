import argparse
import pysam
import pandas as pd
import numpy as np
from interval import interval
from stitcher import *
import sys

def make_assembled_dict(bamfile,contig, read_dict = {}):
    #print(read_dict)
    for read in bamfile.fetch(contig):
        if not read.is_unmapped and read.has_tag('XT'):
            cell = read.get_tag('BC')
            gene = read.get_tag('XT')
            umi = read.get_tag('UB')
            #print(cell,gene,umi)
            if cell in read_dict.keys():
                if gene in read_dict[cell].keys():
                    if umi in read_dict[cell][gene].keys():
                        #print(cell,gene,umi)
                        read_dict[cell][gene][umi].append(read)
                    else:
                        read_dict[cell][gene][umi] = read
                else:
                    read_dict[cell][gene] = {}
                    read_dict[cell][gene][umi] = read
            else:
                read_dict[cell] = {}
                read_dict[cell][gene] = {}
                read_dict[cell][gene][umi] = read
    return read_dict

def make_combined_dict(bamfile,contig, read_dict = {}):
    for read in bamfile.fetch(contig):
        if not read.is_unmapped and read.has_tag('XT'):
            cell = read.get_tag('BC')
            gene = read.get_tag('XT')
            umi = read.get_tag('UB')
            #print(cell,gene,umi)
            if cell in read_dict.keys():
                if gene in read_dict[cell].keys():
                    if umi in read_dict[cell][gene].keys():
                        #print(cell,gene,umi)
                        read_dict[cell][gene][umi].append(read)
                    else:
                        read_dict[cell][gene][umi] = [read]
                else:
                    read_dict[cell][gene] = {}
                    read_dict[cell][gene][umi] = [read]
            else:
                read_dict[cell] = {}
                read_dict[cell][gene] = {}
                read_dict[cell][gene][umi] = [read]
    return read_dict

def find_internal_reads(internal_reads, assembled_read_dict):
    internal_read_names_dict = {}
    for cell in assembled_read_dict.keys():
        internal_read_names_dict[cell] = {}
        for gene in assembled_read_dict[cell].keys():
            #print(gene)
            internal_read_names_dict[cell][gene] = {}
            for umi,mol in assembled_read_dict[cell][gene].items():
                internal_read_names_dict[cell][gene][umi] = []
                if not mol.is_unmapped and gene in internal_reads[cell]:
                    #print('\t', umi,mol.is_reverse,mol.get_tag('NR'),mol.reference_start, mol.reference_end,  mol.reference_end- mol.reference_start)
                    if '' in internal_reads[cell][gene]:
                        for read in internal_reads[cell][gene]['']:
                            if mol.is_reverse:
                                diff = read.reference_end - np.array(mol.get_tag('EL'))
                            else:
                                diff = np.array(mol.get_tag('EL')) - read.reference_start
                            if any(diff == 9):
                                if read.is_paired and not read.mate_is_unmapped:
                                    internal_read_names_dict[cell][gene][umi].append(read.query_name)
                                    #print('\t\t',read.query_name,read.is_read1, read.reference_start, read.next_reference_start, mol.reference_end - read.reference_start)
                    else:
                        continue
    return internal_read_names_dict

def in_other_mol(d, u, readname):
    for k in d.keys():
        if u != k:
            if readname in d[k]:
                return True
    return False

def to_delete(combined_read_dict):
    tuplelist = []
    for cell in combined_read_dict.keys():
        for gene in combined_read_dict[cell].keys():
            for umi in combined_read_dict[cell][gene].keys():
                if len(combined_read_dict[cell][gene][umi]) == 1:
                    tuplelist.append((cell,gene,umi))
    return tuplelist

def make_internal_read_dict(internal_reads, internal_read_names_dict):
    internal_read_dict = {}
    for cell in internal_read_names_dict:
        internal_read_dict[cell] = {}
        for gene in internal_read_names_dict[cell].keys():
            internal_read_dict[cell][gene] = {}
            for umi in internal_read_names_dict[cell][gene].keys():
                if gene in internal_reads[cell]:
                    if '' in internal_reads[cell][gene]:
                        for read in internal_reads[cell][gene]['']:
                            if read.query_name in internal_read_names_dict[cell][gene][umi] and not in_other_mol(internal_read_names_dict[cell][gene], umi, read.query_name):
                                if umi in internal_read_dict[cell][gene].keys():
                                    internal_read_dict[cell][gene][umi].append(read)
                                else:
                                    internal_read_dict[cell][gene][umi]= [read]
    return internal_read_dict

def stitch_internal(read_d):
    master_read = None
    for read in read_d:
        p_x = list(10**-(np.float_(np.array(read.query_alignment_qualities))/10))
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
            #print(read)
            continue
        new_read = {'p_x': p_x, 'seq': seq, 'ref_positions': ref_positions,'skipped_intervals':skipped_intervals}
        if master_read is None:
            master_read = new_read
        master_read = merge_reads(master_read, new_read)
    master_read['SN'] = read.reference_name
    master_read['is_reverse'] = read_d[0].is_reverse
    ends = list(read_d[0].get_tag('EL'))
    ends.extend(list(read_d[1].get_tag('EL')))
    master_read['ends'] = list(set(ends))
    master_read['ref_intervals'] = get_ref_intervals(master_read)
    ref_skip_union = (master_read['ref_intervals'] | master_read['skipped_intervals'])
    master_read['del_intervals'] =  get_del_intervals(ref_skip_union)
    del master_read['ref_positions']
    master_read['phred'] = np.rint(-10*np.log10(master_read['p_x']))
    del master_read['p_x']
    master_read['seq'] = ''.join(master_read['seq'])
    master_read['NR'] = read_d[0].get_tag('NR')+read_d[1].get_tag('NR')
    return master_read

def combine_umi_and_internal(combined_read_dict):
    combined_merged_dict = {}
    for cell in combined_read_dict.keys():
        combined_merged_dict[cell] = {}
        for gene in combined_read_dict[cell].keys():
            combined_merged_dict[cell][gene] = {}
            for umi in combined_read_dict[cell][gene].keys():
                #print(cell,gene,umi)
                combined_merged_dict[cell][gene][umi] = stitch_internal(combined_read_dict[cell][gene][umi])
    return combined_merged_dict

def iter_over_internal_reads(merged_dict):
    for cell in merged_dict:
        #print(cell)
        for gene in merged_dict[cell]:
            #print('\t', gene)
            for umi in merged_dict[cell][gene]:
                #print('\t\t', umi)
                for read in  merged_dict[cell][gene][umi]:
                    yield read, cell, gene, umi

def refine_read_dict(r_d, exclude_d):
    new_d = {}
    for read,cell,gene,umi in iter_over_internal_reads(r_d):
        cell = read.get_tag('BC')
        gene = read.get_tag('XT')
        umi = read.get_tag('UB')
        if cell in new_d.keys() and umi == '' and gene in exclude_d[cell].keys():
            if gene in new_d[cell].keys():
                if umi in new_d[cell][gene].keys():
                    if read.query_name not in [r for l in [v for k,v in exclude_d[cell][gene].items()] for r in l]:
                        new_d[cell][gene][umi].append(read)
                else:
                    if read.query_name not in [r for l in [v for k,v in exclude_d[cell][gene].items()] for r in l]:
                        new_d[cell][gene][umi] = [read]
            else:
                if read.query_name not in [r for l in [v for k,v in exclude_d[cell][gene].items()] for r in l]:
                    new_d[cell][gene] = {}
                    new_d[cell][gene][umi] = [read]
        elif umi == '' and gene in exclude_d[cell].keys():
            if read.query_name not in [r for l in [v for k,v in exclude_d[cell][gene].items()] for r in l]:
                new_d[cell] = {}
                new_d[cell][gene] = {}
                new_d[cell][gene][umi] = [read]
    return new_d

def add_internal_reads(infile, origfile, outfile):
    print('Loading reads')
    assembled_mols = pysam.AlignmentFile(infile)
    orig_reads = pysam.AlignmentFile(origfile)

    assembled_read_dict = make_assembled_dict(assembled_mols, None, read_dict={})
    internal_reads = make_read_dict(orig_reads, None, read_dict={})

    print('Finding overlapping internal reads')
    internal_read_names = find_internal_reads(internal_reads, assembled_read_dict)
    internal_read_dict = make_internal_read_dict(internal_reads, internal_read_names)
    nreads = len(list(iter_over_internal_reads(internal_read_dict)))
    print('Reads found this iteration: {}'.format(nreads))

    print('Merging internal reads')
    internal_merged_dict = make_merged_dict(internal_read_dict, spec_strand=True, mol_dict=assembled_read_dict)

    print('Saving merged internal reads')
    write_sam_file(internal_merged_dict, '{}_internal_mols_0.sam'.format(outfile), assembled_mols)
    internal_mols = pysam.AlignmentFile('{}_internal_mols_0.sam'.format(outfile))

    print('Combining UMI reads and internal reads')
    internal_mol_read_dict = make_combined_dict(internal_mols,None,read_dict={})

    assembled_mols = pysam.AlignmentFile(infile)
    combined_read_dict = make_combined_dict(assembled_mols, None, internal_mol_read_dict)

    for cell,gene, umi in to_delete(combined_read_dict):
        del combined_read_dict[cell][gene][umi]

    combined_merged_dict = combine_umi_and_internal(combined_read_dict)

    write_sam_file(combined_merged_dict, '{}_combined_mols_0.sam'.format(outfile), assembled_mols)

    i = 1
    while nreads > 0:
        
        print('Filtering reads for iteration {}'.format(i))
        internal_reads = refine_read_dict(internal_reads, exclude_d=internal_read_names)
        assembled_mols = pysam.AlignmentFile('{}_combined_mols_{}.sam'.format(outfile,i-1))
        assembled_read_dict = make_assembled_dict(assembled_mols, None, read_dict={})
        
        print('Finding overlapping internal reads')
        internal_read_names = find_internal_reads(internal_reads, assembled_read_dict)
        internal_read_dict = make_internal_read_dict(internal_reads, internal_read_names)
        nreads = len(list(iter_over_internal_reads(internal_read_dict)))
        print('Reads found this iteration: {}'.format(nreads))
        if nreads == 0:
            break
        print('Merging internal reads')
        internal_merged_dict = make_merged_dict(internal_read_dict, spec_strand=True, mol_dict=assembled_read_dict)
        
        print('Saving merged internal reads')
        write_sam_file(internal_merged_dict, '{}_internal_mols_{}.sam'.format(outfile,i), assembled_mols)
        internal_mols = pysam.AlignmentFile('{}_internal_mols_{}.sam'.format(outfile,i))

        print('Combining UMI reads and internal reads')
        internal_mol_read_dict = make_combined_dict(internal_mols,None,read_dict={})
        assembled_mols = pysam.AlignmentFile(infile)
        combined_read_dict = make_combined_dict(assembled_mols, None, internal_mol_read_dict)

        for cell,gene, umi in to_delete(combined_read_dict):
            del combined_read_dict[cell][gene][umi]

        combined_merged_dict = combine_umi_and_internal(combined_read_dict)

        write_sam_file(combined_merged_dict,'{}_combined_mols_{}.sam'.format(outfile,i) , assembled_mols)
        i += 1
    print('Iteration done')
    print('Merging .sam files')
    mol_names = []
    outbam = pysam.AlignmentFile("{}_all.bam".format(outfile), "wb", template=assembled_mols)
    for j in range(i-1, -1, -1):
        print('Adding molecules from {}_combined_mols_{}.sam'.format(outfile,j))
        combined_mols = pysam.AlignmentFile('{}_combined_mols_{}.sam'.format(outfile,j))
        for mol in combined_mols.fetch(None):
            if mol.query_name not in mol_names and not mol.is_unmapped:
                mol_names.append(mol.query_name)
                mol.set_tag('RI',j+1)
                outbam.write(mol)
    print('Adding molecules from {}'.format(infile))
    combined_mols = pysam.AlignmentFile(infile)
    for mol in combined_mols.fetch(None):
        if mol.query_name not in mol_names:
            mol_names.append(mol.query_name)
            mol.set_tag('RI',0)
            outbam.write(mol)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add internal reads to stitched molecules')
    parser.add_argument('--i',metavar='input', type=str, nargs=1, help='Input .bam file with stitched molecules')
    parser.add_argument('--internal', metavar='internal', type=str, nargs=1, help='Input .bam file with internal reads')
    parser.add_argument('--o', metavar='output', type=str, nargs=1, help='Output prefix for .sam files')
    args = parser.parse_args()
    infile = args.i[0]
    outfile = args.o[0]
    origfile = args.internal[0]
    add_internal_reads(infile, origfile,outfile)

