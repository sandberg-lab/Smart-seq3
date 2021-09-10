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


def convert_ref_to_dict(ref):
    
    ref_dict = defaultdict(dict)
    for i in ref.index:
        bool_array = np.zeros(int(ref.iloc[i]['Total_n_exons']))
        ex_idx = [int(ii)-1 for ii in ref.iloc[i]['Exon_Index'].split(',')]
        bool_array[ex_idx] = 1
        
        if 'Transcripts' not in ref_dict[ref.iloc[i]['Gene']].keys():
            ref_dict[ref.iloc[i]['Gene']]['Transcripts'] = defaultdict(dict)
            ref_dict[ref.iloc[i]['Gene']]['Total_n_exons'] = ''
        
        ref_dict[ref.iloc[i]['Gene']]['Transcripts'][ref.iloc[i]['Transcript']] = {'Exon_bool_array': bool_array, 'Exon_Loc': ref.iloc[i]['Exon_Loc'],
                                                                                  'Junction': ref.iloc[i]['Junction']}
        ref_dict[ref.iloc[i]['Gene']]['Total_n_exons'] = int(ref.iloc[i]['Total_n_exons'])
        ref_dict[ref.iloc[i]['Gene']]['chrom'] = ref.iloc[i]['chrom']
        
    return ref_dict

def get_overlapping_genes(gff_merged):
    
    overlaped_df = gff_merged.loc[gff_merged[8]!=gff_merged[17]]
    aa = overlaped_df[8].str.split(';', expand=True)[0].replace('gene_id=','',regex=True).to_list()
    bb = overlaped_df[17].str.split(';', expand=True)[0].replace('gene_id=','',regex=True).to_list()
    df = pd.DataFrame([aa,bb],index=['A','B']).T
    
    overlaped_gene_dict = df.groupby(by='A').apply(lambda x: x['B'].to_list()).to_dict()
    
    return overlaped_gene_dict

def _filter_reads_from_other_gene(raw_aligned, gene_neighbors, indir, gene):
    
    gene_neighbors = gene_neighbors + [gene]
    
    gneighbor = ref_iso.loc[ref_iso['Gene'].isin(gene_neighbors)][['Transcript','Gene','Exon_Loc']]
    gneighbor['Exon_Loc'] = gneighbor['Exon_Loc'].replace(';',',',regex=True)
    gneighbor = gneighbor.assign(Exon_Loc=gneighbor.Exon_Loc.str.split(',')).explode('Exon_Loc')
    
    gneighbor_bed = gneighbor.Exon_Loc.str.split('-', expand=True)
    gneighbor_bed.insert(0,'chrom','chr1')
    gneighbor_bed.insert(3,'Transcript',gneighbor['Transcript'])
    gneighbor_bed.insert(4,'Gene',gneighbor['Gene'])
    gneighbor_bed.to_csv('%s/.tmp/_%s_neighbors.bed' %(indir, gene), sep="\t", index=False)
    gneighbor_bed_obj = pybedtools.BedTool('%s/.tmp/_%s_neighbors.bed' %(indir, gene))
    
    aligned_info = raw_aligned.iloc[:,[0,1,12]]
    aligned_info.columns = ['Read','Flag','Region']
    aligned_info = aligned_info.assign(Region=aligned_info.Region.str.split(',')).explode('Region')
    aligned_info_bed = aligned_info.Region.str.split('-', expand=True)
    aligned_info_bed.insert(0,'chrom','chr1')
    aligned_info_bed.insert(3,'Read',aligned_info['Read'])
    aligned_info_bed.insert(4,'Flag',aligned_info['Flag'])
    aligned_info_bed.to_csv('%s/.tmp/_%s_raw_aligned.bed' %(indir, gene), sep="\t", index=False)
    aligned_info_bed_obj = pybedtools.BedTool('%s/.tmp/_%s_raw_aligned.bed' %(indir, gene))
    
    tmp = aligned_info_bed_obj.intersect(gneighbor_bed_obj, wo=True)
    my_intersect_all = tmp.to_dataframe().drop_duplicates()
    read_overlaps = my_intersect_all.groupby(by=['name','itemRgb','blockCount']).apply(lambda x: np.sum(x['blockSizes']))
    read_overlaps = read_overlaps.reset_index()
    read_overlaps.columns = list(read_overlaps.columns)[:-1] + ['Len']
    
    max_frag_overlap = read_overlaps.groupby(by='name').apply(lambda x: x.loc[x['Len']==x['Len'].max()])
    
    kept_read_names = list(set(max_frag_overlap.loc[max_frag_overlap['blockCount']==gene]['name'].values))
    kept_reads = raw_aligned.loc[raw_aligned[0].isin(kept_read_names)]
    
    return kept_reads

def correct_bool_array(bool_array, junc_list):
    
    multi_ex_juncs = [junc for junc in junc_list if len(junc.split('^'))>2]
    junc_list = [junc for junc in junc_list if junc not in multi_ex_juncs]
    
    if len(junc_list) > 0:
        junc_df = pd.DataFrame(junc_list, columns=['Start']).Start.str.split('^',expand=True)
        junc_df.columns = ['Start','End']
        ambig_junc_same_start = junc_df.groupby(by="Start").apply(lambda x: x.apply(lambda y: '%s^%s' %(y[0],y[1]), axis=1).to_list())
        ambig_juncs1 = ambig_junc_same_start.loc[ambig_junc_same_start.apply(lambda x: len(x)) > 1].to_list()
        ambig_juncs1 = sum(ambig_juncs1,[])
        ambig_junc_same_stop = junc_df.groupby(by="End").apply(lambda x: x.apply(lambda y: '%s^%s' %(y[0],y[1]), axis=1).to_list())
        ambig_juncs2 = ambig_junc_same_stop.loc[ambig_junc_same_stop.apply(lambda x: len(x)) > 1].to_list()
        ambig_juncs2 = sum(ambig_juncs2,[])
        ambig = list(set(ambig_juncs1 + ambig_juncs2))
        junc_list = list(set(junc_list) - set(ambig))
        
    junc_list = junc_list + multi_ex_juncs
    junc_ex_idx = [list(map(int,junc.split('^'))) for junc in junc_list]
    skipped_exon_idx = [list(set(range(np.min(ex_idx),np.max(ex_idx)+1)) - set(ex_idx)) for ex_idx in junc_ex_idx]
    skipped_exon_idx = list(set(sum(skipped_exon_idx,[])))
    
    flag = np.isnan(bool_array)
    skipped_exon_idx = list(set(np.array(range(len(bool_array)))[flag]+1) & set(skipped_exon_idx))
    
    if len(skipped_exon_idx) > 0:
        bool_array[np.array(skipped_exon_idx) - 1] = 0
    
    return bool_array

def _infer_isoform(exon_bool_array, ref):
    
    scores = []
    for trans in ref['Transcripts'].keys():
        curr_ref_iso_bool = ref['Transcripts'][trans]['Exon_bool_array']
        flag = exon_bool_array == curr_ref_iso_bool
        scores.append(len(flag[flag]))
        
    max_score = np.max(scores)
    infered_isoforms = np.array(list(ref['Transcripts'].keys()))[np.where(scores==max_score)]
    n_infered = len(infered_isoforms)
    
    infered_ref_transcript_bool_string = [''.join(list(map(str,map(int,ref['Transcripts'][iso]['Exon_bool_array'])))) for iso in infered_isoforms]
    
    return {'max_score': max_score, 'infered_transcripts': ','.join(infered_isoforms), 'n_infered_transcripts': n_infered,
            'total_n_ref_transcripts': len(ref['Transcripts'].keys()), 'infered_ref_transcript_bool_string': ','.join(infered_ref_transcript_bool_string)}

def _isoform_inference_of_single_molec(aligned_reads_df, ref):
    
    mapped_ex_junc = list(set(aligned_reads_df[14].values))
    exon_idx_list = [str(ii) for ii in mapped_ex_junc if not re.search(',', str(ii))]
    
    junc_list = [ii.replace(',','^') for ii in mapped_ex_junc if re.search(',', str(ii))]
    multi_ex_juncs = [junc for junc in junc_list if len(junc.split('^'))>2]
    multi_ex_juncs_nn = [['^'.join(junc.split('^')[idx:(idx+2)]) for idx in range(len(junc.split('^'))-1)] for junc in multi_ex_juncs]
    multi_ex_juncs_nn = sum(multi_ex_juncs_nn, [])
    junc_list = [junc for junc in junc_list if junc not in multi_ex_juncs]
    junc_list = list(set(junc_list + multi_ex_juncs_nn))
    rm_junc1, rm_junc2 = [[],[]]
    
    if len(junc_list) > 0:
        junc_df = pd.DataFrame(junc_list, columns=['Start']).Start.str.split('^',expand=True)
        junc_df.columns = ['Start','End']
        ambig_junc_same_start = junc_df.groupby(by="Start").apply(lambda x: x.apply(lambda y: '%s^%s' %(y[0],y[1]), axis=1).to_list())
        ambig_juncs = ambig_junc_same_start.loc[ambig_junc_same_start.apply(lambda x: len(x)) > 1].to_list()
        no_ambig_juncs = ambig_junc_same_start.loc[ambig_junc_same_start.apply(lambda x: len(x)) == 1].to_list()
        filtered = [[junc for junc in juncL if junc.split('^')[1] in exon_idx_list] for juncL in ambig_juncs]
        rm_junc1 = [[junc for junc in juncL if junc.split('^')[1] not in exon_idx_list] for juncL in ambig_juncs]
        junc_list = filtered + no_ambig_juncs
        junc_list = sum(junc_list,[])
        
    if len(junc_list) > 0:
        junc_df = pd.DataFrame(junc_list, columns=['Start']).Start.str.split('^',expand=True)
        junc_df.columns = ['Start','End']
        ambig_junc_same_stop = junc_df.groupby(by="End").apply(lambda x: x.apply(lambda y: '%s^%s' %(y[0],y[1]), axis=1).to_list())
        ambig_juncs = ambig_junc_same_stop.loc[ambig_junc_same_stop.apply(lambda x: len(x)) > 1].to_list()
        no_ambig_juncs = ambig_junc_same_stop.loc[ambig_junc_same_stop.apply(lambda x: len(x)) == 1].to_list()
        filtered = [[junc for junc in juncL if junc.split('^')[0] in exon_idx_list] for juncL in ambig_juncs]
        rm_junc2 = [[junc for junc in juncL if junc.split('^')[0] not in exon_idx_list] for juncL in ambig_juncs]
        junc_list = filtered + no_ambig_juncs
        junc_list = sum(junc_list,[])
        
    rm_junc_list = [junc.replace('^',',') for junc in set(sum(rm_junc1+rm_junc2,[]))]
    if len(rm_junc_list)>0:
        flag_df = pd.concat([pd.DataFrame(aligned_reads_df[14].str.match(junc,na=False)) for junc in rm_junc_list], axis=1)
        aligned_reads_df = aligned_reads_df.loc[flag_df.sum(axis=1) == 0]
                 
    ex_from_junc = [junc.split('^') for junc in junc_list]
    ex_from_junc = list(set(sum(ex_from_junc,[])))
    
    exon_idx_list = list(set(exon_idx_list + ex_from_junc))
    exon_idx_list = list(map(int, exon_idx_list))
    exon_idx_list.sort()
    
    read_coord = aligned_reads_df.groupby(by=0).apply(lambda x: '|'.join(x[12]))
    n_fragment = read_coord.shape[0]
    read_coord_list = ';'.join((list(read_coord.values)))
    
    exon_bool_array = np.zeros(int(ref['Total_n_exons']))
    exon_bool_array[:] = np.nan
    exon_bool_array[[ii-1 for ii in exon_idx_list]] = 1
    if len(junc_list)>0:
        exon_bool_array = correct_bool_array(exon_bool_array, junc_list)
    
    infered = _infer_isoform(exon_bool_array, ref)
    exon_bool_string = ''.join(['N' if np.isnan(bb) else str(int(bb)) for bb in exon_bool_array])
    
    if aligned_reads_df.shape[0] == 0: return []
    out = [aligned_reads_df[16].iloc[0], aligned_reads_df[17].iloc[0], 
           n_fragment, ','.join(map(str,exon_idx_list)), ','.join(junc_list),
           read_coord_list, aligned_reads_df[13].iloc[0], ref['Total_n_exons'], infered['total_n_ref_transcripts'],
           infered['infered_transcripts'],
           infered['n_infered_transcripts'], infered['infered_ref_transcript_bool_string'],
           exon_bool_string, infered['max_score']]
    
    return out

def _run_isoform(indir, outdir, ref_iso_dict, kept_cell_BCs, conf, overlaped_gene_dict, gene):
    
    print(gene)
    
    if os.path.exists('%s/%s/%s' %(outdir, ref_iso_dict[gene]['chrom'], gene)): return
    raw_aligned_reads = pd.read_table('%s/keptReads/%s/%s_aligned_reads.csv' %(indir, ref_iso_dict[gene]['chrom'], gene), header=None, index_col=None, sep="\t")
    
    if gene in overlaped_gene_dict.keys():
        aligned_reads = _filter_reads_from_other_gene(raw_aligned_reads, overlaped_gene_dict[gene], indir, gene)
    else:
        aligned_reads = raw_aligned_reads
    
    if aligned_reads.shape[0] == 0: return
    chrom = aligned_reads.iloc[0,2]
    outdir = '%s/%s' %(outdir, chrom)
    if not os.path.exists(outdir): os.makedirs(outdir)
    
    bc = aligned_reads.apply(lambda x: '%s+%s' %(x[16],x[17]), axis=1)
    bc.name = 'BC_UB'
    aligned = pd.concat([aligned_reads, bc], axis=1)
    
    results = aligned.groupby(by='BC_UB').apply(_isoform_inference_of_single_molec, ref_iso_dict[gene])
    df = pd.DataFrame(list(results.values)).dropna()
    df.to_csv('%s/%s' %(outdir, gene), sep="\t", index=False, header=False)
    
    return


def get_junction(ass, trans_df):
    
    tt = pd.DataFrame(ass.coordinates.str.split(';').to_list(), index=pd.MultiIndex.from_frame(ass[['Exon_Idx','flag','Transcripts']])).stack()
    tt = tt.reset_index()
    tt.columns = ['Exon_Idx','flag','Transcripts','rm','coordinates']
    ass = tt[ass.columns]
    
    curr_ass = pd.concat([pd.DataFrame(ass.coordinates.str.split('-').tolist()), pd.DataFrame(ass.Transcripts.str.split(',').tolist())], axis=1)
    curr_ass.columns = ['start','end','transcript','exon_idx']
    
    ass_start = curr_ass.groupby(by="exon_idx").apply(lambda x: len(set(x['start'].values)))
    ass_end = curr_ass.groupby(by="exon_idx").apply(lambda x: len(set(x['end'].values)))
    
    ass_start_exid = list(ass_start[ass_start>1].index)
    ass_end_exid = list(ass_end[ass_end>1].index)
    
    max_n_exons = trans_df['Exon_Idx'].max()
    
    ass_start_exid = [eid for eid in ass_start_exid if eid!='1']
    ass_end_exid = [eid for eid in ass_end_exid if eid!=str(max_n_exons)]
    
    ass_start_junc = None
    ass_end_junc = None
    
    if len(ass_start_exid) > 0:
        tmp = curr_ass.loc[curr_ass["exon_idx"].isin(ass_start_exid)]
        ass_start_junc = tmp.apply(_get_junc_start, axis=1, trans_df=trans_df)
        ass_start_junc = pd.DataFrame(list(ass_start_junc[~ass_start_junc.isnull()].values))
    
    if len(ass_end_exid) > 0:
        tmp = curr_ass.loc[curr_ass["exon_idx"].isin(ass_end_exid)]
        ass_end_junc = tmp.apply(_get_junc_end, axis=1, trans_df=trans_df)
        ass_end_junc = pd.DataFrame(list(ass_end_junc[~ass_end_junc.isnull()].values))
    
    return {'ass_start': ass_start_junc, 'ass_end': ass_end_junc}

def _get_junc_start(x, trans_df):
    
    row_idx = list(trans_df.query('Exon_Idx=="%s" and Transcripts=="%s"' %(x[3], x[2])).index)[0]
    min_ex_idx = trans_df.query('Transcripts=="%s"' %(x[2]))['Exon_Idx'].min()
    if int(x['exon_idx'])==min_ex_idx: return
    
    junc_ex_pos = trans_df.loc[row_idx-1]['coordinates'].split('-')[1]
    junc_idx = '%s^%s' %(trans_df.loc[row_idx-1]['Exon_Idx'], x[3])
    
    return [x[2], x[3], junc_idx, '%s,%s' %(junc_ex_pos, x[0])]

def _get_junc_end(x, trans_df):
    
    row_idx = list(trans_df.query('Exon_Idx=="%s" and Transcripts=="%s"' %(x[3], x[2])).index)[0]
    max_ex_idx = trans_df.query('Transcripts=="%s"' %(x[2]))['Exon_Idx'].max()
    if int(x['exon_idx'])==max_ex_idx: return
    
    junc_ex_pos = trans_df.loc[row_idx+1]['coordinates'].split('-')[0]
    junc_idx = '%s^%s' %(x[3],trans_df.loc[row_idx+1]['Exon_Idx'])
    
    return [x[2], x[3], junc_idx, '%s,%s' %(x[1], junc_ex_pos)]

def isoform_inference_correction_by_ass_v2(expr_indir, ref, outdir, gene_file):
    
    chrom, gene = gene_file.split('/')[-2:]
    print(gene)
    
    outdir = '%s/%s' %(outdir, chrom)
    if not os.path.exists(outdir): os.makedirs(outdir)
    
    if os.stat(gene_file).st_size == 0: return
    initial_infered = pd.read_table(gene_file, header=None, index_col=None, sep="\t")
    initial_infered.index = initial_infered.apply(lambda x: '%s_%s' %(x[0], x[1]), axis=1)
    infered_to_correct = initial_infered.loc[initial_infered[10]>1]
    if initial_infered.iloc[0,8]==1 or infered_to_correct.shape[0]==0:
        initial_infered[14] = ['no' for ii in range(initial_infered.shape[0])]
        initial_infered[[0,1,3,5,9,10,12]].to_csv('%s/%s' %(outdir, gene), sep="\t", index=False, header=False)
        return
    
    trans_list = []
    for trans in ref[gene]['Transcripts'].keys():
        exon_bool = list(ref[gene]['Transcripts'][trans]['Exon_bool_array'])
        exon_idx = [ii+1 for ii in range(len(exon_bool))]
        
        exon_idx = pd.DataFrame([exon_idx, exon_bool], index=['Exon_Idx','flag']).T.query('flag==1')
        exon_idx['Exon_Idx'] = exon_idx['Exon_Idx'].astype(int)
        exon_idx.index = range(exon_idx.shape[0])
        exon_idx = pd.concat([exon_idx, pd.DataFrame(ref[gene]['Transcripts'][trans]['Exon_Loc'].split(','), columns=['coordinates']), pd.DataFrame([trans for ii in range(exon_idx.shape[0])], columns=['Transcripts'])], axis=1)
        trans_list.append(exon_idx)
    
    trans_df = pd.concat(trans_list, axis=0)
    trans_df.index = range(trans_df.shape[0])
    ass = trans_df.groupby(by='Exon_Idx').apply(get_comm_exon_ass).dropna(how='all')
    if ass.shape[0] == 0:
        initial_infered[14] = ['no' for ii in range(initial_infered.shape[0])]
        initial_infered[[0,1,3,5,9,10,12]].to_csv('%s/%s' %(outdir, gene), sep="\t", index=False, header=False)
        return
    
    ass['Exon_Idx'] = ass['Exon_Idx'].astype(int)
    ass['Transcripts'] = ass.apply(lambda x: '%s,%s' %(x[-1],x[0]), axis=1)
    new_ass = pd.DataFrame(ass.coordinates.str.split(';').tolist(), index=ass.Transcripts).stack()
    new_ass = new_ass.reset_index([0, 'Transcripts'])
    ass_exon_reg = pd.DataFrame([['chr1']+coord.split('-')+['.'] for coord in new_ass.iloc[:,1].values])
    ass_exon_reg[4] = new_ass['Transcripts']
    ass_exon_reg.to_csv('%s/../.tmp/_%s_ass' %(outdir, gene), sep="\t", index=False, header=False)
    ass_exon_reg_bed = pybedtools.BedTool('%s/../.tmp/_%s_ass' %(outdir, gene))
    
    ass_junc = get_junction(ass, trans_df)
    
    aligned_list = []
    for idx in infered_to_correct.index:
        region = pd.DataFrame([['chr1']+reg.split('-')+['.'] for reg in initial_infered.loc[idx][5].replace(';',',').replace('|',',').split(',')])
        region[4] = [idx for i in range(region.shape[0])]
        aligned_list.append(region)
    aligned_reg = pd.concat(aligned_list, axis=0)
    aligned_reg.to_csv('%s/../.tmp/_aligned_region_in_%s' %(outdir, gene), sep="\t", index=False, header=False)
    aligned_reg_bed = pybedtools.BedTool('%s/../.tmp/_aligned_region_in_%s' %(outdir, gene))
    
    tmp = aligned_reg_bed.intersect(ass_exon_reg_bed, wo=True)
    if os.stat(tmp.fn).st_size==0:
        initial_infered[14] = ['no' for ii in range(initial_infered.shape[0])]
        initial_infered[[0,1,3,5,9,10,12]].to_csv('%s/%s' %(outdir, gene), sep="\t", index=False, header=False)
        
        return
    
    intersect = tmp.to_dataframe()
    intersect = intersect.drop_duplicates()
    intersect[['trans','exon_idx']] = intersect.iloc[:,-2].str.split(',',expand=True)
       
    trans_idx = intersect.groupby(by='score').apply(_get_max_overlap_transcript, infered_to_correct, ass_junc)
    trans_counts = trans_idx.apply(lambda x: len(x.split(',')))
    
    infered = initial_infered.copy()
    infered.loc[trans_idx.index,9] = trans_idx
    infered.loc[trans_idx.index,10] = trans_counts
    infered[14] = ['no' for i in range(infered.shape[0])]
    infered.loc[trans_idx.index,14] = 'yes'
    infered_out = infered[[0,1,3,5,9,10,12]]
    
    infered_out.to_csv('%s/%s' %(outdir, gene), sep="\t", index=False, header=False)
   
    return


def score_junction_mapping(infered_trans, junc_r, ass_junc):
    
    if len(junc_r) == 0: return infered_trans
    
    ass_start = pd.DataFrame(ass_junc['ass_start'], columns=[0,1,2,3])
    ass_end = pd.DataFrame(ass_junc['ass_end'], columns=[0,1,2,3])
    
    ass_start = ass_start.loc[ass_start[0].isin(infered_trans)]
    ass_end = ass_end.loc[ass_end[0].isin(infered_trans)]
    
    mapping1 = ass_start.loc[ass_start[3].isin(junc_r)]
    mapping2 = ass_end.loc[ass_end[3].isin(junc_r)]
    
    if mapping1.shape[0]==0 and mapping2.shape[0]==0:
        return infered_trans
    
    mapped1 = pd.DataFrame(mapping1.groupby(by=1).apply(lambda x: list(x[0].values)),columns=[0])
    mapped2 = pd.DataFrame(mapping2.groupby(by=1).apply(lambda x: list(x[0].values)),columns=[0])
    
    trans_list = []
    if mapped1.shape[0] > 0:
        trans_list.extend(mapped1[0].sum())
    
    if mapped2.shape[0] > 0:
        trans_list.extend(mapped2[0].sum())
    
    score = pd.Series(trans_list).value_counts()
    infered_trans =  list(score[score==score.max()].index)
    
    return infered_trans

def _get_overlap_len(df):
    
    set_list = [set(range(df.iloc[ii]['start'],df.iloc[ii]['end'])) for ii in range(df.shape[0])]
    base_overlap = len(set.union(*set_list))
    return base_overlap

def _get_max_overlap_transcript(x, infered_to_correct, ass_junc):
    
    infered_trans = list(set(infered_to_correct.loc[x.iloc[0,4]][9].split(',')))
    xx = x.loc[x['trans'].isin(infered_trans)]
    if xx.shape[0] == 0:
        junc_r = [reg for reg in infered_to_correct.loc[x.iloc[0,4]][5].split('-') if re.search(',', reg)]
        corr_trans = score_junction_mapping(infered_trans, junc_r, ass_junc)
        return ','.join(corr_trans)
    
    overlap_len = pd.DataFrame(xx.groupby(by='blockCount').apply(_get_overlap_len), columns=['len'])
    overlap_len = pd.concat([overlap_len, pd.DataFrame(overlap_len.index, index=overlap_len.index)['blockCount'].str.split(',',expand=True)], axis=1)
    overlap_len.columns = ['len','trans','exon_idx']
    
    tmp = overlap_len.groupby(by='exon_idx').apply(lambda x: x.loc[x['len']==np.max(x['len'])])
    tmp.index.names = ['idx1','idx2']
    trans = tmp.groupby(by='exon_idx').apply(lambda x: list(set(x['trans'].values)))
    trans_len = trans.apply(len)
    corr_trans = list(set(trans[trans_len[trans_len==np.min(trans_len)].index].sum()))
    
    if len(corr_trans) > 1:
        junc_r = [reg for reg in infered_to_correct.loc[x.iloc[0,4]][5].split('-') if re.search(',', reg)]
        corr_trans = score_junction_mapping(corr_trans, junc_r, ass_junc)
    
    out = ','.join(corr_trans)
    
    return out


def get_comm_exon_ass(df):
    
    if df.shape[0] == 1: return
    if len(set(df['coordinates'].values)) == 1: return
    
    return df

def get_isoforms(conf, indir, ref):
    
    global ref_iso
    ref_iso = ref
    
    os.system('bedtools intersect -s -a %s/gene.gff -b %s/gene.gff -wo > %s/gene_merged.gff' %(indir, indir, indir))
    df = pd.read_table('%s/gene_merged.gff' %(indir), header=None, sep="\t", index_col=None)
    overlaped_gene_dict = get_overlapping_genes(df)
    
    outdir = '%s/.R1' %(indir)
    if not os.path.exists(outdir): os.makedirs(outdir)
    
    kept_cell_BCs = list(pd.read_table(conf['annotation']['zumi_keptbarcode'],header=0, index_col=0, sep=",").index)
    ref_iso_dict = convert_ref_to_dict(ref_iso)
    
    genes = list(set(ref_iso['Gene'].values))
    gene_files = glob.glob('%s/.R1/*/*' %(indir))
    infered_genes = [val.split('/')[-1] for val in gene_files if not re.search('_log', val)]
    remain_genes = list(set(genes) - set(infered_genes))
    print('%s remaining files' %(len(remain_genes)))
    
    os.system('> %s/_log' %(outdir))
    
    pool = mp.Pool(processes=int(conf['expression']['nproc']))
    func = partial(_run_isoform, indir, outdir, ref_iso_dict, kept_cell_BCs, conf, overlaped_gene_dict)
    pool.map(func, remain_genes, chunksize=1)
    pool.close()
    
    gene_files = glob.glob('%s/.R1/*/*' %(indir))
    infered_gene_paths = [val for val in gene_files if not re.search('_log', val)]
    outdir = '%s/assigned_isoforms' %(indir)
    if not os.path.exists(outdir): os.makedirs(outdir)
    
    infered_gene_paths = [gene_file for gene_file in infered_gene_paths if not os.path.exists(gene_file.replace('.R1','assigned_isoforms')) or os.stat(gene_file.replace('.R1','assigned_isoforms')).st_size==0]
    print('%s remaining files' %(len(infered_gene_paths)))
    
    if not os.path.exists('%s/.tmp' %outdir): os.makedirs('%s/.tmp' %outdir)
    pool = mp.Pool(processes=int(conf['expression']['nproc']))
    func = partial(isoform_inference_correction_by_ass_v2, indir, ref_iso_dict, outdir)
    pool.map(func, infered_gene_paths, chunksize=1)
    pool.close()
        
    p = subprocess.Popen('rm -rf %s/.tmp' %(outdir), shell=True)
    (output, err) = p.communicate()
    
    p = subprocess.Popen('rm -rf %s/.tempDir' %(indir), shell=True)
    (output, err) = p.communicate()
    
    p = subprocess.Popen('rm -rf %s/.R1' %(indir), shell=True)
    (output, err) = p.communicate()
    
    
    return
