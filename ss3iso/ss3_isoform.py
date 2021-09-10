#!/usr/bin/python
# Developer: Ping Chen
# Contact: ping.chen@ki.se
# Date: 2020-01-10
# Version: 0.1.3

# ------------------------------------------------ #
#      SS3 isoform reconstruction pipeline         #
# ------------------------------------------------ #
import os
from optparse import OptionParser
import glob
import configparser
import re
from pyModule.informative_reads import *
from pyModule.reference import *
from pyModule.isoform_reconstruct import *
import pybedtools

def main():
    
    parser=OptionParser()
    
    parser.add_option('-i', '--inputBAM', dest='inputBAM', 
                      help='Aligned BAM from zUMI filtering+mapping steps with cell barcode and umi barcode correction.')
    
    parser.add_option('-c', '--config', dest='config', 
                      help='A configuration file for required files and parameters.')
    
    parser.add_option('-e', '--experiment', dest='experiment', 
                      help='Experiment name.')
    
    parser.add_option('-o', '--outputDir', dest='outputDir', default='ss3rnaseq',
                      help='The output directory for the experiment.')
    
    parser.add_option('-p', '--process', dest='process', default=8,
                      help='The number of processes for parallel computing.')
    
    parser.add_option('-s', '--species', dest='species', default='hg38',
                      help='The species under study.')
    
    parser.add_option("-P", "--Preprocess", action="store_true", dest='preprocess',
                      help="Preprocess the input BAM for downstream analysis.")
    
    parser.add_option("-R", "--Reconstruction", action="store_true", dest='reconstruction',
                      help="Run isoform reconstruction.")


    (op, args) = parser.parse_args()
    inputBAM = op.inputBAM
    conf = op.config
    experiment = op.experiment
    outdir = op.outputDir
    nprocess = int(op.process)

    if op.species == 'hg38' or op.species == 'hg19': species = 'hsa'
    elif op.species == 'mm9' or op.species == 'mm10': species = 'mmu'
    
    config = configparser.ConfigParser()
    config.read(conf)
    conf_data = config._sections
    
    if not os.path.exists(outdir): os.makedirs(outdir)
    if not os.path.exists('%s/%s' %(outdir, species)): os.makedirs('%s/%s' %(outdir, species))
    if not os.path.exists('%s/%s/%s' %(outdir, species, experiment)): os.makedirs('%s/%s/%s' %(outdir, species, experiment))
    
    umi_file_prefix = 'UBfix.sort.bam'
    if op.preprocess:
        print('Preprocessing on input BAM ...')
        preDir = os.path.join(outdir, species, experiment, "preprocess")
        if not os.path.exists(preDir): os.makedirs(preDir)
        
        cmd = 'samtools sort -m %s -O bam -@ %s -o %s/%s %s' %(conf_data['preprocess']['memory'], nprocess, preDir, re.sub(umi_file_prefix,'UBfix.coordinateSorted.bam',os.path.basename(inputBAM)), inputBAM)
        os.system(cmd)
        
        cmd = 'samtools view -b -q 255 %s/%s > %s/%s' %(preDir, re.sub(umi_file_prefix,'UBfix.coordinateSorted.bam',os.path.basename(inputBAM)), preDir, re.sub(umi_file_prefix,'UBfix.coordinateSorted_unique.bam',os.path.basename(inputBAM)))
        os.system(cmd)
        
        cmd = 'samtools view -h %s/%s | awk \'$12 != "NH:i:1"\' | samtools view -bS - > %s/%s' %(preDir, re.sub(umi_file_prefix,'UBfix.coordinateSorted.bam',os.path.basename(inputBAM)), preDir, re.sub(umi_file_prefix,'UBfix.coordinateSorted_multi.bam',os.path.basename(inputBAM)))
        os.system(cmd)
        
        os.system('samtools index %s/%s' %(preDir, re.sub(umi_file_prefix,'UBfix.coordinateSorted_unique.bam',os.path.basename(inputBAM))))
        os.system('samtools index %s/%s' %(preDir, re.sub(umi_file_prefix,'UBfix.coordinateSorted_multi.bam',os.path.basename(inputBAM))))

    if op.reconstruction:
        
        print('Collect informative reads per gene...')
        in_bam_uniq = '%s/%s' %(os.path.join(outdir, species, experiment, "preprocess"), re.sub(umi_file_prefix,'UBfix.coordinateSorted_unique.bam',os.path.basename(inputBAM)))
        in_bam_multi = '%s/%s' %(os.path.join(outdir, species, experiment, "preprocess"), re.sub(umi_file_prefix,'UBfix.coordinateSorted_multi.bam',os.path.basename(inputBAM)))
    
        out_path = os.path.join(outdir, species, experiment, "isoforms_%s" %(conf_data['annotation']['gtf_source']))
        if not os.path.exists(out_path): os.makedirs(out_path)
        
        sys_tmp_dir = '%s/.tmp' %(out_path)
        if not os.path.exists(sys_tmp_dir): os.makedirs(sys_tmp_dir)
        pybedtools.set_tempdir(sys_tmp_dir)
        pybedtools.cleanup(remove_all=True)
        
        fetch_gene_reads(in_bam_uniq, in_bam_multi, conf_data, op.species, out_path)
        
        print('Build reference isoforms...')
        ref = build_reference(conf_data, out_path)
        
        print('Start isoform reconstruction...')
        get_isoforms(conf_data, out_path, ref)
        

if __name__ == '__main__':
        main()

    
