#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2018-08-15 14:58

from __future__ import print_function
import os
import argparse
import ConfigParser
import glob
import collections
import re


def parse_list(list_file, method):
    """parse input list of busco dirs
    
    example:
    Arabidopsis_thaliana	single_copy	/path/to/BUSCO_ath
    Arabidopsis_lyrata	single_copy	/path/to/BUSCO_aly

    the second column could be single_copy, duplicated or both
    """
    with open(list_file) as f:
        dir_dict = collections.OrderedDict()
        for line in f:
            ls = line.strip()
            if not ls:
                continue
            species, select, busco_dir = ls.split()
            if method == 'combined' and select != 'single_copy':
                raise Exception(
                        '"--method combined" only support "single_copy"')
            if species in dir_dict:
                raise Exception('Same species in busco list!')
            dir_dict[species] = (busco_dir, select)
    return dir_dict


def parse_full_table(dir_dict):
    gene_dict = collections.defaultdict(dict)
    for species, (busco_dir, select) in dir_dict.iteritems():
        full_table_list = glob.glob(
                '{0}/run_*/full_table_*.tsv'.format(busco_dir))
        if len(full_table_list) != 1:
            raise Exception('No or more than one full table')
        full_table = full_table_list[0]
        with open(full_table) as f:
            for line in f:
                ls = line.strip()
                if not line.startswith('#') and ls:
                    if select == 'duplicated':
                        filter_set = {'Duplicated'}
                    elif select == 'single_copy':
                        filter_set = {'Complete'}
                    elif select == 'both':
                        filter_set = {'Complete', 'Duplicated'}
                    lss = ls.split('\t')
                    busco_id, status = lss[0:2]
                    if status not in filter_set:
                        continue
                    contig, start, end = lss[2:5]
                    if busco_id in gene_dict[species]:
                        gene_dict[species][busco_id].add(
                                '{0}:{1}-{2}'.format(contig, start, end))
                    else:
                        gene_dict[species][busco_id] = {
                                '{0}:{1}-{2}'.format(contig, start, end)}
    first = True
    for ids in gene_dict.itervalues():
        if first:
            common_ids = set(ids.keys())
            first = False
        common_ids &= set(ids.keys())
    
    common_gene_dict = collections.defaultdict(dict)
    for species, busco_id_dict in gene_dict.iteritems():
        for busco_id, genes in busco_id_dict.iteritems():
            if busco_id in common_ids:
                common_gene_dict[species][busco_id] = genes
    return common_gene_dict


def parse_extracted_seqs(dir_dict, common_gene_dict):
    seq_dict = collections.defaultdict(dict)
    for species, (busco_dir, select) in dir_dict.iteritems():
        for busco_id, genes in common_gene_dict[species].iteritems():
            faa_list = glob.glob(
                    '{0}/run_*/augustus_output/extracted_proteins/'
                    '{1}.faa.*'.format(busco_dir, busco_id))
            n = 0
            for faa_file in faa_list:
                with open(faa_file) as f:
                    for line in f:
                        match = re.match(r'>[\w]+\[(.+)\]\n', line)
                        if match:
                            ID = match.group(1)
                            if ID in common_gene_dict[species][busco_id]:
                                n += 1
                                seq_dict[busco_id][
                                        '{0}_{1} [{2}]'.format(species, n, ID)] = ''
                                keep = True
                            else:
                                keep = False
                        elif keep:
                            seq_dict[busco_id][
                                    '{0}_{1} [{2}]'.format(species, n, ID)] += line.strip()
    return seq_dict


def output_fasta(seq_dict):
    os.mkdir('00.seq')
    for busco_id, genes in seq_dict.iteritems():
        with open('00.seq/{0}.fasta'.format(busco_id), 'w') as f:
            for ID, seq in genes.iteritems():
                f.write('>{0}\n{1}\n'.format(ID, seq))


def run_muscle(seq_dict, muscle_path):
    os.mkdir('01.muscle')
    with open('01.muscle/work.sh', 'w') as f:
        for busco_id in seq_dict:
            f.write('{0} -in ../00.seq/{1}.fasta -out '
                    './{1}.muscle.fasta\n'.format(muscle_path, busco_id))


def run_fasttree(seq_dict, method, fasttree_path, combine_fasta_path):
    os.mkdir('02.fasttree')
    if method == 'combined':
        with open('02.fasttree/work.sh', 'w') as f:
            f.write('{0} {1} --out {2}\n'.format(
                combine_fasta_path, 
                ' '.join(['../01.muscle/{0}.muscle.fasta'.format(
                    busco_id) for busco_id in seq_dict]), 
                'combined.muscle.fasta'))
            f.write('{0} combined.muscle.fasta > '
                    './combined.fasttree.nwk\n'.format(fasttree_path))
    else:
        with open('02.fasttree/work.sh', 'w') as f:
            for busco_id in seq_dict:
                f.write('{0} ../01.muscle/{1}.muscle.fasta > '
                        './{1}.fasttree.nwk\n'.format(fasttree_path, busco_id))

def workflow(mode, method, qsub_sge_path):
    with open('work.sh', 'w') as f:
        if mode == 'sge':
            f.write('cd 01.muscle && {0} work.sh '
                    '--lines 10 --cpu 50 --threads 1 '
                    '--vf 1G --convert && cd ..\n'.format(qsub_sge_path))
            if method == 'separate':
                f.write('cd 02.fasttree && {0} work.sh '
                        '--lines 10 --cpu 50 --threads 1 '
                        '--vf 1G --convert && cd ..\n'.format(qsub_sge_path))
            else:
                f.write('cd 02.fasttree && {0} work.sh '
                        '--lines 2 --cpu 1 --threads 1 '
                        '--vf 1G --convert && cd ..\n'.format(qsub_sge_path))
        else:
            f.write('cd 01.muscle && sh work.sh && cd ..\n')
            f.write('cd 02.fasttree && sh work.sh\n')

def parse_config(cfg_file):
    config = ConfigParser.ConfigParser()
    config.read(cfg_file)
    muscle_path = config.get('software', 'muscle')
    fasttree_path = config.get('software', 'fasttree')
    qsub_sge_path = config.get('software', 'qsub_sge')
    return muscle_path, fasttree_path, qsub_sge_path


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('busco_list',
            help="input a list of BUSCO dirs")
    parser.add_argument('--method', 
            choices={'combined', 'separate'}, default='combined',
            help='construct phylogenetic tree using gene supermatrix '
            '(combined, only support "single_copy" mode) or using each '
            'gene separately (separate) [default: %(default)s]')
    parser.add_argument('--mode', 
            choices={'sge', 'local'}, default='sge',
            help='select running mode')
    return parser.parse_args()


def main():
    args = parse_arguments()
    real_path = os.path.split(os.path.realpath(__file__))[0]
    combine_fasta_path = real_path + '/combine_fasta.py'
    cfg_file = real_path + '/busco_tree.cfg'
    muscle_path, fasttree_path, qsub_sge_path = parse_config(cfg_file)
    
    dir_dict = parse_list(args.busco_list, args.method)
    common_gene_dict = parse_full_table(dir_dict)
    seq_dict = parse_extracted_seqs(dir_dict, common_gene_dict)
    output_fasta(seq_dict)
    run_muscle(seq_dict, muscle_path)
    run_fasttree(seq_dict, args.method, fasttree_path, combine_fasta_path)
    workflow(args.mode, args.method, qsub_sge_path)

if __name__ == '__main__':
    main()
