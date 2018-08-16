#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2018-08-16 18:29

from __future__ import print_function
import argparse
import collections

def combine_fasta(alignment, out):
    combined_dict = collections.defaultdict(str)
    for align_file in alignment:
        with open(align_file) as fin:
            for line in fin:
                if not line.strip():
                    continue
                if line.startswith('>'):
                    species = line.split()[0][1:]
                else:
                    combined_dict[species] += line.strip()
    with open(out, 'w') as fout:
        for species, seq in combined_dict.iteritems():
            fout.write('>{0}\n{1}\n'.format(species, seq))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('alignment', nargs='+',
            help='alignment files for combination')
    parser.add_argument('--out', type=str, 
            default='combined.muscle.fasta',
            help='output combined alignment')
    args = parser.parse_args()
    
    combine_fasta(args.alignment, args.out)


if __name__ == '__main__':
    main()
