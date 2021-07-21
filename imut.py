#!/Users/ysh/anaconda3/envs/imut/bin/python

import os
import sys
import pysam
import argparse

print(pysam.__version__)

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('func', choices=['snv','indel'], nargs=1, help='"snv" or "indel"')
    parser.add_argument('bam', nargs=1, help='bam file')
    parser.add_argument('ref', nargs=1, help='reference fasta, like: hg19.fa')
    parser.add_argument('loc', nargs=1, help='like: chr1:10000')
    parser.add_argument('--output', nargs='?', help='text file write into')

    args = vars(parser.parse_args())

    for k,v in args.items():
        print(k, v)

    bam_file = args['bam'][0]
    chrm, pos = args['loc'][0].split(':')
    pos = int(pos)

    bam = pysam.AlignmentFile(bam_file, 'rb')
    for ppc in bam.pileup(chrm, pos-2, pos, truncate=True):
        print('ppc.get_mapping_qualities()', ppc.get_mapping_qualities())
        print('ppc.get_num_aligned()',       ppc.get_num_aligned())
        print('ppc.get_query_names()',       ppc.get_query_names())
        print('ppc.get_query_positions()',   ppc.get_query_positions())
        print('ppc.get_query_sequences()',   ppc.get_query_sequences())


if __name__ == '__main__':
    main()
