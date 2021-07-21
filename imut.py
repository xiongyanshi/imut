#!/Users/ysh/anaconda3/envs/imut/bin/python

import os
import sys
import pysam
import argparse

print(pysam.__version__)

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('func', choices=['snv','indel'], nargs=1,
                          help='"snv" or "indel"')
    parser.add_argument('bam', nargs=1, help='bam file')
    parser.add_argument('ref', nargs=1, help='reference fasta, like: hg19.fa')
    parser.add_argument('loc', nargs=1, help='like: chr1:10000')
    parser.add_argument('--output', nargs='?', help='text file write into')

    args = vars(parser.parse_args())

    for k,v in args.items():
        print(k, v)

    bam_file = args['bam'][0]
    ref = args['ref'][0]
    chrm, pos = args['loc'][0].split(':')
    pos = int(pos)

    bam = pysam.AlignmentFile(bam_file, 'rb', reference_filename=ref)
    for ppc in bam.pileup(chrm, pos-1, pos, truncate=True):
        for base, baseq, basep, mapq, read in zip(ppc.get_query_sequences(mark_matches=True),
                ppc.get_query_qualities(), ppc.get_query_positions(),
                ppc.get_mapping_qualities(), ppc.pileups):
            print(base, baseq, basep, mapq, read)


if __name__ == '__main__':
    main()
