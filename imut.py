import os
import sys
import pysam
import argparse

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('func', choices=['snv','indel'], nargs=1)
    parser.add_argument('bam', nargs=1)
    parser.add_argument('ref', nargs=1)
    parser.add_argument('loc', nargs=1)
    parser.add_argument('--output', nargs='?')

    args = vars(parser.parse_args())

    for k,v in args.items():
        print(k, v)

    '''
    bam_file = sys.argv[1]
    chrm = 'chr17'
    pos = 37872045

    samfile = pysam.AlignmentFile(bam_file, 'rb')
    for ppc in samfile.pileup(chrm, pos-1, pos, truncate=True):
        print(ppc.get_mapping_qualities())
        print(ppc.get_num_aligned())
        print(ppc.get_query_names())
        print(ppc.get_query_positions())
        print(ppc.get_query_sequences())
    '''

if __name__ == '__main__':
    main()
