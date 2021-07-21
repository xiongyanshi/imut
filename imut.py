#!/usr/bin/env python

import os
import sys
import re
import pysam
import argparse
import subprocess as sps

#print(pysam.__version__)


def find_mark(n_line, b_line, target):
    '''
    find the index to place a ^ mark.
    '''
    match = re.search(r'\d+', n_line)
    n1 = int(match.group(0))       # number start the line.
    n1_index = match.start(0)

    n = n1
    index = n1_index
    for b in b_line[n1_index+1:]:
        index += 1
        if b in "acgtnACGTN":      # skip '*' in reference line.
            n += 1
        if n == target:
            break

    #print(' '*index + '^')         # 0-based index.
    #print(index)
    return index

def line_core(line, index):
    '''
    keep only target index position overlaped reads
    '''
    l = index
    r = index
    while l > 0:
        l -= 1
        if line[l].isspace():
            break
    while r < len(line):
        if line[r].isspace():
            break
        r += 1
    #print(l, r)
    return ' '*l + line[l:r] + ' '*(len(line)-r)

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('func', choices=['snv','indel'], nargs=1,
                          help='"snv" or "indel"')
    parser.add_argument('bam', nargs=1, help='bam file')
    parser.add_argument('ref', nargs=1, help='reference fasta, like: hg19.fa')
    parser.add_argument('loc', nargs=1, help='like: chr1:10000')
    parser.add_argument('--samtools', nargs='?', help='samtools bin path')
    parser.add_argument('--output', nargs='?', help='text file write into')

    args = vars(parser.parse_args())

    bam = args['bam'][0]
    ref = args['ref'][0]
    chrm, pos = args['loc'][0].split(':')
    pos = int(pos)

    if args['samtools']:
        SAMTOOLS = args['samtools'][0]
    else:
        SAMTOOLS = sps.run('which samtools',
                            shell=True,
                            capture_output=True).stdout.decode('utf-8').strip()


    cmd = 'export COLUMNS=201; {samtools_bin} tview -d T -p {chrm}:{start} {bam} {ref}'.format(
            samtools_bin=SAMTOOLS, chrm=chrm, start=pos-100, bam=bam, ref=ref)
    run = sps.run(cmd, shell=True, capture_output=True)
    run.check_returncode()
    tv_lines = run.stdout.decode('utf-8').strip('\n').split('\n')

    mark_index = find_mark(tv_lines[0], tv_lines[1], pos)
    print(tv_lines[0])
    print(tv_lines[1])
    print(tv_lines[2])
    print(' '*mark_index + '^')
    base_ref = tv_lines[1][mark_index]
    for line in tv_lines[3:]:
        base_query = line[mark_index]
        if base_query in '.,' or base_query == base_ref:
            continue
        line = line_core(line, mark_index)
        if re.match('^\s+$', line):
            continue
        print(line)


if __name__ == '__main__':
    main()
