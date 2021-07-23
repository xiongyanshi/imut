#!/usr/bin/env python3

import os
import sys
import re
import argparse
import subprocess as sps


def find_mark(n_line, b_line, target):
    '''
    find the index to place a ^ mark, 0-based.
    input:  number line, reference bases line and interested position
    '''
    match = re.search(r'\d+', n_line)
    n1 = int(match.group(0))        # number start the line.
    n1_index = match.start(0)

    n = n1
    index = n1_index
    for b in b_line[n1_index+1:]:
        index += 1
        if b in "acgtnACGTN":       # skip '*' in reference line.
            n += 1
        if n == target:
            break
    #print(' '*index + '^')         # 0-based index.

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


def snv(tv_lines, mark_index, view_all):

    print(tv_lines[0])
    print(tv_lines[1])
    print(tv_lines[2])
    print(' '*mark_index + '^')
    base_ref = tv_lines[1][mark_index]
    base_n = {'A':0, 'C':0, 'G':0, 'T':0,
              'a':0, 'c':0, 'g':0, 't':0,
              'N':0, '*':0}
    for line in tv_lines[3:]:
        base_query = line[mark_index]

        if base_query == ' ':
            continue

        if base_query == '.':
            base_real = base_ref.upper()
        elif base_query == ',':
            base_real = base_ref.lower()
        else:
            base_real = base_query

        base_n[base_real] += 1

        if (not view_all) and (base_query in '.,' or base_query == base_ref):
            continue

        line = line_core(line, mark_index)
        print(line)

    template = '{} {} {} {} {}: {:4d} {:4d} {:4d} {:4d} {:4d}\n' + \
               '{} {} {} {} {}: {:4d} {:4d} {:4d} {:4d} {:4d}'
    print(template.format(
           'A','C','G','T','N',
           base_n['A'], base_n['C'], base_n['G'], base_n['T'], base_n['N'],
           'a','c','g','t','*',
           base_n['a'], base_n['c'], base_n['g'], base_n['t'], base_n['*'],
           ))

    return 0


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('func', choices=['snv','indel'], nargs=1,
                        help='"snv" or "indel"')
    parser.add_argument('bam', nargs=1, help='bam file')
    parser.add_argument('ref', nargs=1, help='reference fasta, like: hg19.fa')
    parser.add_argument('loc', nargs=1, help='like: chr1:10000')
    parser.add_argument('-a', '--all', action='store_true',
                        help='enable to view alt (default) and ref reads')
    parser.add_argument('-s', '--samtools', nargs='?', help='samtools bin')
    parser.add_argument('--output', nargs='?', help='text file write into')

    args = vars(parser.parse_args())

    func = args['func'][0]
    bam = args['bam'][0]
    ref = args['ref'][0]
    chrm, pos = args['loc'][0].split(':')
    pos = int(pos)
    view_all = args['all']    # default = False

    if args['samtools']:
        SAMTOOLS = args['samtools'][0]
    else:
        which_samtools = sps.run('which samtools', shell=True, stdout=sps.PIPE)
        which_samtools.check_returncode()
        SAMTOOLS = which_samtools.stdout.decode('utf-8').strip()


    cmd = 'export COLUMNS=201; '
    cmd += '{samtools} tview -d T -p {chrm}:{start} {bam} {ref}'.format(
           samtools=SAMTOOLS, chrm=chrm, start=pos-100, bam=bam, ref=ref)
    run = sps.run(cmd, shell=True, stdout=sps.PIPE)
    run.check_returncode()
    tv_lines = run.stdout.decode('utf-8').strip('\n').split('\n')
    mark_index = find_mark(tv_lines[0], tv_lines[1], pos)

    if func == 'snv':
        snv(tv_lines, mark_index, view_all)
    elif func == 'indel':
        indel(tv_lines, mark_index, view_all)


if __name__ == '__main__':
    main()
