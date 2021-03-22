#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pysam
import sys
import statistics

from collections import defaultdict

def usage():
    p = argparse.ArgumentParser(   
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True,
        usage=argparse.SUPPRESS,
        description="""Description:
    This script will allow you to filter a bam file based on certain flags.
    """,
        epilog="""Examples:
    python bamParser.py -bam sample.sortedByCoord.bam -id 70 -aln_len 80 -out filtered.bam
    or
    python bamParser.py -bam sample.sortedByCoord.bam -id 70 -aln_len 80 -reads readnames.txt -out filtered.bam
    or
    python bamParser.py -bam sample.sortedByCoord.bam -id 70 -aln_len 80 -reads readnames.txt -out filtered.tsv -out_fmt tsv
    """)
    # Required
    p.add_argument('-bam', dest='bamfile', action='store', type=str, required = True,
                    help='coord sorted bam file.')
    p.add_argument('-out', dest='output', action='store', type=str, required = True,
                    help='output file name')

    p.add_argument('-id', dest='min_id', action='store', type=float,
                    help='minimum percent nucleotide identity (inclusive).', default = 0)
    p.add_argument('-aln_len', dest='min_aln_len', action='store', type=float,
                    help='minimum percent read alignment length (inclusive).', default = 0)
    p.add_argument('-reads', dest='read_list', action='store', type=str,
                    help='only print reads mentioned in this file. 3 tab sep columns w/o header: read_name, organism_name, read_length')
    p.add_argument('-out_fmt', dest='output_fmt', action='store', type=str, default = 'bam',
                    help='output format. Choose either "tsv" or "bam". Default is bam')
    p.add_argument('-plot', dest='plot', action='store_true',
                    help='plot percent id and alignment length histograms', default = False)

    return vars(p.parse_args())

def parse_aln_record(aln):
    edit_dist = dict(aln.tags)['NM']
    query_len = aln.query_length
    ref_start = aln.reference_start
    ref_end = aln.reference_end

    pid = (query_len - edit_dist)*100/query_len
    aln_len = aln.get_overlap(ref_start, ref_end)*100/query_len

    return (pid, aln_len)

def acceptable_alignment(aln_pid, aln_len, min_pid, min_paln):
    # https://www.biostars.org/p/106126/
    # return ((aln_pid >= min_pid) and (aln_len >= min_paln))
    return (min_pid <= aln_pid <= 100) and (min_paln <= aln_len <= 100)

def get_series_stats(given_series):
    given_series = np.array(given_series)
    series_mean = np.mean(given_series)
    series_median = np.median(given_series)
    series_sd = np.std(given_series)
    series_var = np.var(given_series)
    series_coeff_var = series_sd / series_mean

    return (series_mean, series_median, series_sd, series_var, series_coeff_var)

def filter_bam(bamfile_name, min_id, min_aln_len,
              output_file_name, output_fmt, read_file, read_filter, plot):
    # open the BAM file for reading
    bam_in = pysam.AlignmentFile(bamfile_name, mode = 'rb')
    # figure out the output file handle
    if output_fmt == 'tsv':
        out_file = open(output_file_name, 'w')
    else:
        out_file = pysam.AlignmentFile(output_file_name, "wb", template=bam_in)

    aln_kept = 0
    aln_total = 0
    # create the read filter dict
    reads = defaultdict(int)
    if read_filter:
        with open(read_file, 'r') as readFile:
            for line in readFile:
                read, subject, read_bp = line.strip().split('\t')
                reads[read] += 1
    pid_series = []
    alnLen_series = []
    # Go through the bam file, line-by-line till the end of file
    for aln in bam_in.fetch(until_eof=True):
        aln_total += 1
        # Continue until you find a read that's in the read filter
        if read_filter and not reads[aln.query_name]:
            continue
        if not aln.has_tag('NM'):
            continue
        if not aln.query_length:
            continue

        # Write to appropriate output 
        # if alignment matches all the minimum alignment criteria?
        aln_pid, aln_aln_len = parse_aln_record(aln)
        if acceptable_alignment(aln_pid, aln_aln_len, min_id, min_aln_len):
            aln_kept += 1
            if read_filter:
                out_file.write(f'{aln.query_name}\t{aln.reference_name}\n')
            else:
                out_file.write(aln)
            
            if plot:
                pid_series.append(aln_pid)
                alnLen_series.append(aln_aln_len)

    if len(pid_series) > 0:
        plt.hist(pid_series, bins=100)
        plt.xlabel("Percent Identity")
        plt.savefig(f"{output_file_name}.perc_id.png")
        # (pid_mean, pid_median, pid_sd, pid_var, pid_coeff_var) = get_series_stats(pid_series)
        # print(f'{pid_mean}\t{pid_median}\t{pid_sd}\t{pid_var}\t{pid_coeff_var}')

    if len(alnLen_series) > 0:
        plt.hist(alnLen_series, bins=100)
        plt.ylabel("Alignment Length")
        plt.savefig(f"{output_file_name}.alnLen.png")
        # (alnLen_mean, alnLen_median, alnLen_sd, alnLen_var, alnLen_coeff_var) = get_series_stats(alnLen_series)
        # print(f'{alnLen_mean}\t{alnLen_median}\t{alnLen_sd}\t{alnLen_var}\t{alnLen_coeff_var}')

    # close the files
    out_file.close()
    bam_in.close()

    return (aln_kept, aln_total)

if __name__ == "__main__":
    args = usage()

    bamfile_name = args['bamfile']
    output_file_name = args['output']
    output_fmt = args['output_fmt']

    min_id = args['min_id']
    min_aln_len = args['min_aln_len']
    read_filter = False
    read_file = ''
    if args['read_list']:
        read_filter = True
        read_file = args['read_list']
    else:
        output_fmt = 'bam'
    plot = args['plot']

    (aln_remaining, total_aln) = filter_bam(bamfile_name,min_id, min_aln_len,output_file_name,output_fmt,read_file,read_filter, plot)
    print(f'Retained {aln_remaining} out of {total_aln}. [~{round(aln_remaining*100/total_aln, 3)}%]')
    sys.exit(0)