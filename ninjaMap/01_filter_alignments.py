#!/usr/bin/env python3
'''
Purpose:
    filter a BAM file based on minimum alignment criteria
'''
import argparse
import logging
import os
import sys
from time import perf_counter as timer

import ninjamap.Utils as nm

# Stopgap setup until I figure out how to setup a python module properly.
# sys.path.insert(0, '/Users/sunit.jain/GitHub/czbiohub/ninjaMap')
# sys.path.insert(0, '/mnt')

###############################################################################
# Setup Input and Script Usage
###############################################################################
def usage():
    p = argparse.ArgumentParser(   
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True,
        usage=argparse.SUPPRESS,
        description="""Description:
    This script will calculate the abundance of a strain in a defined microbial community. 
    Usage: ninjaMap.py -bam name_sorted.bam -bin contig_strain_assignments.tsv -out abundance_table_output.tsv
    """,
        epilog="""Examples:
    python ninjaMap.py -bin contig_names_bin_map.txt -bam Bacteroides-sp-9-1-42FAA/Bacteroides-sp-9-1-42FAA.processed.sortedByCoord.bam -prefix Bacteroides-sp-9-1-42FAA    
    """)
    # Required
    p.add_argument('-bam', dest='bamfile', action='store', type=str, required = True,
                    help='name sorted bam file.')
    p.add_argument('-fasta', dest='fastafile', action='store', type=str, required = True,
                    help='database fasta file')
    p.add_argument('-bin', dest='binmap', action='store', type=str, required = True,
                    help='tab-delimited file with Col1= contig name and Col2=Bin/Strain name')
    p.add_argument('-outdir', dest='outdir', action='store', type=str, required = True,
                    help='output directory')
    # Optional
    p.add_argument('-prefix', dest='prefix', action='store', type=str,
                    help='output prefix')
    p.add_argument('-ccs', dest='ccs', action='store_true', default=False,
                help='Alignment used Pacbio CCS reads')
    p.add_argument('-se', dest='single', action='store_true', default=False,
                help='Alignment used Single ended illumina reads')
    p.add_argument('-min_pid', dest='min_pid', action='store', type=float, default=100,
                    help='minimum alignment percent identity')
    p.add_argument('-min_paln', dest='min_paln', action='store', type=float, default=100,
                    help='minimum percent alignment length')

    return vars(p.parse_args())
# DEBUG
# bamfile_name = "/Users/sunit.jain/Research/SyntheticCommunities/ReadAlignment/Testing/Mismaps/Bacteroides-coprophilus-DSM-18228/Bacteroides-coprophilus-DSM-18228.processed.bam"
# abundance_output_file = 'B_coprophilius.ninjaMap.v1.abundance.tsv'
if __name__ == '__main__':
    start = timer()
    args = usage()
    bamfile_name = args['bamfile']
    fastafile_name = args['fastafile']
    binmap_file = args['binmap']
    output_dir = args['outdir']

    os.makedirs(output_dir, exist_ok=True)

    if not args['prefix']:
        default_prefix = os.path.basename(bamfile_name).split('.')[0]
        prefix = os.path.join(output_dir, default_prefix)
    else:
        prefix = os.path.join(output_dir, args['prefix'])
    
    paired = True
    if args['single'] or args['ccs']:
        paired = False

    minPercId = args['min_pid']
    minPercAln = args['min_paln']
    if args['ccs']:
        minPercId = 100
        minPercAln = 100

    max_perc_indel = 0.1
    min_f_score =  90

    logging.basicConfig(
        # filename=logfile, 
        # filemode='w+', 
        level=logging.DEBUG,
        format='%(asctime)s\t[%(levelname)s]:\t%(message)s')

    logging.info('Started')
    ###############################################################################
    # Parse the Contig to Bin/Strain name map file.
    ###############################################################################
    logging.info('Processing the Bin Map file: %s ...', binmap_file)

    bins, all_strain_obj = nm.parse_db_metadata(binmap_file, fastafile_name)

    logging.info('\t%d contigs assigned to %d strains', len(bins.keys()), len(all_strain_obj.keys()))
    ###############################################################################
    # Parse the BAM file
    ###############################################################################
    logging.info('Processing the BAM file: %s ...', bamfile_name)

    (reads_w_valid_alignments, 
    total_reads_aligned, 
    total_alignments, 
    valid_alignments, 
    bad_aln) = nm.filter_bam(bamfile_name, 
                            minPercId, minPercAln,
                            max_perc_indel, min_f_score,
                            bins, all_strain_obj, prefix, paired)

    ###############################################################################
    # Run Summary
    ###############################################################################
    if total_reads_aligned == 0:
        logging.critical(f'{total_reads_aligned} reads aligned to the reference database. Please check the BAM file or the database used.',)
        sys.exit(1)

    read_passed_percent = round(reads_w_valid_alignments*100/total_reads_aligned,3)
    aln_passed_percent = round(valid_alignments*100/total_alignments,3)
    logging.info(f'\tFound {total_alignments} alignments for {total_reads_aligned} reads. Of which, {bad_aln} were invalid alignments')
    logging.info(f'\t{valid_alignments} alignments (~{aln_passed_percent}%) passed the acceptable alignment criteria.')
    logging.info(f'\tAcceptable alignments belonged to {reads_w_valid_alignments} reads (~{read_passed_percent}%).')
    ###############################################################################
    # The End
    ###############################################################################
    end = timer()
    logging.info('Completed in %s (d:h:m:s)', nm.human_time(end - start))
