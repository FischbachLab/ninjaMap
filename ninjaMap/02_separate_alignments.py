#!/usr/bin/env python3
'''
Purpose:
    Split the filtered bamfile into Singular and Escrow alignment bamfiles
'''
import argparse
import logging
import os
import sys
from time import perf_counter as timer

import ninjamap.Utils as nm
# from ninjautils.Reads import Reads
# from ninjautils.Strains import Strains

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
                    help='Filtered and name sorted bam file.')
    p.add_argument('-ftab', dest='filter_table', action='store', type=str, required = True,
                    help='Alignment table resulting from 01_filter_alignments')
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
    # p.add_argument('-min_pid', dest='min_pid', action='store', type=float, default=100,
    #                 help='minimum alignment percent identity')
    # p.add_argument('-min_paln', dest='min_paln', action='store', type=float, default=100,
    #                 help='minimum percent alignment length')
    p.add_argument('-cores', dest='num_cores', action='store', type=int, default=4,
                    help='number of threads')
    p.add_argument('-mem', dest='mem_per_core', action='store', type=int, default=4,
                    help='memory per thread (in GB)')

    return vars(p.parse_args())
# DEBUG
# bamfile_name = "/Users/sunit.jain/Research/SyntheticCommunities/ReadAlignment/Testing/Mismaps/Bacteroides-coprophilus-DSM-18228/Bacteroides-coprophilus-DSM-18228.processed.bam"
# abundance_output_file = 'B_coprophilius.ninjaMap.v1.abundance.tsv'
if __name__ == '__main__':
    start = timer()
    args = usage()
    bamfile_name = args['bamfile']
    binmap_file = args['binmap']
    filtered_table_file = args['filter_table']
    output_dir = args['outdir']
    num_cores = args['num_cores']
    mem_per_core = args['mem_per_core']

    os.makedirs(output_dir, exist_ok=True)

    if not args['prefix']:
        default_prefix = os.path.basename(bamfile_name).split('.')[0]
        prefix = os.path.join(output_dir, default_prefix)
    else:
        prefix = os.path.join(output_dir, args['prefix'])

    paired = True
    if args['single'] or args['ccs']:
        paired = False

    # minPercId = args['min_pid']
    # minPercAln = args['min_paln']
    # if args['ccs']:
    #     minPercId = 100
    #     minPercAln = 100

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

    bins, all_strain_obj = nm.parse_db_metadata(binmap_file)

    logging.info('\t%d contigs assigned to %d strains', len(bins.keys()), len(all_strain_obj.keys()))
    ###############################################################################
    # Parse the filtered alignmant table to get Primary and Escrow Reads
    ###############################################################################
    logging.info('Separating data into Singular and Escrow Data frames ...')
    primary_reads_list, escrow_reads_list = nm.identify_primary_and_escrow_reads(filtered_table_file, 
                                                                        prefix = prefix, paired = paired)

    logging.info('Separating data into Singular and Escrow BAM files ...')
    primary_bam, escrow_bam = nm.split_bam_by_readtype(bamfile_name, bins, primary_reads_list, escrow_reads_list, 
                                                            prefix = prefix, paired = paired, by='coord',
                                                            cores = num_cores, memPerCore = mem_per_core, log = logging)

    # Primary/Singular Coverage
    logging.info('Sorting Singular BAM files by coordinates and indexing ...')
    # Calculate coverages for bamfile, using the function calculate_coverage
    #   cov_df columns {'contig','pos', 'depth'}
    # pandas_write_csv compressed
    logging.info('Calculating Contig coverage using Singular reads ...')
    nm.calculate_coverage(primary_bam, output_dir).to_csv(f'{prefix}.primary_coverage.csv.gz', index=False)

    # Escrow Coverage
    logging.info('Sorting Escrow BAM files by coordinates and indexing ...')
    # sorted_escrow_bam = Utils.sort_and_index(escrow_bam, cores = num_cores, memPerCore = mem_per_core)
    logging.info('Calculating Contig coverage using Escrow reads ...')
    # pandas_write_csv compressed
    nm.calculate_coverage(escrow_bam, output_dir).to_csv(f'{prefix}.escrow_coverage.csv.gz', index=False)
