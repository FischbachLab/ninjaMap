#!/usr/bin/env python3

import gzip
import os
import sys
from collections import Counter, defaultdict
import math
import numpy as np
import pandas as pd
# import pybedtools
import pysam
from Bio import SeqIO
from functools import partial, reduce 
import operator

import ninjamap.Tools as tools
from ninjamap.Alignments import Alignments
from ninjamap.Reads import Reads
from ninjamap.Strains import Strains

sys.path.insert(0, '/mnt')

# def read_bam_file(bamfile_name, bins):
#     total_reads = 0
#     read_objects = defaultdict()
#     valid_alignment = defaultdict(lambda: defaultdict(int))
#     discarded_reads = set()

#     # Read the BAM file
#     bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
#     for aln in bamfile.fetch(until_eof=True):
#         read_name, mates_name, read_length, template_length = extract_read_info(aln)
#         if is_valid_alignment(aln):
#             # reads_aligned[read_name] = (mates_name, read_length, template_length)
#             read = Reads(read_name, mates_name, read_length, template_length)
#             read_objects[read_name] = read
#             strain_name = bins[aln.reference_name] # bins[contig_name] --> strain name
#             valid_alignment[read_name][strain_name] += 1
#         else:
#             discarded_reads.add(read_name)
        
#     total_reads = len(discarded_reads) + len(read_objects.keys())
#     bamfile.close()
#     return(total_reads, valid_alignment, read_objects)

def choose_primary_candidate(read, mate):
    if read.mate_has_valid_match or mate.mate_has_valid_match :
        common_strains_list = intersection(read.mapped_strains.keys(), mate.mapped_strains.keys())
        if len(common_strains_list) == 1:
            return common_strains_list[0].name
        else:
            return None

def check_for_fraud(votes_file, discard_pile=1):
    fraud_status = True
    df = pd.read_csv(votes_file, index_col='Read_Name', usecols=['Read_Name', 'cSingular_Vote', 'cEscrow_Vote'])
    df['Total_Votes'] = round((df['cSingular_Vote'] + df['cEscrow_Vote']), 5)
    df = df.drop(['cSingular_Vote',  'cEscrow_Vote'], axis = 1)
    votes_df = df.groupby('Read_Name').sum()
    gt_row, gt_col = votes_df[(votes_df.Total_Votes > 1.001)].shape
    mid_row, mid_col = votes_df[(votes_df.Total_Votes > 0.001) & (votes_df.Total_Votes < 0.999)].shape
    if gt_row == 0 and mid_row == 0:
        fraud_status = False

    return fraud_status

def human_time(time):
    time = abs(time)
    day = time // (24 * 3600)
    time = time % (24 * 3600)
    hour = time // 3600
    time %= 3600
    minutes = time // 60
    time %= 60
    seconds = time
    time_str = format('%02d:%02d:%02d:%02d'%(day,hour,minutes,seconds))
    return time_str

def intersection(list1, list2):
    return list(set(list1) & set(list2))

def get_series_stats(given_series):
    given_series = np.array(given_series)
    series_mean = np.mean(given_series)
    series_median = np.median(given_series)
    series_sd = np.std(given_series)
    series_var = np.var(given_series)
    series_coeff_var = series_sd / series_mean

    return (series_mean, series_median, series_sd, series_var, series_coeff_var)

##################################
# Static methods from ninjaIndex
##################################

# def create_db_metadata(fasta_list, fasta_filename):
#     # formerly called get_db_metadata
#     '''
#     read the fasta files
#     return bins, all_strain_obj, strains_list, concatenated fasta file
#     '''
#     bins = defaultdict()
#     all_strain_obj = defaultdict()

#     with open(fasta_filename, "w") as fasta_output:
#         for input_file in fasta_list:
#             # fasta_sequences = SeqIO.parse(open(input_file),'fasta')
#             tmp_strain_name = os.path.basename(input_file)
#             strain_name = os.path.splitext(tmp_strain_name)[0]
#             strain = Strains(strain_name)
#             all_strain_obj[strain_name] = strain
#             # strains_list.append(strain_name)
#             with open(input_file, "r") as fasta_input:
#                 for fasta in SeqIO.parse(fasta_input,'fasta'):
#                     contig_name, sequence = fasta.id, str(fasta.seq)
#                     contig_len = len(sequence)
#                     strain.add_contig(contig_name, contig_len)
#                     bins[contig_name] = strain_name
#                     Strains.total_genome_size += contig_len

#                     SeqIO.write(fasta, fasta_output , "fasta")

#     return (bins, all_strain_obj, fasta_filename)

# def create_bin_map(all_strain_obj, binmap_file):
#     binmap = open(binmap_file, 'w')
#     # Header
#     # strain_name, strain_wt_1, strain_wt_2, strain_wt_3, strain_score, strain_unique_bases, contig_name, contig_length
#     binmap.write('Strain_Name,Strain_Weight_Andres,Strain_Weight_Original, Strain_Weight_UScore,Strain_Uniqueness_Score,Strain_Absolute_Unique_Bases,Contig_Name,Contig_Length\n')
#     for name, strain in all_strain_obj.items():
#         line = ''
#         first_half = strain.name +','+ \
#             str(strain.adj_primary_wt) +','+ \
#             str(strain.adj_primary_wt_2) +','+ \
#             str(strain.adj_primary_wt_3) +','+ \
#             str(strain.uniqueness_score) +','+ \
#             str(strain.uniquely_covered_bases)
#         for contig_name, contig_len in strain.contigs.items():
#             line += str(first_half)  +','+ \
#                     contig_name +','+ \
#                     str(contig_len) +'\n'
#         binmap.write(line)
#     binmap.close()

##################################
# 01_filter_alignments.py
##################################
def parse_db_metadata(binmap_file, fastafile_name=None):
    """
    Function that parses user inputs to create the db metadata
    Parameters
    ----------
    binmap_file: textIO
        the path to depth file for each sample
    Returns
    -------
    depth_dic_coverage: dict
            dictionary with the coverage per position for each plasmid
    """
    try:
        # to handle the new binmap format from indexed database
        bins, all_strain_obj = Strains.get_indexed_db_metadata(binmap_file)
    except:
        bins, all_strain_obj = Strains.get_db_metadata(binmap_file, fastafile_name)

    return (bins, all_strain_obj)

def filter_bam(bamfile_name, min_id, min_aln_len,
                max_perc_indel, min_f_score,
                bins, all_strain_obj, prefix, paired=False):
    """
    Function to parse a bam file and filter out alignments that do not fit the specified thresholds.
    Creates a new filtered bam file and an accompanying alignment stats file.

    Parameters
    ----------
    bamfile_name:
        name of bam file
    min_id:  
    min_aln_len:  
    max_perc_indel:  
    min_f_score:  
    bins:  
    all_strain_obj:  
    prefix:  
    paired: boolean (default=False)  

    Returns
    ----------
    num_unique_valid_alignments: int  
    num_unique_reads_aligned: int  
    total_alignments: int  
    total_valid_alignments: int  
    total_bad_aln: int  
    """
    # Write intermediate files
    filtered_bamfile_name, filtered_bamfile_handle = tools.get_bam_filehandle(bamfile_name, prefix, 'filtered')
    aln_stats_file = prefix +'.ninjaMap.aln_stats.csv.gz'
    aln_stats = gzip.open(aln_stats_file, mode = 'wt')
    aln_stats.write("read_basename,orientation,read_name,read_length,read_mate,"
                    "contig_name,contig_length,genome_name,genome_len,"
                    "percent_id,perc_aln,edit_distance,alignment_length,matches,"
                    "insertions,deletions,precision,recall,f2_score\n")
    
    # Read the BAM file
    bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb',threads = 8)
    # for each contig
    for contig in bins.keys():
        # for each alignment on the contig;
        # removed until_eof due to pysam trying to keep bam file in memory.
        for alignment in bamfile.fetch(contig=contig):
            aln = Alignments(alignment, paired)
            # Don't bother is the NM tag was not found.
            if not aln.is_aln:
                continue

            if aln.is_valid_alignment(min_perc_id = min_id, min_perc_aln = min_aln_len, max_perc_indel = max_perc_indel, min_f_score = min_f_score):
            # for aln in Alignments.fetch_next_alignment(bamfile_name, filtered_bamfile_handle, paired, max_perc_indel, min_f_score):
                contig_name = aln.reference_name
                strain_name = bins[contig_name] # bins[contig_name] --> strain name
                strain = all_strain_obj[strain_name]
                contig_length = strain.contigs[contig_name]
                genome_size = strain.genome_size

                # write to a new BAM file
                filtered_bamfile_handle.write(alignment)
                
                # write a read-strain stats data frame
                aln_stats.write(f"{aln.read_base},{aln.orientation},{aln.read_name},{aln.read_length},{aln.mate_name},"
                                f"{contig_name},{contig_length},{strain_name},{genome_size},"
                                f"{aln.pid},{aln.p_aln_len},{aln.edit_dist},{aln.aln_len},{aln.match},"
                                f"{aln.insertions},{aln.deletions},{aln.precision},{aln.recall},{aln.f_score}"
                                "\n")

    bamfile.close()
    aln_stats.close()
    filtered_bamfile_handle.close()
    filtered_sorted_bamfile_name = tools.sort_and_index(filtered_bamfile_name)

    return (len(Alignments.reads_w_valid_alignments), len(Alignments.total_reads), Alignments.total_alignments, Alignments.total_valid_alignments, Alignments.total_bad_aln)

##################################
# 02_separate_alignments.py
##################################
def filter_by_readnames(whole_df, subset_list):
    return(whole_df[whole_df.read_basename.isin(subset_list)])

def mate_recruitment(df):
    """
    1) Recruitment: If one read of a paired end read set finds an exact match to one and only one strain,
    while the other read in a pair is unassigned, assume that the unassigned read would have
    aligned to the same strain as it's pair had the reference genome/sequencing been of sufficient quality. An example,
    of why this could happen is when the reference has long stretches of Ns.

    2) Compromise: If both reads in a pair are assigned to multiple strains. Identify a set of strains, 
    that are common between the pair and discard the remaining strain assignments not in the intersection.
    If a common set of strains cannot be found, assign the read pair to the union of the strain assignments.
    """
    intersect = partial(reduce, operator.and_)
    df = df.reset_index()
    # First, create a set of genomes for each read separately. reset_index()
    # Second, calculate the intersection of genomes between fwd and reverse.
    set_df = (df[['read_basename','orientation','genome_name']]
        .groupby(['read_basename','orientation'])
        .agg(
            genomes_set = ('genome_name', lambda x: set(x))
        )
        .reset_index()
        .groupby('read_basename')
        .agg(
            genomes_intersection = ('genomes_set', intersect)
        )
        .reset_index()
    )
    set_df['num_genomes'] = set_df['genomes_intersection'].apply(len)

    return set_df[['read_basename','num_genomes']]

def identify_primary_and_escrow_reads(whole_table_file, prefix, paired):
    df = pd.read_csv(whole_table_file)
    required_columns = ['read_basename', 'orientation','genome_name']
    assert(all(col in df.columns for col in required_columns)),"Missing required column"

    if paired:
        agg_df = mate_recruitment(df)
    else:
        agg_df = (df[['read_basename', 'genome_name']]
                    .groupby("read_basename")
                    .agg(
                        num_genomes = ('genome_name', 'nunique')
                    ).reset_index()
                )

    primary_df = filter_by_readnames(df, agg_df['read_basename'][agg_df.num_genomes == 1])
    escrow_df = filter_by_readnames(df, agg_df['read_basename'][agg_df.num_genomes > 1])

    primary_df.to_csv(f'{prefix}.ninjaMap.primary.csv.gz', index=False)
    escrow_df.to_csv(f'{prefix}.ninjaMap.escrow.csv.gz', index=False)

    primary_reads_series = set(primary_df.read_basename.unique())
    escrow_reads_series = set(escrow_df.read_basename.unique())

    common_reads_error = intersection(primary_reads_series, escrow_reads_series)
    if len(common_reads_error) > 0:
        print(f'[FATAL] {common_reads_error} reads were present in both Primary and Escrow.\nPlease report this incidence to the developer.')
        sys.exit(1)

    return(primary_reads_series, escrow_reads_series)

def write2sam(row, sam_handle):
    # remove empty elements from the list
    # convert each element to string and join with tabs
    clean_row = '\t'.join(str(item) for item in row if item is not None)
    # write to SAM file
    sam_handle.write(f'{clean_row}\n')
    return

# TESTED - WORKS
def split_bam_by_readtype(bamfile_name, bin, primary_list, escrow_list, 
                prefix, paired, log, by='coord', cores=4, memPerCore=4):
    log.info(f'\tConverting BAM file, {bamfile_name} to SAM')
    samfile = tools.bam2sam(bamfile_name, cores = cores)
    
    log.info(f'\tSplitting SAM file, {samfile} to into Primary and Escrow SAM files')
    primary_sam = f'{prefix}.ninjaMap.primary.sam'
    primary_samhandle = open(primary_sam, 'w')

    escrow_sam = f'{prefix}.ninjaMap.escrow.sam'
    escrow_samhandle = open(escrow_sam, 'w')

    # Write headers
    log.info(f'\t\tCreating Primary and Escrow SAM headers from {samfile} ...')
    num_cols = 0
    last_line = 0
    with open(samfile, "r") as samhandle:
        for lnum, line in enumerate(samhandle):
            # columns = line.split('\t')
            # num_cols = len(columns)
            last_line = lnum

            # if num_cols < 11:
            if line.startswith('@'):
                primary_samhandle.write(line)
                escrow_samhandle.write(line)
            else:
                break

            if lnum % 2000 == 0:
                log.info(f'\t\t\tProcessed {lnum} lines as SAM header ...')

    log.info(f'\t\tProcessed a total of {lnum} lines as SAM header')
    # Use pandas to read the dataframe
    # Get SAM num of default columns:
    # May 08,2020: http://samtools.github.io/hts-specs/SAMv1.pdf Sec 1.4
    # 11 mandatory tab-delimited fields, more possible
    # sam_colnames = ['query_name','bit_flag','ref_name', 
    #                 'pos', 'map_qual','cigar','mate_name',
    #                 'mate_pos','template_len','query_seq','query_qual',
    #                 'optional_1', 'optional_2', ... , 'optional_N']
    # https://stackoverflow.com/questions/55129640/read-csv-into-a-dataframe-with-varying-row-lengths-using-pandas
    sam_df = pd.read_csv(samfile, header=None, sep='\n', skiprows = last_line)
    sam_df = sam_df[0].str.split('\t', expand=True)
    num_rows, num_cols = sam_df.shape
    sam_colnames = [f'col_{x + 1}' for x in range(num_cols)]
    sam_df.columns = sam_colnames
    log.info(f'\t\t[TMP] Retained {num_rows} rows and {num_cols} cols from SAM file without headers.')
    
    for row in sam_df.itertuples(index=False):
        if row.col_1 in primary_list:
            # read is primary
            write2sam(row,primary_samhandle)
        elif row.col_1 in escrow_list:
            # read is escrow
            write2sam(row,escrow_samhandle)
        else:
            log.critical(f'\t\tRead could not be found in Primary or Escrow...')
            sys.exit(1)
    primary_samhandle.close()
    escrow_samhandle.close()

    log.info(f'\tConverting Primary and Escrow SAM files to BAM ...')
    primary_bam = tools.sam2bam(primary_sam, sort_by = by, cores = cores)
    escrow_bam = tools.sam2bam(escrow_sam, sort_by = by, cores = cores)
    return (primary_bam, escrow_bam)

def calculate_coverage(bamfile_name, output_dir):
    return tools.calculate_coverage(bamfile_name, output_dir)

##################################
# 03_calculate_abundance.py
##################################

# def get_coverage(binmap_df, # contig mappings to strain names and strain weights from the binmap file
#                 bamfile_name,
#                 prefix, output_dir):
#     # Calculate coverages for bamfile, using the function calculate_coverage
#         # cov_df columns {'contig','pos', 'depth'}
#     sorted_bamfile = sort_and_index(bamfile_name)
#     cov_df = calculate_coverage(sorted_bamfile, output_dir)
    
#     # contigs_df is cov_df subsetted by genomes from the binmap_df

#     # Breadth of coverage > 0.
#     coverage_breadth = contigs_df.pos[contigs_df.depth > 0].count()
#     coverage = coverage_breadth * 100 / self.genome_size
#     depth = contigs_df.depth.sum() / self.genome_size
#         # Pad depth_array with 0 to make it equal to the genome size before calculationf SD and Coeff of Var.
#     depth_array = np.concatenate([contigs_df.depth.to_numpy(),np.zeros(self.genome_size - coverage_breadth)]) * scale
#     depth_sd = np.std(depth_array)
#     depth_variance = depth_sd/np.mean(depth_array) # coeff of variation
#     bases_covered = set(contigs_df.unique_name[contigs_df.depth > 0])
#     (coverage, wt_depth, depth_sd, depth_variance) = strain.calculate_genome_coverage(singular_bed_coverage)
#     # Return Coverage DF
#     return(subset_cov_df)
