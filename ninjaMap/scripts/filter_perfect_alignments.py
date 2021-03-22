#!/usr/bin/env python3

# Each read gets 1 vote for candidate strains. 
# The more strains it maps to the lower the value of it's vote.
# Note: we are not calculating coverage, just read assignments.
# MAJOR ASSUMPTION 1: We will always have a complete genome for organisms (exact strain) that we seek in the sample.

import argparse
import gzip
import logging
import os
import pysam
import pandas as pd
import re
import sys

from collections import defaultdict, Counter
from time import perf_counter as timer

start = timer()
###############################################################################
# Functions for making this more robust
###############################################################################
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

# def bam_is_empty(fn):
#     if os.path.getsize(fn) > 1000000:
#         return False

#     bam = pysam.Samfile(fn, check_sq=False)
#     try:
#         bam.next()
#         return False
#     except StopIteration:
#         return True

# def sort_and_index(file_name, sorted_prefix=None, cores=1):
#     """ Sorts and indexes a bam file by coordinates.
#     """
#     if sorted_prefix is None:
#         sorted_prefix = file_name.replace('.bam', '') + '_sorted'

#     sorted_name = sorted_prefix + '.bam'
#     pysam.sort('-@',cores, file_name, sorted_prefix)
#     pysam.index(sorted_name)

#     return pysam.Samfile(sorted_name, 'rb')

# def sort_by_name(file_name, sorted_name=None, cores=1):
#     """ Sorts a bam file by the read name, for paired-end
#     """
#     if sorted_name is None:
#         sorted_name = file_name.replace('.bam', '') + '_namesorted.bam'

#     pysam.sort('-@',cores,'-n', '-o' , sorted_name, file_name)

#     return pysam.Samfile(sorted_name, 'rb')

###############################################################################
# Side Project
###############################################################################

# def predict_alignment():
#     feature_cols = ["align_len","query_len","aln_cov","quality","perc_id","aln_score","mate_score","mismatches","gap_open","gap_ext","is_dup","is_primary","is_supp"]
#     read the model pickle, predict and write outcome.

###############################################################################
# Setup Input and Script Usage
###############################################################################

usage = """
    USAGE: 
    python ninjaMap.py \
-bam input_bamfile \
-bin tab-delimited file with Col1= contig name and Col2=Bin/Strain name \
-out abundance table output
-log logfile.txt
    """

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
p.add_argument('-debug', dest='debug', action='store_true', default=False,    
                help='save intermediate false positives bam file')
p.add_argument('-truth', dest='truth', action='store', default=False,    
                help='If using debug, please provide one strain name that you would like to track.')
p.add_argument('-mbq', dest='min_base_qual', action='store', default=20, type=int,    
                help='minimum read base quality to consider for coverage calculations.')
p.add_argument('-pacbio', dest='pacbio', action='store_true', default=False,    
                help='Alignment used Pacbio reads')

args = vars(p.parse_args())
# DEBUG
# bamfile_name = "/Users/sunit.jain/Research/SyntheticCommunities/ReadAlignment/Testing/Mismaps/Bacteroides-coprophilus-DSM-18228/Bacteroides-coprophilus-DSM-18228.processed.bam"
# abundance_output_file = 'B_coprophilius.ninjaMap.v1.abundance.tsv'
bamfile_name = args['bamfile']
fastafile_name = args['fastafile']
binmap_file = args['binmap']
output_dir = args['outdir']

os.makedirs(output_dir, exist_ok=True)

minPercAln = 100
minPercId = 100
if args['pacbio']:
    minPercAln = 99
    minPercId = 98

if not args['prefix']:
    default_prefix = os.path.basename(bamfile_name).split('.')[0]
    prefix = os.path.join(output_dir, default_prefix)
else:
    prefix = os.path.join(output_dir, args['prefix'])

abundance_output_file = prefix +'.ninjaMap.abundance.csv'
stats_file = prefix +'.ninjaMap.read_stats.csv'
vote_file = prefix +'.ninjaMap.votes.csv.gz'
fraud_file = prefix +'.ninjaMap.fraud_votes.txt'
strain_stats_file = prefix +'.ninjaMap.strain_stats.csv'
logfile = prefix +'.ninjaMap.log.txt'

if args['min_base_qual']:
    min_base_qual = args['min_base_qual']
else:
    min_base_qual = 20

tmp_bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
true_positives = pysam.AlignmentFile(prefix + '.ninjaMap.true_positives.bam', "wb", template=tmp_bamfile)

logging.basicConfig(
    # filename=logfile, 
    # filemode='w+', 
    level=logging.DEBUG,
    format='%(asctime)s\t[%(levelname)s]:\t%(message)s')

logging.info('Started')
###############################################################################
# Classes
###############################################################################
class Strains:
    total_genome_size = 0
    total_strains = 0
    
    def __init__(self, strain_name):
        self.name = strain_name
        self.num_singular_reads = 0
        self.num_escrow_reads = 0
        self.num_contigs = 0
        self.genome_size = 0
        self.cum_primary_votes = 0
        self.cum_escrow_votes = 0
        self.uniqueness_score = 0
        
        self.total_covered_bases = 0
        self.total_covered_depth = 0
        self.uniquely_covered_bases = 0
        self.uniquely_covered_depth = 0
        self.adj_primary_wt = 0
        self.weighted_base_depth = 0
        
        self.escrow_covered_bases = 0
        self.escrow_covered_depth = 0
        
        self.aln_norm_abundance = 0
        self.genome_norm_abundance = 0
        self.adjusted_votes = 0
        self.escrow_vote_conversion_rate = 0
        self.singular_vote_conversion_rate = 0
        
        self.covered_bases = set()
        self.contigs = defaultdict(int)
        self.singular_bin = defaultdict(int)
        self.escrow_bin = defaultdict(int)

    def __hash__(self):
        return hash(str(self.name))

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not(self == other)
        
    def add_contig(self, contig_name, contig_length):
        self.contigs[contig_name] = contig_length
        self.genome_size += contig_length
        self.num_contigs += 1
        
    def add_paired_singular_vote(self, read_name, mate_name, read_template_len, read_vote = 1):
#         print(str(self) +'\t:\t' + self.name +'\t:\t' + read.unique_name +'\t:\t'+ str(read_vote))
        self.uniquely_covered_bases += read_template_len
        # for R1
        self.singular_bin[read_name] += read_vote
        self.cum_primary_votes += read_vote 
        self.num_singular_reads += 1

        # for R2
        self.singular_bin[mate_name] += read_vote
        self.cum_primary_votes += read_vote
        self.num_singular_reads += 1
        
    
    def add_escrow_vote(self, read_name, read_template_len, read_vote):
        if read_name in self.singular_bin.keys():
            return

        self.escrow_bin[read_name] += read_vote
        self.cum_escrow_votes += read_vote
        self.num_escrow_reads += 1
        self.escrow_covered_bases += read_template_len
    
    def _add_covered_base(self, contig_name, contig_pos):
        unique_contig_name = contig_name+'_'+str(contig_pos)
        self.covered_bases.add(unique_contig_name)
        
    def calculate_singular_coverage (self, bamfile_name, fasta_file):
        if self.num_singular_reads == 0:
            return 0
        
        cov_bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
        fasta = pysam.FastaFile(fasta_file)
        
        for contig_name in self.contigs.keys():
            for pileupcolumn in cov_bamfile.pileup(contig = contig_name, stepper = 'all', min_base_quality = 20):
                base_cov_contribution = 0
                if pileupcolumn.nsegments > 0:
                    base_cov_contribution = self._calc_base_depth(pileupcolumn, self.singular_bin)

                    if base_cov_contribution > 0:
                        # self.uniquely_covered_bases += 1
                        self._add_covered_base(contig_name, pileupcolumn.reference_pos)
                        self.uniquely_covered_depth += base_cov_contribution
        cov_bamfile.close()
        fasta.close()
        
        return self.uniquely_covered_depth/self.genome_size
    
    def calculate_escrow_coverage (self, bamfile_name, fasta_file):
        if self.num_escrow_reads == 0:
            return 0
        
        cov_bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
        fasta = pysam.FastaFile(fasta_file)
        for contig_name in self.contigs.keys():
#             for pileupcolumn in cov_bamfile.pileup(contig = contig_name, stepper = 'samtools', fastafile = fasta, min_base_quality = min_base_qual):
            for pileupcolumn in cov_bamfile.pileup(contig = contig_name, stepper = 'all', min_base_quality = 20):
                # For this base/column on the contig
                # print ("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
                base_cov_contribution = 0
                if pileupcolumn.nsegments > 0:
                    base_cov_contribution = self._calc_base_depth(pileupcolumn, self.escrow_bin)

                    if base_cov_contribution > 0:
                        # self.escrow_covered_bases += 1
                        self._add_covered_base(contig_name, pileupcolumn.reference_pos)
                        self.escrow_covered_depth += base_cov_contribution

        cov_bamfile.close()
        return self.escrow_covered_depth/self.genome_size

    def _calc_base_depth(self, pileupcolumn, bin_dict):
        wt_base_depth = 0
        for pileupread in pileupcolumn.pileups:
            # For all reads aligned to this base/column on the contig
            # print ('\tbase in read %s = %s' %
            #       (pileupread.alignment.query_name,
            #        pileupread.alignment.query_sequence[pileupread.query_position]))
            read_unique_name = Reads.get_unique_read_name(pileupread.alignment)
            if read_unique_name in bin_dict.keys():
                # Only calculate coverage from reads with perfect alignment(singular or escrow).
                wt_base_depth += bin_dict[read_unique_name]
                self.total_covered_depth += 1
                self.weighted_base_depth += wt_base_depth

        return wt_base_depth
    
    def _normalize_votes(self):
        # Read to Vote conversion ratio
        self.total_covered_bases = len(self.covered_bases)
        
        if self.num_singular_reads > 0:
            self.singular_vote_conversion_rate = self.cum_primary_votes/self.num_singular_reads
        else:
            self.singular_vote_conversion_rate = 0
            
        if self.num_escrow_reads > 0:
            self.escrow_vote_conversion_rate = self.cum_escrow_votes/self.num_escrow_reads
        else:
            self.escrow_vote_conversion_rate = 0

    def compile_general_stats(self):
        '''
        for each object of this class, return a pandas data frame with 
        strains as rows and number of singular and escrow votes
        '''
        if (self.uniquely_covered_bases == 0) and (self.escrow_covered_bases == 0):
            return None
        
        singular_depth = 0
        escrow_depth = 0
        # self.total_covered_bases = len(self.covered_bases)
        self._normalize_votes()
        
        if self.uniquely_covered_bases > 0:
            singular_depth = (self.uniquely_covered_depth/self.uniquely_covered_bases)
            
        if self.escrow_covered_bases > 0:
            escrow_depth = (self.escrow_covered_depth/self.escrow_covered_bases)

        frac_escrow_reads = 0
        if Reads.total_escrow_reads > 0:
            frac_escrow_reads = self.num_escrow_reads/Reads.total_escrow_reads

        frac_singular_reads = 0
        if Reads.total_singular_reads > 0:
            frac_singular_reads = self.num_singular_reads/Reads.total_singular_reads
        
        # Dataframe
        return pd.DataFrame(
            index = [self.name],
            data  = {
                'Genome_Size' : self.genome_size,
                'Percent_Coverage' : self.total_covered_bases * 100 / self.genome_size,
                'Total_Bases_Covered' : self.total_covered_bases,
                'Coverage_Depth' : self.total_covered_depth/self.genome_size,
                'Read_Fraction' : self.read_fraction,
                'Singular_Strain_Weight' : self.adj_primary_wt,
                'Total_Singular_Reads' : self.num_singular_reads,
                'Total_Singular_Votes' : self.cum_primary_votes,
                'Singular_Read_Vote_Ratio' : self.singular_vote_conversion_rate,
                'Singular_Fraction_of_Singular_Reads' : frac_singular_reads,
                'Singular_Coverage' : self.uniquely_covered_bases,
                'Singular_Depth' : singular_depth,
                'Total_Escrow_Reads' : self.num_escrow_reads,
                'Total_Escrow_Votes' : self.cum_escrow_votes,
                'Escrow_Read_Vote_Ratio' : self.escrow_vote_conversion_rate,
                'Fraction_of_all_Escrow_Reads' : frac_escrow_reads,
                'Escrowed_Cov' : self.escrow_covered_bases,
                'Escrowed_Depth' : escrow_depth
            }
        )
    
    def compile_by_abundance(self):
        '''
        return a pandas data frame with 1 row x 4 columns. 
                    relative_abundance, percent_coverage, wt_coverage_depth
        strain_name
        '''
        if self.read_fraction == 0:
            return None

        percent_coverage = self.total_covered_bases * 100 / self.genome_size
        coverage_depth = self.total_covered_depth/self.genome_size

        wt_coverage_depth = self.weighted_base_depth / self.genome_size

        return pd.DataFrame(
                    index = [self.name],
                    data  = {
                        'Read_Fraction' : self.read_fraction,
                        'Percent_Coverage' : percent_coverage,
                        'Coverage_Depth' : coverage_depth,
                        'Weighted_Coverage_Depth' : wt_coverage_depth
                        }
                    )
    
    # def calculate_adj_primary_wt_uniqueness_score(self):
    #     # Sunit's interpretation of Mike drop adjustment
    #     self.adj_primary_wt = (1 - self.uniqueness_score) * self.cum_primary_votes / Reads.total_reads_aligned
    #     return self.adj_primary_wt

    # def calculate_andres_adjustment(self):
    #     # Andres's interpretation of Mike drop adjustment.
    #     # Actually in his interpretation, 'unique # bases' needed to be calculated for all the strains
    #     # a read is assigned to => Size of unique region in a genome compared to other genomes that this read maps.
    #     # This is a VERY costly operation to repeat millions of times, since it cannot be pre-calculated.
    #     # I have adjusted it to use the genome's absolute unique region compared to all genomes in the database.
    #     self.adj_primary_wt = self.cum_primary_votes / self.indexed_uniquely_covered_bases
    #     return self.adj_primary_wt

    def calculate_sunits_original_adjustment(self):
        self.adj_primary_wt = self.cum_primary_votes / Reads.total_reads_aligned
        return self.adj_primary_wt
    
    # def calculate_sunits_simplified_adjustment(self):
    #     self.adj_primary_wt = self.cum_primary_votes / self.uniquely_covered_bases
    #     return self.adj_primary_wt

    # def calculate_mike_drop_penalty(self):
    #     self.uniqueness_score = self.indexed_uniquely_covered_bases / Strains.total_indexed_uniquely_covered_bases
    #     return self.uniqueness_score

    def calculate_read_fraction(self):
        self.read_fraction = (self.cum_escrow_votes + self.cum_primary_votes) * 100 / Reads.total_reads_aligned
        return self.read_fraction

class Reads:
    total_reads_aligned = 0
    reads_w_perfect_alignments = 0
    total_singular_reads = 0
    total_escrow_reads_kept = 0
    total_escrow_reads_discarded = 0
    total_escrow_reads = 0
    total_singular_reads_in_pairs = 0

    def __init__(self, name, mate_name, read_length, template_length):
        self.name = name
        self.unique_name = name
        self.mates_unique_name = mate_name
        self.read_length = abs(read_length)
        self.template_length = abs(template_length)

        self.cum_vote = 0
        self.has_voted = False
        self.in_singular_bin = False
        self.mate_has_perfect_match = False

        self.mapped_strains = defaultdict()

    def __hash__(self):
        return hash(str(self.unique_name))

    def __eq__(self, other):
        return self.unique_name == other.unique_name

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not(self == other)

    def add_exact_match(self, strain):
        self.mapped_strains[strain] = strain.name
    
    def put_pair_in_singular_bin(self, mate):
        self.in_singular_bin = True
        mate.in_singular_bin = True

    def add_vote(self, vote_value):
        # vote_value = round(vote_value, 7)
        self.cum_vote += vote_value
        self.has_voted = True

    def is_fraud(self):
        '''
        for each object of this class, return True if cumulative votes > 1
        '''
        fraud = True
        # approx 1
        if (self.cum_vote < 1.001) and (self.cum_vote > 0.999):
            fraud = False

        # approx 0
        if (self.cum_vote < 0.001):
            fraud = False

        return fraud

    def get_voting_details(self, approved_strain_list):
        '''
        returns a list of 3 element lists, each containing: strain_name, singular vote value and escrow vote value
        '''
        vote_list = list()
        for strain in approved_strain_list:
            strain_name = ''
            escrow_votes = 0
            singular_votes = 0
            cumulative_vote = 0
            
            if self.cum_vote is not None:
                cumulative_vote = self.cum_vote

            if strain.name is not None:
                strain_name = strain.name

            if self.unique_name in strain.escrow_bin.keys():
                escrow_votes = strain.escrow_bin[self.unique_name]

            if self.unique_name in strain.singular_bin.keys():
                singular_votes = strain.singular_bin[self.unique_name]

            vote_list.append([strain_name, singular_votes, escrow_votes, cumulative_vote])
        return vote_list
    
    @staticmethod
    def choose_primary_candidate(read, mate):
        # if (read.template_length == 0) or (mate.template_length == 0):
        #     return None

        if read.mate_has_perfect_match or mate.mate_has_perfect_match :
            common_strains_list = intersection(read.mapped_strains.keys(), mate.mapped_strains.keys())
            if len(common_strains_list) == 1:
                return common_strains_list[0].name
            else:
                return None
        #     if read.in_singular_bin and mate.in_singular_bin:
        #         # they both match a single strain
        #         read_strain = list(read.mapped_strains.keys())[0] # basically the only one.
        #         mate_strain = list(mate.mapped_strains.keys())[0] # basically the only one.
        #         if read_strain.name == mate_strain.name:
        #             # the strains are the same
        #             return read_strain.name
        #     elif read.in_singular_bin and not mate.in_singular_bin:
        #         # R1 matches a single strain, but R2 matches multiple
        #         read_strain = list(read.mapped_strains.keys())[0] # basically the only one.
        #         if read_strain.name in mate.mapped_strains.values():
        #             # If there is an overlap between strain matches, R1 recruits R2 for it's strain match
        #             return read_strain.name
        #     elif mate.in_singular_bin and not read.in_singular_bin:
        #         # R2 matches a single strain, but R1 matches multiple
        #         mate_strain = list(mate.mapped_strains.keys())[0] # basically the only one.
        #         if mate_strain.name in read.mapped_strains.values():
        #             # If there is an overlap between strain matches, R2 recruits R1 for it's strain match
        #             return mate_strain.name
        # return None

    @staticmethod
    def is_perfect_alignment(aln, min_perc_id = 100, min_perc_aln = 100):
        if not aln.has_tag('NM'):
            return None
        
        edit_dist = aln.get_tag('NM')
        query_len = aln.query_length
        ref_start = aln.reference_start
        ref_end = aln.reference_end
        
        if not query_len:
            return None

        pid = (query_len - edit_dist)*100/query_len
        aln_len = aln.get_overlap(ref_start, ref_end)*100/query_len

        # https://www.biostars.org/p/106126/
        # return ((edit_dist == 0) and (query_len == aln.get_overlap(ref_start, ref_end)))
        return (min_perc_id <= pid <= 100) and (min_perc_aln <= aln_len <= 100)
    
    @staticmethod
    def parse_read_name(aln):
        '''
        Accept: AlignmentFile object from PySam
        if read name has a '/', this is the old format. 
        strip the content after the '/', return remaining
        else, return it as is.
        '''
        try:
            key, value = aln.query_name.split("/")
        except ValueError:
            return str(aln.query_name)
        else:
            return str(key)
        
    @staticmethod
    def get_unique_read_name(aln):
        orientation = ''
        if aln.is_read1:
            orientation =  'fwd'
        else:
            orientation =  'rev'
            
        return Reads.parse_read_name(aln) +'__'+ orientation
    
    @staticmethod
    def get_unique_mate_name(aln):
        orientation = ''
        if aln.is_read1:
            orientation =  'rev'
        else:
            orientation =  'fwd'
            
        return Reads.parse_read_name(aln) +'__'+ orientation
    
    @staticmethod
    def extract_read_info(aln):
        read_name = Reads.get_unique_read_name(aln)
        mate_name = Reads.get_unique_mate_name(aln)
        read_length = aln.reference_length
        template_length = aln.template_length

        return (read_name, mate_name, read_length, template_length)
    
    @staticmethod
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

###############################################################################
# Parse the Contig to Bin/Strain name map file.
###############################################################################

logging.info('Processing the Bin Map file: %s ...', binmap_file)
# Header = ('Strain_Name,Contig_Name,Contig_Length\n')
all_strains = defaultdict(list)
with open(binmap_file, "r") as binmap:
    for line in binmap:
        line=line.rstrip()
        contig_name, strain_name = line.split('\t')
        # all_strains[strain_name].append(line)
        all_strains[strain_name].append(contig_name)

all_strain_obj = defaultdict(int)
strains_list = sorted(list(all_strains.keys()), key=str.lower)
fasta = pysam.FastaFile(fastafile_name)
bins = defaultdict()
for strain_name in strains_list:
    strain = Strains(strain_name)
    all_strain_obj[strain_name] = strain
    for contig_name in all_strains[strain_name]:
    # for line in all_strains[strain_name]:
        # binmap_strain_name, strain.uniqueness_score, contig_name, contig_length = line.split('\t')
        strain.add_contig(contig_name, fasta.get_reference_length(contig_name))
        bins[contig_name] = strain_name
        Strains.total_genome_size += fasta.get_reference_length(contig_name)

Strains.total_strains = len(all_strains.keys())
logging.info('\t%d contigs assigned to %d strains', len(bins.keys()), Strains.total_strains)
fasta.close()

del all_strains
###############################################################################
# Parse the BAM file
###############################################################################
logging.info('Processing the BAM file: %s ...', bamfile_name)

total_reads = set()
perfect_alignment = defaultdict(lambda: defaultdict(list))
read_info = defaultdict()

# Read the BAM file
bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
for aln in bamfile.fetch(until_eof=True):
    # read_name = Reads.get_unique_read_name(aln)
    # mates_name = Reads.get_unique_mate_name(aln)
    read_name, mates_name, read_length, template_length = Reads.extract_read_info(aln)
    read_info[read_name] = (mates_name, read_length, template_length)

    if Reads.is_perfect_alignment(aln, min_perc_id = minPercId, min_perc_aln = minPercAln):
        true_positives.write(aln)
        strain_name = bins[aln.reference_name] # bins[contig_name] --> strain name
        perfect_alignment[read_name][strain_name].append(aln)

    total_reads.add(read_name)

bamfile.close()

Reads.total_reads_aligned = len(total_reads)
if Reads.total_reads_aligned == 0:
    logging.critical(f'0 reads aligned to the reference database. Please check the BAM file or the database used.',)
    sys.exit(1)

Reads.reads_w_perfect_alignments = len(perfect_alignment.keys())
logging.info('\tUsed %d reads with perfect alignments, out of %d (%7.3f%%).', 
    Reads.reads_w_perfect_alignments,
    Reads.total_reads_aligned,
    Reads.reads_w_perfect_alignments*100/Reads.total_reads_aligned
    )
del total_reads

tmp_bamfile.close()
true_positives.close()

del perfect_alignment
del read_info