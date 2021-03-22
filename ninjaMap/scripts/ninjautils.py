#!/usr/bin/env python3

import argparse
import gzip
import logging
import numpy as np
import os
import pybedtools
import pysam
import pandas as pd
import re
import sys

from collections import defaultdict, Counter
from time import perf_counter as timer

class Strains:
    total_genome_size = 0
    total_strains = 0
    total_uniquely_covered_bases = 0
    net_promiscuity = 0
    
    def __init__(self, strain_name):
        self.name = strain_name
        self.num_singular_reads = 0
        self.num_escrow_reads = 0
        self.num_contigs = 0
        self.genome_size = 0
        self.cum_primary_votes = 0
        self.cum_escrow_votes = 0
        self.uniqueness_score = 0
        self.percent_coverage = 0
        self.depth_variance = 0
        
        self.breadth_of_coverage = 0
        self.depth_of_coverage = 0

        self.uniquely_covered_bases = set()
        self.num_uniquely_covered_bases = 0
        self.uniquely_covered_depth = 0
        self.adj_primary_wt = 0
        self.singular_depth_var = 0
        
        self.escrow_covered_bases = set()
        self.num_escrow_covered_bases = 0
        self.escrow_covered_depth = 0
        self.escrow_depth_var = 0
        
        self.aln_norm_abundance = 0
        self.genome_norm_abundance = 0
        self.adjusted_votes = 0
        self.escrow_vote_conversion_rate = 0
        self.singular_vote_conversion_rate = 0
        
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
        # self.uniquely_covered_bases += read_template_len
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
        # self.escrow_covered_bases += read_template_len

    def calculate_genome_coverage(self, df, singular=True):
        if self.num_singular_reads == 0:
            return (0,0,0,0)

        coverage, depth, depth_sd, depth_variance = (0,0,0,0)

        # Scaling factor for escrow depth values
        scale = 1
        if not singular:
            scale = self.adj_primary_wt / Strains.net_promiscuity
            self.escrow_covered_bases = set()
            self.escrow_depth = 0
            self.num_escrow_covered_bases = 0
            self.escrow_depth_var = 0

        # filter by contig names for this strain
        contigs_df = df[df.contig.isin(self.contigs.keys())].copy()
        num_bases = contigs_df.shape[0]
        if num_bases > 0:
            contigs_df['unique_name'] = contigs_df['contig'] + '_' + contigs_df['pos'].astype(str)
            # contigs_df['unique_name'] = contigs_df[['contig', 'pos']].apply(lambda x: '_'.join(str(x)), axis = 1)
            # print(f'{contigs_df.shape}')
            # coverage = contigs_df.num_bases[contigs_df.depth > 0].sum() * 100 / self.genome_size
            # wt_depth = (contigs_df.depth * contigs_df.num_bases).sum() / self.genome_size
            # depth_variance = (contigs_df.depth * contigs_df.num_bases).var() * 100 / self.genome_size
            
            # Breadth of coverage > 0.
            coverage_breadth = contigs_df.pos[contigs_df.depth > 0].count()
            coverage = coverage_breadth * 100 / self.genome_size
            depth = contigs_df.depth.sum() / self.genome_size
            # Pad depth_array with 0 to make it equal to the genome size before calculationf SD and Coeff of Var.
            depth_array = np.concatenate([contigs_df.depth.to_numpy(),np.zeros(self.genome_size - coverage_breadth)]) * scale
            depth_sd = np.std(depth_array)
            depth_variance = depth_sd/np.mean(depth_array) # coeff of variation
            bases_covered = set(contigs_df.unique_name[contigs_df.depth > 0])
            if singular:
                self.uniquely_covered_bases = bases_covered
                self.singular_depth = depth
                self.num_uniquely_covered_bases = coverage_breadth
                self.singular_depth_var = depth_variance
                Strains.total_uniquely_covered_bases += coverage_breadth
            else:
                self.escrow_covered_bases = bases_covered
                self.escrow_depth = depth
                self.num_escrow_covered_bases = coverage_breadth
                self.escrow_depth_var = depth_variance
        
        # Allow this temporary dataframe to be garbage collected.
        del contigs_df

        return (coverage, depth, depth_sd, depth_variance)

    def _normalize_votes(self):
        # Read to Vote conversion ratio
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
        if (self.num_uniquely_covered_bases == 0) and (self.num_escrow_covered_bases == 0):
            return None
        
        # Depth of cov should be added since escrow depth is already adjusted by escrow read weights.
        self.depth_of_coverage = self.escrow_depth + self.singular_depth
        
        # Breadth should be the union of Escrow and Singular coverage
        self.breadth_of_coverage = len(self.uniquely_covered_bases.union(self.escrow_covered_bases))
        self.percent_coverage = self.breadth_of_coverage * 100 / self.genome_size

        # Right now, just prioritizing reads with exclusive matches to the genome and their variance. 
        # Might not be the best way forward for cases where the exact genome is missing but multiple similar strains are present.
        self.depth_variance = self.singular_depth_var 
        self._normalize_votes()

        frac_escrow_reads = 0
        if Reads.total_escrow_reads > 0:
            frac_escrow_reads = self.num_escrow_reads/Reads.total_escrow_reads

        frac_singular_reads = 0
        if Reads.total_singular_reads_after_recruitment > 0:
            frac_singular_reads = self.num_singular_reads/Reads.total_singular_reads_after_recruitment
        
        # Dataframe
        return pd.DataFrame(
            index = [self.name],
            data  = {
                'Genome_Size' : self.genome_size,
                'Percent_Coverage' : self.percent_coverage,
                'Total_Bases_Covered' : self.breadth_of_coverage,
                'Coverage_Depth' : self.depth_of_coverage,
                'Depth_Variation' : self.depth_variance,
                'Read_Fraction' : self.read_fraction,
                'Singular_Strain_Weight' : self.adj_primary_wt,
                'Total_Singular_Reads' : self.num_singular_reads,
                'Total_Singular_Votes' : self.cum_primary_votes,
                'Singular_Read_Vote_Ratio' : self.singular_vote_conversion_rate,
                'Singular_Fraction_of_Singular_Reads' : frac_singular_reads,
                'Singular_Coverage' : self.num_uniquely_covered_bases,
                'Singular_Depth' : self.singular_depth,
                'Total_Escrow_Reads' : self.num_escrow_reads,
                'Total_Escrow_Votes' : self.cum_escrow_votes,
                'Escrow_Read_Vote_Ratio' : self.escrow_vote_conversion_rate,
                'Fraction_of_all_Escrow_Reads' : frac_escrow_reads,
                'Escrowed_Cov' : self.num_escrow_covered_bases,
                'Escrowed_Depth' : self.escrow_depth
            }
        )

    def compile_by_abundance(self):
        '''
        return a pandas data frame with 1 row x 4 columns. 
            Strain_Name,Read_Fraction, Percent_Coverage, Coverage_Depth, Depth_Coeff_Variation
        '''
        if self.read_fraction == 0:
            return None

        return pd.DataFrame(
                    index = [self.name],
                    data  = {
                        'Read_Fraction' : self.read_fraction,
                        'Percent_Coverage' : self.percent_coverage,
                        'Coverage_Depth' : self.depth_of_coverage,
                        'Depth_Coeff_Variation' : self.depth_variance
                        }
                    )
    def strain_promiscuity_adjustment(self):
        # self.adj_primary_wt = calculate_sunits_original_adjustment(self)
        self.adj_primary_wt = self.beta_adjustment()
        Strains.net_promiscuity += self.adj_primary_wt
        return self.adj_primary_wt

    # The OG
    def calculate_sunits_original_adjustment(self):
        return self.cum_primary_votes / Reads.total_reads_aligned

    def beta_adjustment(self):
        if self.num_uniquely_covered_bases == 0:
            return 0

        # 2019-11-07 Sunit's interpretation 2 of mike drop: This should be the rate of aggregation of reads for 1 strain compared to the others
        # Borrowed from the beta measure of how well a stock does compared to the rest of the sector.
        strains_reads_per_base = self.num_singular_reads / self.num_uniquely_covered_bases / self.genome_size
        
        if (Reads.total_singular_reads_after_recruitment == self.num_singular_reads):
            # all singular reads assigned were assigned to this strain.
            others_singular_reads_aligned = self.num_singular_reads
            others_uniquely_covered_bases = self.num_uniquely_covered_bases
        elif(Reads.total_singular_reads_after_recruitment > self.num_singular_reads):
            others_singular_reads_aligned = Reads.total_singular_reads_after_recruitment - self.num_singular_reads
            others_uniquely_covered_bases = Strains.total_uniquely_covered_bases - self.num_uniquely_covered_bases
        else:
            sys.exit(f"""
            {self.name} others_singular_reads_aligned can't be negative!
            Total singular reads:{Reads.total_singular_reads_after_recruitment}
            Total singular bases:{Strains.total_uniquely_covered_bases}
            Total bases:{Strains.total_genome_size}

            Strain singular reads:{self.num_singular_reads}
            Strain singular bases:{self.num_uniquely_covered_bases}
            Strain Genome Size:{self.genome_size}
            """)
        
        others_total_genomes_size = Strains.total_genome_size - self.genome_size

        others_reads_per_base = others_singular_reads_aligned / others_uniquely_covered_bases / others_total_genomes_size

        adj_primary_wt = strains_reads_per_base / others_reads_per_base
        # Case 1: self.adj_primary_wt = 1
        #       - rate of read recruitment for this strain is equal to the rate for the rest of the genomes.
        #       - No biases.
        # Case 2: self.adj_primary_wt > 1
        #       - rate of recruitment for this strain is greater than the rate for the rest of the genomes.
        #       - reads preferentially map to this strain over others.
        # Case 3: self.adj_primary_wt < 1
        #       - rate of recruitment for this strain is less than the rate for the rest of the genomes.
        #       - reads preferentially map to other strains over this strain.
        if adj_primary_wt < 0:
            sys.exit(f"""
            {self.name} adj_primary_wt can't be negative!
            Total singular reads:{Reads.total_singular_reads}
            Total singular bases:{Strains.total_uniquely_covered_bases}
            Total bases:{Strains.total_genome_size}
            
            Strain singular reads:{self.num_singular_reads}
            Strain singular bases:{self.num_uniquely_covered_bases}
            Strain Genome Size:{self.genome_size}
            """)
        return adj_primary_wt

    def calculate_read_fraction(self):
        if (self.num_uniquely_covered_bases == 0) and (self.num_escrow_covered_bases == 0):
            self.read_fraction = 0
        else:
            self.read_fraction = (self.cum_escrow_votes + self.cum_primary_votes) * 100 / Reads.total_reads_aligned
        return self.read_fraction

class Reads:
    total_reads_aligned = 0
    reads_w_perfect_alignments = 0
    total_singular_reads = 0
    total_escrow_reads_kept = 0
    total_escrow_reads_discarded = 0
    total_escrow_reads = 0
    total_singular_reads_after_recruitment = 0

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
            common_strains_list = Utils.intersection(read.mapped_strains.keys(), mate.mapped_strains.keys())
            if len(common_strains_list) == 1:
                return common_strains_list[0].name
            else:
                return None

    @staticmethod
    def is_perfect_alignment(aln):
        edit_dist = aln.get_tag('NM')
        query_len = aln.query_length
        ref_start = aln.reference_start
        ref_end = aln.reference_end
        
        # https://www.biostars.org/p/106126/
        return ((edit_dist == 0) and (query_len == aln.get_overlap(ref_start, ref_end)))
    
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

class Utils():
    @staticmethod
    def calculate_coverage(bamfile_name, output_dir):
        pybedtools.set_tempdir(f'{output_dir}/tmp')
        a = pybedtools.BedTool(bamfile_name)
        df = a.genome_coverage(dz = True).to_dataframe(names=['contig','pos', 'depth'])
        pybedtools.cleanup()
        return df

    @staticmethod
    def get_aln_quality(aln):
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

        return (pid, aln_len)

    @staticmethod
    def is_perfect_alignment(pid, aln_len, min_perc_id = 100, min_perc_aln = 100):
        # (pid, aln_len) = get_aln_quality(aln)
        # https://www.biostars.org/p/106126/
        # return ((edit_dist == 0) and (query_len == aln.get_overlap(ref_start, ref_end)))
        if (pid is None):
            return None
        
        if (aln_len is None):
            return None

        return (min_perc_id <= pid <= 100) and (min_perc_aln <= aln_len <= 100)

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

    @staticmethod
    def calculate_coverage(bamfile_name):
        a = pybedtools.BedTool(bamfile_name)
        # b = a.genome_coverage()
        # df = b.to_dataframe(names=['contig', 'depth', 'num_bases', 'contig_size', 'fraction'])
        b = a.genome_coverage(d = True)
        df = b.to_dataframe(names=['contig','pos', 'depth'])
        return df

    @staticmethod
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
    
    @staticmethod
    def intersection(list1, list2):
        return list(set(list1) & set(list2))

    # @staticmethod
    # def bam_is_empty(fn):
    #     if os.path.getsize(fn) > 1000000:
    #         return False

    #     bam = pysam.Samfile(fn, check_sq=False)
    #     try:
    #         bam.next()
    #         return False
    #     except StopIteration:
    #         return True

    @staticmethod
    def sort_and_index(file_name, cores=4, by='Coord'):
        """ Sorts and indexes a bam file by coordinates.
        """
        if by == 'Coord':
            sorted_name = file_name.replace('.bam', '') + '.sortedByCoord.bam'
            # pysam sort multithreading support doesn't work
            # pysam.sort('-@',cores,'-o',sorted_name, file_name)
            pysam.sort('-o',sorted_name, file_name)
            
        elif by == 'Name':
            sorted_name = file_name.replace('.bam', '') + '.sortedByName.bam'
            pysam.sort('-n','-o',sorted_name, file_name)
        else:
            raise Exception("Bam file can only be sorted by 'Coord' or 'Name'.")

        pysam.index(sorted_name)
        os.remove(file_name)
        return sorted_name