#!/usr/bin/env python3

import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam

class Strains:
    total_genome_size = 0
    total_strains = 0
    total_uniquely_covered_bases = 0
    net_promiscuity = 0
    strains_list = list()
    
    def __init__(self, strain_name):
        self.name = strain_name

        self.db_weight = 0
        self.indexed_uniquely_covered_bases = 0
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
        
    def add_paired_singular_vote(self, read_name, mate_name, read_vote = 1):
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
        
    
    def add_escrow_vote(self, read_name, read_vote):
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

    def compile_general_stats(self, total_escrow_reads, total_singular_reads_after_recruitment):
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
        # if Reads.total_escrow_reads > 0:
        #     frac_escrow_reads = self.num_escrow_reads/Reads.total_escrow_reads
        if total_escrow_reads > 0:
            frac_escrow_reads = self.num_escrow_reads/total_escrow_reads

        frac_singular_reads = 0
        # if Reads.total_singular_reads_after_recruitment > 0:
        #     frac_singular_reads = self.num_singular_reads/Reads.total_singular_reads_after_recruitment
        if total_singular_reads_after_recruitment > 0:
            frac_singular_reads = self.num_singular_reads/total_singular_reads_after_recruitment
        
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
    def calculate_sunits_original_adjustment(self, total_reads_aligned):
        # return self.cum_primary_votes / Reads.total_reads_aligned
        return self.cum_primary_votes / total_reads_aligned

    def beta_adjustment(self, total_singular_reads_after_recruitment):
        # where, total_singular_reads_after_recruitment is a class varaiable from Reads class
        if self.num_uniquely_covered_bases == 0:
            return 0

        # 2019-11-07 Sunit's interpretation 2 of mike drop: This should be the rate of aggregation of reads for 1 strain compared to the others
        # Borrowed from the beta measure of how well a stock does compared to the rest of the sector.
        strains_reads_per_base = self.num_singular_reads / self.num_uniquely_covered_bases / self.genome_size
        
        if (total_singular_reads_after_recruitment == self.num_singular_reads):
            # all singular reads assigned were assigned to this strain.
            others_singular_reads_aligned = self.num_singular_reads
            others_uniquely_covered_bases = self.num_uniquely_covered_bases
        elif(total_singular_reads_after_recruitment > self.num_singular_reads):
            others_singular_reads_aligned = total_singular_reads_after_recruitment - self.num_singular_reads
            others_uniquely_covered_bases = Strains.total_uniquely_covered_bases - self.num_uniquely_covered_bases
        else:
            sys.exit(f"""
            {self.name} others_singular_reads_aligned can't be negative!
            Total singular reads:{total_singular_reads_after_recruitment}
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
            Total singular reads:{total_singular_reads_after_recruitment}
            Total singular bases:{Strains.total_uniquely_covered_bases}
            Total bases:{Strains.total_genome_size}
            
            Strain singular reads:{self.num_singular_reads}
            Strain singular bases:{self.num_uniquely_covered_bases}
            Strain Genome Size:{self.genome_size}
            """)
        return adj_primary_wt

    def calculate_read_fraction(self, total_reads_aligned):
        # where total_reads_aligned is a class variable from Reads class
        if (self.num_uniquely_covered_bases == 0) and (self.num_escrow_covered_bases == 0):
            self.read_fraction = 0
        else:
            self.read_fraction = (self.cum_escrow_votes + self.cum_primary_votes) * 100 / total_reads_aligned
        return self.read_fraction

    def add_metadata(self, strain_wt, strain_score, strain_unique_bases):
        self.db_weight = strain_wt
        self.uniqueness_score = strain_score
        self.indexed_uniquely_covered_bases = strain_unique_bases

    @classmethod
    def get_indexed_db_metadata(cls, binmap_file):
        '''
        purpose
            Instantiates Strain objects
        accepts
            bin file created using the ninjaIndex.py script on your database.
        returns two dictionaries
            1) contig name --> strain name
            2) strain name --> strain object
        '''
        bins = defaultdict()
        all_strain_obj = defaultdict()
        # all_strains = defaultdict(list)
        with open(binmap_file, "r") as binmap:
            next(binmap) # skip first/header line.
            for line in binmap:
                line=line.rstrip()
                strain_name, strain_wt_1, strain_wt_2, strain_wt_3, strain_score, strain_unique_bases, contig_name, contig_length = line.split(',')
                # The chosen one
                strain_wt = strain_wt_2

                strain_name = os.path.basename(strain_name)
                if strain_name in all_strain_obj.keys():
                    strain = all_strain_obj[strain_name]
                else:
                    strain = Strains(strain_name)
                    strain.add_metadata(float(strain_wt), float(strain_score), float(strain_unique_bases))
                    all_strain_obj[strain_name] = strain
                    cls.total_strains += 1

                strain.add_contig(contig_name, int(contig_length))

                bins[contig_name] = strain_name
                cls.total_genome_size += int(contig_length)
        
        cls.strains_list = sorted(list(all_strain_obj.keys()), key=str.lower)
        return bins, all_strain_obj

    @classmethod
    def get_db_metadata(cls, binmap_file, fastafile_name):
        '''
        purpose
            Instantiates Strain objects
        accepts binmap file and fasta file
            bin map file contains two columns and no header:
                1) Contig name as it appears in the fasta file (without the '>')
                2) strain name
        returns two dictionaries
            1) contig name --> strain name
            2) strain name --> strain object
        '''
        all_strain_obj = defaultdict()
        bins = defaultdict()
        fasta = pysam.FastaFile(fastafile_name)

        with open(binmap_file, "r") as binmap:
            for line in binmap:
                line=line.rstrip()
                contig_name, strain_name = line.split('\t')
                strain_name = os.path.basename(strain_name)
                bins[contig_name] = strain_name

                if strain_name in all_strain_obj.keys():
                    strain = all_strain_obj[strain_name]
                else:
                    strain = Strains(strain_name)
                    all_strain_obj[strain_name] = strain
                    cls.total_strains += 1
                
                strain.add_contig(contig_name, fasta.get_reference_length(contig_name))
                
                cls.total_genome_size += fasta.get_reference_length(contig_name)

        cls.strains_list = sorted(list(all_strain_obj.keys()), key=str.lower)
        return bins, all_strain_obj
