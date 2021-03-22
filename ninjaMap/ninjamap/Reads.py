#!/usr/bin/env python3

# import pysam
import sys

from collections import defaultdict

class Reads:
    total_reads_aligned = 0
    reads_w_valid_alignments = 0
    total_singular_reads = 0
    total_escrow_reads_kept = 0
    total_escrow_reads_discarded = 0
    total_escrow_reads = 0
    total_singular_reads_after_recruitment = 0

    def __init__(self, name, mate_name, read_length):
        self.name = name
        self.unique_name = name
        self.mates_unique_name = mate_name
        self.read_length = abs(read_length)

        self.cum_vote = 0
        self.has_voted = False
        self.in_singular_bin = False
        self.mate_has_valid_match = False

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
