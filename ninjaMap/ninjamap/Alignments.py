#!/usr/bin/env python3

# import logging
import os
import sys
import pysam
# import subprocess

class Alignments():
    total_alignments = total_valid_alignments = total_bad_aln = 0
    total_reads = set()
    reads_w_valid_alignments = set()
    def __init__(self, aln, paired):
        self.pid = self.p_aln_len = self.edit_dist = self.aln_len = float(0)
        self.match = self.insertions = self.deletions = float(0)
        self.precision = self.recall = self.f_score = float(0)
        self.read_length = self.perc_indel = float(0)
        self.read_base = self.read_name = self.mate_name = self.orientation = self.reference_name = ''
        self.is_valid_aln = self.is_aln = False
        Alignments.total_alignments += 1
        self.paired = paired
        self.extract_read_names(aln)
        Alignments.total_reads.add(self.read_name)
        if aln.has_tag('NM'):
            self.is_aln = True
            self.reference_name = aln.reference_name
            (self.read_length, self.pid, self.p_aln_len, self.edit_dist,
            self.aln_len, self.match, self.insertions, self.deletions,
            self.precision, self.recall, self.f_score) = self.extract_aln_info(aln)
        else:
            self.is_aln = False

    def is_valid_alignment(self, min_perc_id, min_perc_aln, max_perc_indel, min_f_score):
        """
        method to check if alignment is valid using the given thresholds.

        Parameters
        ----------
        min_perc_id (float),
        min_perc_aln (float),
        max_perc_indel (float),
        min_f_score (float)

        Returns
        -------
        boolean
        """
        self.is_valid_aln = False
        if not self.is_aln:
            Alignments.total_bad_aln += 1
            return self.is_valid_aln

        if self.paired:
            if (self.pid < min_perc_id) or (self.p_aln_len < min_perc_aln):
                Alignments.total_bad_aln += 1
                return self.is_valid_aln
        else:
            self.perc_indel = (self.insertions+self.deletions)*100/self.aln_len
            if self.perc_indel > max_perc_indel:
                Alignments.total_bad_aln += 1
                return self.is_valid_aln
            
            if self.f_score < min_f_score:
                Alignments.total_bad_aln += 1
                return self.is_valid_aln

        self.is_valid_aln = True
        Alignments.total_valid_alignments += 1
        Alignments.reads_w_valid_alignments.add(self.read_name)
        return self.is_valid_aln

    def extract_read_names(self,aln):
        if self.paired:
            self.read_base = Alignments._parse_read_name(aln)
            self.read_name = self._get_unique_read_name(aln)
            self.mate_name = self._get_unique_mate_name(aln)
        else:
            self.read_base = self.read_name = Alignments._parse_read_name(aln)
            self.orientation = 'fwd'
            self.mate_name = ''
        
        # return (read_name, mate_name)

    def extract_aln_info(self, aln):
        poor_aln = [None] * 11
        read_length = aln.infer_read_length()
        # if not aln.has_tag('NM'):
        #     return poor_aln
        
        # http://samtools.github.io/hts-specs/SAMv1.pdf (page 9; CIGAR String)
        # Sum of  aln_match + insertions + BAM_CSOFT_CLIP + match + mismatch = Aligned Segment Length
        (aln_match, insertions, deletions, 
        BAM_CREF_SKIP, BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CPAD, 
        match, mismatch, BAM_CBACK, edit_dist) = aln.get_cigar_stats()[0]

        # query_len = get_read_length(aln)
        qlen = aln_match + insertions + BAM_CSOFT_CLIP + BAM_CHARD_CLIP + match + mismatch
        if qlen != read_length:
            print(f"Calculated query length({qlen}) is not equal to the given query length:{read_length}")
            sys.exit(1)

        aln_len = aln.query_alignment_length
        # if match != aln_len:
        #     print(f"Calculated alignment length({match}) is not equal to the bam alignment length:{aln_len}")
        #     sys.exit(1)

        # v1
        # pid = (query_len - edit_dist)*100/query_len
        # p_aln_len = aln_len*100/query_len

        # from IGGSearch (https://github.com/snayfach/IGGsearch/blob/master/iggsearch/search.py#L250)
        # and NinjaMap
        pid = 100*(aln_len - edit_dist)/float(aln_len)
        p_aln_len = 100*aln_len/float(read_length)

        # v2
        # pid : num of residues that are identical between query and subject
        # aln_len : length of the aligned region (incl. match, mismatch, ins, del, etc)
        # Hypothesis: 
        #   pid ~= precision: num of true matches/num of aligned bases
        #   qcov ~= recall: num of true matches/max possible matches (query length)
        precision = recall = f_score = 0
        if match > 0:
            precision = 100*match/float(aln_len)
            recall = 100*match/float(read_length)
            f_score = Alignments.f_beta_score(precision,recall,2) # f2 score to weigh recall higher
        # print(f'{aln_match}, {insertions}, {deletions}, {BAM_CREF_SKIP}, {BAM_CSOFT_CLIP}, {BAM_CHARD_CLIP}, {BAM_CPAD}, {match}, {mismatch}, {BAM_CBACK}, {edit_dist}')
        # print(f'{read_length}, {pid}, {p_aln_len}, {edit_dist}, {aln_len}, {match}, {insertions}, {deletions}, {precision}, {recall}, {f_score}')
        # if Alignments.total_alignments > 10:
        #     sys.exit(1)
        return(read_length, pid, p_aln_len, edit_dist, aln_len, match, insertions, deletions, precision, recall, f_score)

    def _get_unique_read_name(self, aln):
        if aln.is_read1:
            self.orientation =  'fwd'
        else:
            self.orientation =  'rev'
            
        return self.read_base +'__'+ self.orientation

    def _get_unique_mate_name(self, aln):
        orientation = ''
        if aln.is_read1:
            orientation =  'rev'
        else:
            orientation =  'fwd'
            
        return self.read_base +'__'+ orientation

    @staticmethod
    def _parse_read_name(aln):
        '''
        Accept: AlignmentFile object from PySam
        if read name has a '/', this is the old illumina format. 
        strip the content after the '/', return remaining
        else, return it as is.
        '''
        if ' ' in str(aln.query_name):
            space_delim_name, *value = aln.query_name.split(' ')
        else:
            space_delim_name = str(aln.query_name)
        
        if '/' in space_delim_name:
            inferred_query_name, *value = space_delim_name.split('/')
        else:
            inferred_query_name = space_delim_name

        return inferred_query_name

    @staticmethod
    def f_beta_score(precision, recall, beta = 1):
        beta_sq = (beta)**2
        f_beta = (1 + beta_sq) * (precision * recall)/float((beta_sq * precision) + recall)
        return f_beta

    # @classmethod
    # def fetch_next_alignment(cls, bamfile_name, filtered_bamfile_handle,paired, max_perc_indel, min_f_score):
    #     bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
    #     for alignment in bamfile.fetch(until_eof=True):
    #         aln = Alignments(alignment, paired)
    #         cls.total_reads.add(aln.read_name)
    #         # pId, pAln, edit_dist, aln_len, insertions, deletions = get_aln_quality(aln)
    #         # if aln.is_aln:
    #         #     aln.extract_read_info(alignment)
    #         if aln.is_valid_alignment(max_perc_indel, min_f_score):
    #             # write to a new BAM file
    #             filtered_bamfile_handle.write(alignment)
    #             yield aln
    #     bamfile.close()
    #####################################
    #####################################
    # def get_aln_quality(aln):
    #     poor_aln = (None, None, None, None)
    #     if not aln.has_tag('NM'):
    #         return poor_aln
        
    #     edit_dist = aln.get_tag('NM')
    #     query_len = aln.query_length
    #     # ref_start = aln.reference_start
    #     # ref_end = aln.reference_end
        
    #     if not query_len:
    #         return poor_aln

    #     # aln_len = aln.get_overlap(ref_start, ref_end)
    #     aln_len = aln.query_alignment_length
    #     pid = (query_len - edit_dist)*100/query_len
    #     p_aln_len = aln_len*100/query_len

    #     return (pid, p_aln_len, edit_dist, aln_len)
