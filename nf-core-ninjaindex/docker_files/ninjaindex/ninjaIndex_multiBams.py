#!/usr/bin/env python3

import argparse
import gzip
import logging
import os
import pysam
import pandas as pd
import re
import sys

from Bio import SeqIO
from collections import defaultdict, Counter
from time import perf_counter as timer
import multiprocessing as mp
from multiprocessing import Process, Manager
from multiprocessing import Pool
import threading
#import psutil

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

###############################################################################
# Setup Input and Script Usage
###############################################################################

usage = """
    USAGE:
    python ninjaIndex.py \
-selfbamdir folder with each individually aligned bamfile\
-fastadir folder with each genome in a separate file in fasta format. \
-prefix [Output] output prefix for a comma-delimited binmap file and a combined fasta file
    """

p = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter,
    add_help=True,
    usage=argparse.SUPPRESS,
    description="""Description:
This script will calculate the strain weights in a defined microbial community.
Usage: ninjaIndex.py -selfbamdir bam_dir -fastadir fasta_dir -prefix index_name
""",
    epilog="""
""")
#p.add_argument('-bam', dest='bamfile', action='store', type=str, #required = True,
#                help='name sorted bam file.')
p.add_argument('-selfbamdir', dest='selfbamdir', action='store', type=str, required = True,
                help='dir with self aligned genomes in separate bam files.')
p.add_argument('-fastadir', dest='fastadir', action='store', type=str, required = True,
                help='dir with genomes in separate fasta files.')
p.add_argument('-"fasta_ext"', dest='ext', action='store', type=str, default = "fna",
                help='file file extension; all files must have the same extension.')
p.add_argument('-prefix', dest='prefix', action='store', type=str,
                help='[Output] output prefix. Watch out for a comma-delimited binmap file and a concatenated fasta file for ninjaMap')

# Optional
p.add_argument('-outdir', dest='outdir', action='store', type=str,
                help='output directory')
p.add_argument('-threads', dest='threads', action='store', type=str,
                help='number of threads available for this job and subprocesses')
p.add_argument('-debug', dest='debug', action='store_true', default=False,
                help='save intermediate false positives bam file')

args = vars(p.parse_args())

#bamfile_name = args['bamfile']
fastafile_dir = args['fastadir']
fasta_ext = args['ext']

selfbamfile_dir = args['selfbamdir']

if not args['threads']:
    cpus = mp.cpu_count()
else:
    cpus = args['threads']

if not args['outdir']:
    output_dir = "db"
else:
    output_dir = args['outdir']

os.makedirs(output_dir, exist_ok=True)

if not args['prefix']:
    prefix = os.path.join(output_dir, args['prefix'])
#    default_prefix = os.path.basename(bamfile_name).split('.')[0]
#    prefix = os.path.join(output_dir, default_prefix)
else:
    prefix = args['prefix']

binmap_file = prefix +'.ninjaIndex.binmap.csv'
fasta_file = prefix + ".ninjaIndex.fasta"

logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s\t[%(levelname)s]:\t%(message)s')

logging.info('Started')
###############################################################################
# Classes
###############################################################################
class Strains():
    total_genome_size = 0
    total_strains = 0
    total_uniquely_covered_bases = 0

    def __init__(self, strain_name):
        self.name = os.path.basename(strain_name)
        self.num_singular_reads = 0
        self.num_contigs = 0
        self.genome_size = 0
        self.cum_primary_votes = 0
        self.uniqueness_score = 0

        self.total_covered_bases = 0
        self.uniquely_covered_bases = 0
        self.adj_primary_wt = 0
        self.adj_primary_wt_2 = 0
        self.adj_primary_wt_3 = 0

        self.aln_norm_abundance = 0
        self.genome_norm_abundance = 0
        self.adjusted_votes = 0
        self.singular_vote_conversion_rate = 0

        # reads aligned to its genome itself
        self.self_aligned_reads_number = 0

        self.contigs = defaultdict(int)
        self.singular_bin = defaultdict(int)

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
        # for R1
        self.singular_bin[read_name] += read_vote
        self.cum_primary_votes += read_vote
        self.num_singular_reads += 1

        # for R2
        self.singular_bin[mate_name] += read_vote
        self.cum_primary_votes += read_vote
        self.num_singular_reads += 1

    def _calc_base_depth(self, pileupcolumn, bin_dict):
        wt_base_depth = 0
        for pileupread in pileupcolumn.pileups:
            # For all reads aligned to this base/column on the contig
            # print ('\tbase in read %s = %s' %
            #       (pileupread.alignment.query_name,
            #        pileupread.alignment.query_sequence[pileupread.query_position]))
            read_unique_name = Parser.get_unique_read_name(pileupread.alignment)
            if read_unique_name in bin_dict.keys():
                # Only calculate depth from reads with perfect alignment(singular or escrow).
                wt_base_depth += bin_dict[read_unique_name]

        return wt_base_depth

    def calculate_singular_coverage (self, bamfile_name): #, cov_bamfile):
        #print('Openning bam file name : ', bamfile_name)
        if self.num_singular_reads == 0:
            return

        cov_bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')

        for contig_name in self.contigs.keys():

             #pattern = re.compile(self.name)
             #m = pattern.match(contig_name)
             #if m:
            #    print('Match found: ', m.group())

            for pileupcolumn in cov_bamfile.pileup(contig = contig_name, stepper = 'all', min_base_quality = 20):
                base_cov_contribution = 0
                if pileupcolumn.nsegments > 0:
                    base_cov_contribution = self._calc_base_depth(pileupcolumn, self.singular_bin)

                    if base_cov_contribution > 0:
                        self.uniquely_covered_bases += 1
                        Strains.total_uniquely_covered_bases += 1
             #else:
            #     print('No match on ', contig_name)
        cov_bamfile.close()
        #print('*** Done : ', bamfile_name)
        return (self.uniquely_covered_bases, self.name, Strains.total_uniquely_covered_bases)

    #def run(self, bamfile):
    #        t = threading.Thread(target=self.calculate_singular_coverage, args=(bamfile,))
            #t.start()
            #t.join()
            #print()

    def calculate_andres_adjustment(self):
        # Andres's interpretation of Mike drop adjustment.
        # Actually in his interpretation, 'unique # bases' needed to be calculated for all the strains
        # a read is assigned to => Size of unique region in a genome compared to other genomes that this read maps.
        # This is a VERY costly operation to repeat millions of times, since it cannot be pre-calculated.
        # I have adjusted it to use the genome's absolute unique region compared to all genomes in the database.
        #print(self.genome_size, self.uniquely_covered_bases)  self.genome_size -
        #print (self.name, self.genome_size, self.uniquely_covered_bases)
        self.adj_primary_wt = self.cum_primary_votes / (self.uniquely_covered_bases) #if self.uniquely_covered_bases != 0 else 0
        return self.adj_primary_wt

    def calculate_sunits_original_adjustment(self):
        self.adj_primary_wt_2 = self.cum_primary_votes / Reads.total_reads_aligned #if Reads.total_reads_aligned != 0 else 0
        return self.adj_primary_wt_2

    def calculate_mike_drop_penalty(self):
        self.uniqueness_score = self.uniquely_covered_bases / Strains.total_uniquely_covered_bases  #if Strains.total_uniquely_covered_bases != 0 else 0
        return self.uniqueness_score

    def calculate_adj_primary_wt_uniqueness_score(self):
        # Sunit's interpretation of Mike drop adjustment
        self.adj_primary_wt_3 = (1 - self.uniqueness_score) * self.cum_primary_votes / Reads.total_reads_aligned # if Reads.total_reads_aligned != 0 else 0
        return self.adj_primary_wt_3
class Reads:
    total_reads_aligned = 0
    reads_w_perfect_alignments = 0
    total_singular_reads = 0
    total_singular_reads_in_pairs = 0

    def __init__(self, name, mate_name, read_length, template_length):
        self.name = name
        self.unique_name = name
        self.mates_unique_name = mate_name
        self.read_length = read_length
        self.template_length = template_length

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
        vote_value = round(vote_value, 5)
        self.cum_vote += vote_value
        self.has_voted = True

    def is_fraud(self):
        '''
        for each object of this class, return True if cumulative votes > 1
        '''
        return (round(self.cum_vote, 5) > 1)

class Parser:
    @staticmethod
    def choose_primary_candidate(read, mate):
        if (read.template_length == 0) or (mate.template_length == 0):
            return None

        if read.mate_has_perfect_match or mate.mate_has_perfect_match :
            common_strains_list = intersection(read.mapped_strains.keys(), mate.mapped_strains.keys())
            if len(common_strains_list) == 1:
                return common_strains_list[0].name
            else:
                return None

    @staticmethod
    def is_perfect_alignment(aln):
        edit_dist = dict(aln.tags)['NM']
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

        return Parser.parse_read_name(aln) +'__'+ orientation

    @staticmethod
    def get_unique_mate_name(aln):
        orientation = ''
        if aln.is_read1:
            orientation =  'rev'
        else:
            orientation =  'fwd'

        return Parser.parse_read_name(aln) +'__'+ orientation

    @staticmethod
    def extract_read_info(aln):
        read_name = Parser.get_unique_read_name(aln)
        mate_name = Parser.get_unique_mate_name(aln)
        read_length = aln.reference_length
        template_length = aln.template_length

        return (read_name, mate_name, read_length, template_length)

    @staticmethod
    def read_bam_file(bamfile_name, bins0): #, perfect_alignment, read_objects):
        total_reads = 0
        read_objects = defaultdict()
        perfect_alignment = defaultdict(lambda: defaultdict(int))
        discarded_reads = set()

        # Read the BAM file
        bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
        for aln in bamfile.fetch(until_eof=True):
            read_name, mates_name, read_length, template_length = Parser.extract_read_info(aln)
            if Parser.is_perfect_alignment(aln):
                # read_info[read_name] = (mates_name, read_length, template_length)
                read = Reads(read_name, mates_name, read_length, template_length)
                read_objects[read_name] = read
                strain_name = bins[aln.reference_name] # bins[contig_name] --> strain name
                perfect_alignment[read_name][strain_name] += 1
            else:
                discarded_reads.add(read_name)

        total_reads = len(discarded_reads) + len(read_objects.keys())
        bamfile.close()
        return(total_reads, perfect_alignment, read_objects)


    @staticmethod
    def read_selfbam_file(bamfile_name, bins):
        total_reads = 0
        read_objects = defaultdict()
        perfect_alignment = defaultdict(lambda: defaultdict(int))
        discarded_reads = set()

        # Read the BAM file
        bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
        for aln in bamfile.fetch(until_eof=True):
            read_name, mates_name, read_length, template_length = Parser.extract_read_info(aln)
            if Parser.is_perfect_alignment(aln):
                # read_info[read_name] = (mates_name, read_length, template_length)
                read = Reads(read_name, mates_name, read_length, template_length)
                read_objects[read_name] = read
                strain_name = bins[aln.reference_name] # bins[contig_name] --> strain name
                perfect_alignment[read_name][strain_name] += 1
            else:
                discarded_reads.add(read_name)

        total_reads = len(discarded_reads) + len(read_objects.keys())
        bamfile.close()
        return(total_reads, strain_name, perfect_alignment, read_objects)

    @staticmethod
    def get_db_metadata(fasta_list, fasta_filename):
        '''
        read the fasta files
        return bins, all_strain_obj, strains_list, concatenated fasta file
        '''
        bins = defaultdict()

        #manager=Manager()
        #all_strain_obj=manager.dict()
        #global all_strain_obj
        all_strain_obj = defaultdict()

        with open(fasta_filename, "w") as fasta_output:
            for input_file in fasta_list:
                # fasta_sequences = SeqIO.parse(open(input_file),'fasta')
                tmp_strain_name = os.path.basename(input_file)
                strain_name = os.path.splitext(tmp_strain_name)[0]
                strain = Strains(strain_name)
                all_strain_obj[strain_name] = strain
                # strains_list.append(strain_name)
                with open(input_file, "r") as fasta_input:
                    for fasta in SeqIO.parse(fasta_input,'fasta'):
                        contig_name, sequence = fasta.id, str(fasta.seq)
                        contig_len = len(sequence)
                        strain.add_contig(contig_name, contig_len)
                        #print ("contig_name, contig_len: ", contig_name, contig_len)
                        bins[contig_name] = strain_name
                        Strains.total_genome_size += contig_len

                        SeqIO.write(fasta, fasta_output , "fasta")

        return (bins, all_strain_obj, fasta_filename)

    @staticmethod
    def create_bin_map(all_strain_obj, binmap_file):
        binmap = open(binmap_file, 'w')
        # Header
        # strain_name, strain_wt_1, strain_wt_2, strain_wt_3, strain_score, strain_unique_bases, contig_name, contig_length
        binmap.write('Strain_Name,Strain_Weight_Andres,Strain_Weight_Original, Strain_Weight_UScore,Strain_Uniqueness_Score,Strain_Absolute_Unique_Bases,Contig_Name,Contig_Length\n')
        for name, strain in all_strain_obj.items():
            line = ''
            first_half = strain.name +','+ \
                str(strain.adj_primary_wt) +','+ \
                str(strain.adj_primary_wt_2) +','+ \
                str(strain.adj_primary_wt_3) +','+ \
                str(strain.uniqueness_score) +','+ \
                str(strain.uniquely_covered_bases)
                #strain.genome_size-
            for contig_name, contig_len in strain.contigs.items():
                line += str(first_half)  +','+ \
                        contig_name +','+ \
                        str(contig_len) +'\n'
            binmap.write(line)
        binmap.close()
###############################################################################
# Parse the genome fasta files.
###############################################################################
logging.info('Processing the Fasta files from: %s ...', fastafile_dir)

fastafile_names = list()
strains_list = list()
for filename in os.listdir(fastafile_dir):
    full_filename = os.path.join(fastafile_dir, filename)
    if os.path.isfile(full_filename) and filename.endswith(fasta_ext):
        fastafile_names.append(full_filename)

#manager=Manager()
#all_strain_obj=manager.dict()

bins, all_strain_obj, cat_fasta_filename = Parser.get_db_metadata(fastafile_names, fasta_file)
strains_list = all_strain_obj.keys()
Strains.total_strains = len(all_strain_obj.keys())
logging.info('\t%d contigs assigned to %d strains', len(bins.keys()), Strains.total_strains)
###############################################################################
# Parse the individually aligned refernece genome BAM file
###############################################################################
selfbamfile_names = list()
bam_total_reads = 0
#read_objects1 = defaultdict()
all_perfect_alignment =  defaultdict(lambda: defaultdict(int))
all_read_objects = defaultdict()
for filename in os.listdir(selfbamfile_dir):
    full_filename = os.path.join(selfbamfile_dir, filename)
    if os.path.isfile(full_filename) and filename.endswith(".bam"):
        selfbamfile_names.append(full_filename)
        logging.info('Processing the individually aligned BAM file: %s ...', filename)
        read_objects1 = defaultdict()
        mytotal_reads_aligned, strain_name, perfect_alignment1, read_objects1 = Parser.read_selfbam_file(full_filename, bins )
        all_perfect_alignment.update(perfect_alignment1)
        all_read_objects.update(read_objects1)
        #print ("Total reads aligned: ", mytotal_reads_aligned)
        Strains(strain_name).self_aligned_reads_number = mytotal_reads_aligned
        Reads.total_reads_aligned += mytotal_reads_aligned

Reads.reads_w_perfect_alignments = len(all_perfect_alignment.keys())
logging.info('\tUsed %d reads with perfect alignments, out of %d (%7.3f%%).',
    Reads.reads_w_perfect_alignments,
    Reads.total_reads_aligned,
    Reads.reads_w_perfect_alignments*100/Reads.total_reads_aligned if Reads.total_reads_aligned !=0 else 0
    )
#print (selfbamfile_names)
###############################################################################
# Parse the merged BAM file
###############################################################################
if False: '''
...statements...

logging.info('Processing the merged BAM file: %s ...', bamfile_name)
read_objects = defaultdict()
bam_total_reads = 0
bam_total_reads, perfect_alignment, read_objects = Parser.read_bam_file(bamfile_name, bins) #, perfect_alignment1, read_objects1)

Reads.total_reads_aligned += bam_total_reads

all_perfect_alignment.update(perfect_alignment)
all_read_objects.update(read_objects )

#print (" # of read objects", len(all_read_objects.keys()) )
Reads.reads_w_perfect_alignments = len(all_perfect_alignment.keys())
logging.info('\tUsed %d reads with perfect alignments, out of %d (%7.3f%%).',
    Reads.reads_w_perfect_alignments,
    Reads.total_reads_aligned,
    Reads.reads_w_perfect_alignments*100/Reads.total_reads_aligned if Reads.total_reads_aligned !=0 else 0
    )
'''
###############################################################################
# Separate the Primary from the Escrow alignments
# Calculate the Strain abundance distribution based on the primary alignments.
###############################################################################

logging.info('Separating the Primary from the Escrow alignments ...')
for read_name in all_perfect_alignment.keys():
    # mate_name, read_length, template_length = read_info[read_name]
    # read = Reads(read_name, mate_name, read_length, template_length)
    # read_objects[read_name] = read

    read = all_read_objects[read_name]

    if read.mates_unique_name in all_perfect_alignment.keys():
        # This means the mate had a perfect match too.
        read.mate_has_perfect_match = True

    read.num_strains = len(all_perfect_alignment[read_name].keys())

    for strain_name in all_perfect_alignment[read_name].keys():
        read.add_exact_match(all_strain_obj[strain_name])

    if read.num_strains == 1:
        # Singular
        # True hit. This read gets 1 whole vote to assign to a Strain.
        read.in_singular_bin = True
        Reads.total_singular_reads += 1

del all_perfect_alignment
del perfect_alignment1
#del perfect_alignment

logging.info('\tUsed %d reads for primary distribution, out of %d (%7.3f%%) reads with perfect alignments or %7.3f%% of total.',
    Reads.total_singular_reads,
    Reads.reads_w_perfect_alignments,
    Reads.total_singular_reads*100/Reads.reads_w_perfect_alignments if Reads.reads_w_perfect_alignments !=0 else 0,
    Reads.total_singular_reads*100/Reads.total_reads_aligned if Reads.total_reads_aligned !=0 else 0
    )

if len(all_read_objects.keys()) != Reads.reads_w_perfect_alignments:
    logging.critical('Read %d reads with perfect alignments, but created %d read objects. There is something fishy going on here...',
    Reads.reads_w_perfect_alignments,
    len(all_read_objects.keys())
    )
    sys.exit()
###############################################################################
# Processing the reads that mapped exclusively to a single strain
###############################################################################
for name, read in all_read_objects.items():
    if read.has_voted:
        continue

    # if read.mates_unique_name in read_objects.keys():
    if read.mate_has_perfect_match:
        # Means both reads in a mate are prefect alignments
        mate = all_read_objects[read.mates_unique_name]
        strain_name = Parser.choose_primary_candidate(read, mate)

        if strain_name is not None:
            strain = all_strain_obj[strain_name]
            strain.add_paired_singular_vote(read.unique_name, mate.unique_name, 1)
            read.add_vote(1)
            mate.add_vote(1)

            Reads.total_singular_reads_in_pairs += 2


logging.info('\t%d reads will be used for singular alignment strain abundance out of %d (%7.3f%%) reads with perfect alignments or %7.3f%% of total aligned.',
    Reads.total_singular_reads_in_pairs,
    Reads.reads_w_perfect_alignments,
    Reads.total_singular_reads_in_pairs*100/Reads.reads_w_perfect_alignments if Reads.reads_w_perfect_alignments !=0 else 0,
    Reads.total_singular_reads_in_pairs*100/Reads.total_reads_aligned if Reads.total_reads_aligned !=0 else 0)

del all_read_objects
#del read_objects
###############################################################################
# Calculate unique number of bases covered for each genome
###############################################################################
#print("Number of cpu : ", mp.cpu_count())

logging.info('Computing coverage for each strain in the database based on singular alignments ...')
logging.info("Number of cpus is used : %s", cpus)
pool = mp.Pool(processes=int(cpus))
result = []

#Parse by strain name
i = 0
for name, strain in all_strain_obj.items():
    i += 1
    logging.info("\t[%d/%d] Searching for exclusive support for :\t%s",i, Strains.total_strains, name)
    #singular_reads_set = set(strain.singular_bin.keys())
    #strain.calculate_singular_coverage(bamfile_name)
    #logging.info("\tDone search for :\t%s", bamfile_name)
    threads = []
    for selfbamfile in selfbamfile_names:
        #strain.calculate_singular_coverage(selfbamfile)
        #p = Process(target=strain.calculate_singular_coverage, args=(selfbamfile,))
        #th = threading.Thread(target=strain.calculate_singular_coverage, args=(selfbamfile, ))
        #threads.append(th)
        #th.start()
    #for th in threads:
    #    th.join()
        result.append(pool.apply_async(strain.calculate_singular_coverage, args=(selfbamfile,)))

pool.close()
pool.join()
for res in result:
    #print (res.get())
    (my_uniquely_covered_bases, mystrain_name, my_total_uniquely_covered_bases) =res.get()
    all_strain_obj[mystrain_name].uniquely_covered_bases += my_uniquely_covered_bases
    Strains.total_uniquely_covered_bases += my_total_uniquely_covered_bases

    #pool.map(strain.calculate_singular_coverage, selfbamfile_names)
    #pool.close()
    #pool.join()

    #for selfbamfile in selfbamfile_names:
    #    strain.calculate_singular_coverage(selfbamfile)
        #logging.info("\tDone search for :\t%s", selfbamfile)
    #logging.info("\tDone search for :\t%s", name)
logging.info("\tDone search for all strains")

if False: '''
# Parse by BAM files

for selfbamfile in selfbamfile_names:
    #cov_bamfile = pysam.AlignmentFile(selfbamfile, mode = 'rb')
    #procs = []
    threads = []
    for strain_name in strains_list:
        strain = all_strain_obj[strain_name]
        #strain.calculate_singular_coverage(selfbamfile, cov_bamfile )
        #strain.calculate_singular_coverage(selfbamfile)
        # calculage each strain in parallel
        #print(strain_name)
        #pool.apply(strain.calculate_singular_coverage, args=(selfbamfile,))
        #pool.map(strain.calculate_singular_coverage, selfbamfile)
        #th = Process(target=strain.calculate_singular_coverage, args=(selfbamfile,))
        #th = threading.Thread(target=strain.calculate_singular_coverage, args=(selfbamfile, ))
        #threads.append(th)
        #th.start()
        result.append(pool.apply_async(strain.calculate_singular_coverage, args=(selfbamfile,)))
pool.close()
pool.join()
for res in result:
    #print (res.get())
    (my_uniquely_covered_bases, mystrain_name, my_total_uniquely_covered_bases) =res.get()
    all_strain_obj[mystrain_name].uniquely_covered_bases += my_uniquely_covered_bases
    Strains.total_uniquely_covered_bases += my_total_uniquely_covered_bases

        #procs.append(p)
        #p.start()
#for th in threads:
#        th.join()

logging.info("\tDone search all stains") # in :\t%s", selfbamfile)
      # complete the processes
    #print("The length of list is: ", len(procs))
    #for p in procs:
    #    p.join()
    #cov_bamfile.close()


    #logging.info("\tDone search all stains in :\t%s", selfbamfile)
'''
if False: '''
    for selfbamfile in selfbamfile_names:
        pattern = re.compile(name)
        m = pattern.match(selfbamfile)
        if m:
            logging.info("\tStart search for :\t%s", selfbamfile)
            strain.calculate_singular_coverage(selfbamfile)
            logging.info("\tDone search for :\t%s", selfbamfile)
        else:
            pass
'''
# del read_info
###############################################################################
# Calculating Strain Weights
###############################################################################
logging.info('Calculating strain weights based on singular alignments...')
for name, strain in all_strain_obj.items():
    strain.calculate_mike_drop_penalty()
    strain.calculate_andres_adjustment()
    strain.calculate_sunits_original_adjustment()
    strain.calculate_adj_primary_wt_uniqueness_score()

logging.info('Creating the bin map file...')
Parser.create_bin_map(all_strain_obj, binmap_file)

###############################################################################
# The End
###############################################################################
end = timer()
logging.info('Completed in %s (d:h:m:s)', human_time(end - start))
