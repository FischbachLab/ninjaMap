#!/usr/bin/env python3

import argparse
# import logging
import os
import pybedtools
import pysam
import subprocess

# def bam_is_empty(fn):
#     if os.path.getsize(fn) > 1000000:
#         return False

#     bam = pysam.Samfile(fn, check_sq=False)
#     try:
#         bam.next()
#         return False
#     except StopIteration:
#         return True

# def sort_and_index(file_name, cores=4, by='coord'):
#     """ Sorts and indexes a bam file by coordinates.
#     """
#     if by.lower() == 'coord':
#         sorted_name = file_name.replace('.bam', '') + '.sortedByCoord.bam'
#         # pysam sort multithreading support doesn't work
#         # pysam.sort('-@',cores,'-o',sorted_name, file_name)
#         pysam.sort('-o',sorted_name, file_name)
#     elif by.lower() == 'name':
#         sorted_name = file_name.replace('.bam', '') + '.sortedByName.bam'
#         pysam.sort('-n','-o',sorted_name, file_name)
#     else:
#         raise Exception("Bam file can only be sorted by 'coord' or 'name'.")

#     pysam.index(sorted_name)
#     os.remove(file_name)
#     return sorted_name

def sort_and_index(file_name, by='coord', cores=4, memPerCore=2):
    """
    Sorts and indexes a bam file by coordinates.
    """
    sorted_name = None
    if by.lower() == 'coord':
        sorted_name = file_name.replace('.bam', '') + '.sortedByCoord.bam'

        sort_command = f'samtools sort -@ {cores} -m {memPerCore}G -o {sorted_name} {file_name}'
        subprocess.run(sort_command, shell=True, check=True)

        subprocess.run(f'samtools index -@ {cores} {sorted_name}', shell=True, check=True)
    elif by.lower() == 'name':
        sorted_name = file_name.replace('.bam', '') + '.sortedByName.bam'

        sort_command = f'samtools sort -n -@ {cores} -m {memPerCore}G -o {sorted_name} {file_name}'
        subprocess.run(sort_command, shell=True, check=True)
    else:
        raise Exception("Bam file can only be sorted by 'coord' or 'name'.")

    os.remove(file_name)
    return sorted_name

def get_bam_filehandle(bamfile_name, prefix, tag):
    tmp_bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
    new_bamfile_name = f'{prefix}.ninjaMap.{tag}.bam'
    new_bamfile_handle = pysam.AlignmentFile(new_bamfile_name, "wb", template=tmp_bamfile)
    tmp_bamfile.close()
    return (new_bamfile_name, new_bamfile_handle)

def bam2sam(bamfile_name, cores = 4):
    samfile_name = bamfile_name.replace('.bam', '.sam')
    if not os.path.exists(samfile_name):
        subprocess.run(f'samtools view -h -@ {cores} -o {samfile_name} {bamfile_name}', shell=True, check=True)
    return samfile_name

def sam2bam(samfile_name, sort_by='coord', cores = 4):
    bamfile_name = samfile_name.replace('.sam', '.bam')
    subprocess.run(f'samtools view -h -b -@ {cores} -o {bamfile_name} {samfile_name}', shell=True, check=True)
    sorted_bam = sort_and_index(bamfile_name, by=sort_by)
    # TODO: return a header file and a dataframe of SAM file
    # os.remove(samfile_name)
    return sorted_bam

def calculate_coverage(bamfile_name, output_dir):
    os.makedirs(f'{output_dir}/tmp', exist_ok=True)
    pybedtools.set_tempdir(f'{output_dir}/tmp')
    bed = pybedtools.BedTool(bamfile_name)
    df = bed.genome_coverage(dz = True).to_dataframe(names=['contig','pos', 'depth'])
    pybedtools.cleanup()
    return df

def calculate_detailed_coverage(bamfile_name, outputdir):
    # sambamba depth \
    #     base \
    #     -o sambamba_depth.singular.out \
    #     -t 15 \
    #     -c 0 \
    #     -q 10 \
    #     -m \
    #     sortedByCoord.bam
    # sambamba output: REF     POS     COV     A       C       G       T       DEL     REFSKIP SAMPLE
    pass