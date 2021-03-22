#!/usr/bin/env python

import argparse
import gzip
import os
import re
import sys

import pandas as pd
from Bio import SeqIO

## Functions
def usage():
    p = argparse.ArgumentParser(   
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True,
        usage=argparse.SUPPRESS,
        description="""Description:
    Given a list of genome names and accessible URLs (HTTP/FTP/S3), 
    - download the genome, 
    - standardize the genome name and sequence headers, 
    - generate basic stats; and 
    - concatenate and create a db
    Usage: python create_ninjamap_db.py -list genomeName_url.list -db database_name
    """,
        epilog="""Examples:
    python create_ninjamap_db.py  -list scv1_1.s3paths.list -db scv1_1
    """)
    # Required
    p.add_argument('-list', dest='db_ref_file', action='store', type=str, required = True,
                    help='2 columns ("Strain_Name" and "FTP_link") tab-sep or comma-sep file, one genome per line')
    p.add_argument('-db', dest='db_name', action='store', type=str, required = True,
                    help='output database file prefix')
    # Optional
    p.add_argument('-no_cleanup', dest='no_cleanup', action='store_true', default=False,
                help='Do not remove intermediate metadata files (default=True)')
    
    return vars(p.parse_args())

def process_raw_seedfile(seedfile):
    try:
        df = pd.read_csv(seedfile).dropna(axis=0, subset=['Strain_Name', 'FTP_link'])
        df['Strain_Name']
    except:
        df = pd.read_table(seedfile).dropna(axis=0, subset=['Strain_Name', 'FTP_link'])
        df['Strain_Name']
    # replace all punctuations, special characters and spaces with '-'
    df['Strain_Name'] = df['Strain_Name'].apply(lambda x: re.sub('[^A-Za-z0-9]+', '-', x))
    # remove '-' if it appears at the end of the string
    df['Strain_Name'] = df['Strain_Name'].apply(lambda x: re.sub('\-+$', '', x))

    # get reference IDs from the FTP path.
    # This makes the script very specific to the current NCBI FTP paths.
    df['refIDs'] = df['FTP_link'].apply(lambda x: os.path.basename(os.path.dirname(x)))

    return df

def download_genome(genome_name, reference_id, ftp_link, fasta_dir, metadata_dir):
    # max_retries = 5
    
    if ftp_link.startswith('s3://'):
        localname = f'{fasta_dir}/{genome_name}.fna.gz'
        import boto3
        s3 = boto3.client('s3')
        bucket = ftp_link.split('/')[2]
        key = '/'.join(ftp_link.split('/')[3:])
        local_filename = localname
        print(f'B:{bucket}\tK:{key}')
        with open(localname, 'wb') as f:
            try:
                s3.download_fileobj(bucket, key, f)
            except:
                # error_msg = "\n".join(sys.exc_info())
                # print(f'[ERROR] Unable to retrieve reference from URL {ftp_link}\nThe following error occurred: {error_msg}')
                print(f'[ERROR] Unable to retrieve reference from URL {ftp_link}')
                sys.exit(1)
    else:
        localname = f'{fasta_dir}/{genome_name}__{reference_id}.fna.gz'
        try:
            import urllib
            local_filename, _ = urllib.request.urlretrieve(ftp_link, localname)
        except: # catch *all* exceptions
            error_msg = "\n".join(sys.exc_info())
            print(f'[ERROR] Unable to retrieve reference from URL {ftp_link}\nThe following error occurred: {error_msg}')
            sys.exit(1)
    
    outfasta, outbin, outmeta = standardize_sequence_names(name = genome_name, 
                                                            reference_id = reference_id, 
                                                            localfile = local_filename,
                                                            fasta_dir = fasta_dir,
                                                            metadata_dir = metadata_dir)

    print(f'Downloaded {genome_name}')
    return (outfasta, outbin, outmeta)

def standardize_sequence_names(name, reference_id, localfile, fasta_dir, metadata_dir):
    fasta_name  = os.path.join(f'{fasta_dir}/{name}.fna')
    bin_name    = os.path.join(f'{metadata_dir}/{name}.bin.tsv')
    meta_name   = os.path.join(f'{metadata_dir}/{name}.metadata.tsv')

    outfasta = open(fasta_name,'w')
    outbin = open(bin_name,'w')
    outmeta = open(meta_name,'w')
    
    outmeta.write('AssemblyID\tGenomeName\tNewSeqHeaders\tOriginalSeqHeaders\tSeqLen\tNum_Ns\n')

    try: 
        with gzip.open(localfile, "rt") as handle:
            for index, record in enumerate(SeqIO.parse(handle, "fasta")):
                header = f'{name}_Node_{index}'
                # Metadata
                # Cols (tab-delim): AssemblyID, GenomeName, NewSeqHeaders, OriginalSeqHeaders, SeqLen, Num_Ns
                outmeta.write(f'{reference_id}\t{name}\t{header}\t{record.id}\t{len(record.seq)}\t{record.seq.count("N")}\n')
                outbin.write(f'{header}\t{name}\n') # contig_name\tstrain_name
                outfasta.write(f'>{header}\n{record.seq}\n') # reference fasta
    except:
        error_msg = "\n".join(sys.exc_info())
        print(f'[ERROR] Could not process fasta file for {name} from {localfile}\n{error_msg}')
    else:
        os.remove(localfile)
    finally:
        outfasta.close()
        outbin.close()
        outmeta.close()

    return (fasta_name, bin_name, meta_name)

def concatenate_files(output_file, file_list, has_header, cleanup):
    header_saved=False
    with open(output_file,'w') as outfile:
        for filename in file_list:
            with open(filename) as infile:
                if has_header:
                    header = next(infile)
                    if not header_saved:
                        outfile.write(header)
                        header_saved = True
                for line in infile:
                    outfile.write(line)
            if cleanup:
                os.remove(filename)
    return

def generate_db_summary(genome_metadata_file, db_metadata_file):
    metadata = pd.read_table(genome_metadata_file).groupby('GenomeName').agg(
        Num_contigs = ('NewSeqHeaders', 'nunique'),
        Genome_size = ('SeqLen', 'sum'), 
        Total_Ns = ('Num_Ns', 'sum')
    )
    metadata['N_frac'] = metadata['Total_Ns']/metadata['Genome_size']
    metadata.sort_values('N_frac', ascending=False).to_csv(db_metadata_file)
    # genome_qc
    return

def genome_qc():
    pass

if __name__ == '__main__':
    # db_ref_file example (see sheet "v2"): 
    # https://docs.google.com/spreadsheets/d/1rZUfisKKmPfLkZP1eEsNFyCKUyMm6nHQQAhNxnJJYXE/edit#gid=0
    # This file needs to have the following 2 columns: "Strain_Name" and "FTP_link"
    args = usage()
    db_ref_file = args['db_ref_file']
    db_name = args['db_name']
    cleanup = not args['no_cleanup']
    output_dir = os.path.join(db_name,'db')
    fasta_dir = os.path.join(db_name,'fasta')
    metadata_dir = os.path.join(db_name,'metadata')

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(fasta_dir, exist_ok=True)
    os.makedirs(metadata_dir, exist_ok=True)

    genomes_df = pd.DataFrame()
    db_files = pd.DataFrame()

    # Read input
    genomes_df = process_raw_seedfile(db_ref_file)
    
    # Download Genomes
    try:
        db_files = genomes_df.apply(lambda row: download_genome(genome_name = row['Strain_Name'], 
                                                                reference_id = row['refIDs'], 
                                                                ftp_link = row['FTP_link'],
                                                                fasta_dir = fasta_dir,
                                                                metadata_dir = metadata_dir),
                                                                axis=1)
        db_files=db_files.apply(pd.Series).dropna()
        db_files.columns=['fastafiles', 'binfiles', 'metafiles']
    except:
        # error_msg = "\n".join(sys.exc_info())
        # print(f'[FATAL] {error_msg}')
        sys.exit(1)
    
    # Concatenate data
    concatenate_files(output_file = f'{output_dir}/{db_name}.fna', file_list = db_files['fastafiles'], has_header=False, cleanup = False)
    concatenate_files(output_file = f'{output_dir}/{db_name}.nm_mates_bin.tsv', file_list = db_files['binfiles'], has_header=False, cleanup = cleanup)
    concatenate_files(output_file = f'{output_dir}/{db_name}.per_genome_metadata.tsv', file_list = db_files['metafiles'], has_header=True, cleanup = cleanup)
    
    generate_db_summary(f'{output_dir}/{db_name}.per_genome_metadata.tsv', f'{output_dir}/{db_name}.db_metadata.tsv')

    if cleanup:
        os.rmdir(metadata_dir)

    print('Done!')
