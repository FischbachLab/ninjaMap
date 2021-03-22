#!/usr/bin/env python

import os
import sys
import s3fs
import numpy as np
import pandas as pd

input_dir = sys.argv[1].rstrip('/')
dbname = sys.argv[2]

# input_dir = 's3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Index/SCv1_1_20200422'
output_dir = input_dir
# output_dir = '/Users/sunit.jain/Research/NinjaMap/Database/SCv1_1_20200422'
report_csv = f'{dbname}.report.csv'

checkm_summary_file = f'{input_dir}/qc/checkm/qc.tsv'
gtdbtk_summary_file = f'{input_dir}/qc/GTDBtk/gtdbtk.bac120.summary.tsv'
genome_summary_file = f'{input_dir}/qc/stats/fasta_stats.txt'
gene_summary_file  = f'{input_dir}/qc/orfs/num_genes.txt'

database_basic_stat = f'{input_dir}/db/{dbname}.db_metadata.tsv'

def checkm_score(row):
    # From GTDB:
    # Completeness - (5 * Contamination)
    # Should be less than or equal to 50
    return row['Completeness'] - (5 * row['Contamination'])

## GTDB
def classify_method_stats(row):
    row['bin_name'] = row['user_genome']
    if row['classification_method']=="ANI/Placement":
        row['ani']=row['fastani_ani']
        row['aln_frac']=row['fastani_af']
        row['closest_ref']=row['fastani_reference']
    elif row['classification_method']=="Placement":
        row['ani']=row['closest_placement_ani']
        row['aln_frac']=row['closest_placement_af']
        row['closest_ref']=row['closest_placement_reference']
    else:
        row['ani']=None
        row['aln_frac']=None
        row['closest_ref']=None
    return row

## CheckM
def parse_checkm_output(checkm_summary):
    df = pd.read_table(checkm_summary).rename(columns={"Bin Id": "bin_name", "Strain heterogeneity":"Strain_heterogeneity"})
    df['checkm_score'] = df.apply(checkm_score, axis = 1)
    keep_cols = ['bin_name','Completeness', 'Contamination', 'Strain_heterogeneity', 'checkm_score']
    return df[keep_cols].sort_values('checkm_score', ascending = True).reset_index(drop = True)

# GTDB
def parse_gtdb_outputs(gtdb_output):
    df = pd.read_table(gtdb_output)
    df = df.apply(classify_method_stats, axis = 1)
    keep_col = ['bin_name','classification','classification_method','closest_ref', 'ani', 'aln_frac','aa_percent','red_value', 'note', 'warnings']
    return df[keep_col].set_index('bin_name')

# Genes
def parse_gene_stats(genes_file):
    return pd.read_table(genes_file, names=['bin_name','num_genes'], 
                        index_col='bin_name')

# SeqKit
def parse_seqkit_stats(stats_file):
    df = pd.read_table(stats_file)
    df['bin_name'] = df['file'].apply(lambda x: os.path.basename(os.path.splitext(x)[0]))
    keep_cols=['bin_name','min_len','avg_len','max_len']
    return df[keep_cols].set_index('bin_name')

if __name__ == '__main__':
    print('Parsing CheckM output ...')
    checkm_df = parse_checkm_output(checkm_summary_file).set_index('bin_name')
    print('Parsing Genome stats ...')
    db_stats_df = pd.read_csv(database_basic_stat).rename(columns={'GenomeName':'bin_name'}).set_index('bin_name')
    print('Parsing GTDBtk output ...')
    gtdb_df = parse_gtdb_outputs(gtdbtk_summary_file)
    print('Parsing Seqkit stats ...')
    bin_stats_df = parse_seqkit_stats(genome_summary_file)
    print('Parsing Genes output ...')
    genes_df = parse_gene_stats(gene_summary_file)

    ## Aggregate
    print('Aggregating data ...')
    bins_df = pd.concat([db_stats_df,
                        bin_stats_df,
                        genes_df,
                        gtdb_df, 
                        checkm_df],
                        axis=1, sort=False).sort_values('checkm_score', ascending = True).reset_index().rename(columns={'index':'bin_name'})
    print('Writing output ...')
    bins_df.to_csv(os.path.join(output_dir,report_csv), index = False)
    print('Done.')
