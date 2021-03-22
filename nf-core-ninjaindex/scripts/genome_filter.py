#!/usr/bin/env python

import sys
from Bio import SeqIO

def list_headers():
    """
    Return a set containing the headers presented in a genome,
    line by line, starting with ">"
    """

    headers = set()

    with open(sys.argv[1], 'r') as target_fasta:
        records = SeqIO.parse(target_fasta, 'fasta')
        for record in records:
        #for line in fi:
           # line = line.strip()
           # print (line)
            headers.add(str(record.description).replace(">", ""))
    
    return headers

def filter():
    """
    Writes a filtered file containing only the sequences with headers NOT
    present in a set of headers
    """

    headers = list_headers()

    with open(sys.argv[2], 'r') as all_fasta, open(sys.argv[3], 'w') as filtered_fasta:
        records = SeqIO.parse(all_fasta, 'fasta')
        for record in records:
            #print (record.description)
            if record.description not in headers:
                SeqIO.write(record, filtered_fasta, 'fasta')

if __name__ == '__main__':
    if len(sys.argv) != 4 :
       print ('USAGE: ' + sys.argv[0] + ' target_genome.fa all_genomes.fa filtered_genomes.fa')
    else: 
       filter()

