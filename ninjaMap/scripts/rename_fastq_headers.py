#!/usr/bin/env python3
from Bio import SeqIO
import sys 
import os

def is_number(n):
    try:
        float(n)   # Type-casting the string to `float`.
                   # If string is not a valid `float`, 
                   # it'll raise `ValueError` exception
    except ValueError:
        return False
    return True

def parse_header(seq_record_id):
    header_number = None
    header_name = None
    delimiter = None
    # is there a space in the name?
    if " " in seq_record_id:
        # split at space
        delimiter = " "
    # is there a / in the name?
    elif "/" in seq_record_id:
        # split at '/'
        delimiter = "/"
    # if split at _ is the last element a number?
    elif "_" in seq_record_id:
        delimiter = "_"

    if delimiter is not None:
        header_pieces = seq_record_id.split(delimiter)
        header_number = header_pieces.pop()
        if is_number(header_number):
            header_name = delimiter.join(header_pieces)
    
    return (header_name, header_number, delimiter)

def main(seq_file, prefix, outdir):
    # grabbing the file and the name 
    labels = os.path.basename(seq_file).split(".")
    output_fastq = f'{outdir}/{prefix}.fastq'
    try:
        out_fastq = open(output_fastq, 'a+')
        with open(seq_file,"r") as in_handle:
            for seq_record in SeqIO.parse(in_handle,"fastq"):
                current_header_name, current_header_number, delimiter = parse_header(seq_record.id)
                seq_record.id = f'{labels[0]}{delimiter}{current_header_number}'
                SeqIO.write(seq_record, out_fastq, "fastq")
        out_fastq.close()
    except:
        print(sys.exc_info())
        os.remove(output_fastq)


if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2],sys.argv[3])