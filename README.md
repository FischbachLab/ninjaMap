# NinjaMap

ninjaMap is a software tool to calculate strain abundance for a given microbial database.

This tool runs in two steps, ninjaIndex and ninjaMap. It will accept a directory of your reference genomes (one genome per file). It calculate the uniqueness of the genome in the database along with other contigs related metadata, and return a binmap file along with a concatenated fasta file of your references.

## Requirements

1. Nextflow (https://www.nextflow.io/)
2. Python 3.7
3. samtools 1.9
4. pysam
5. bowtie2
6. bbmap


## Pipeline overview

Step 1. The ninjaIndex pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [ART] - Generate synthetic short reads for each genome
* [Bowtie2] - align reads to all reference genomes
* [ninjaIndex] - generate ninja index for a given synthetic community

Step 2. The ninjaMap pipeline accurately quantify a strain with abundance.


## Quick usage

### Generate a ninjaIndex file on aws
```bash

nextflow run ./nf-core-ninjaindex/main.nf --genomes 's3://bucket/input/*.fna' --outdir 's3://bucket/output/' -profile aws
```

### Run ninjaMap with an existing ninjaIndex using 16 threads:
```pyhton
./ninjaMap/scripts/ninjaMap.py -threads 16 -bam sample_sorted.bam -bin binmap.tsv -prefix mycommunity
```

## Docker
The generic command to run a ninjaMap docker container:

```bash
docker container run \
    -v /host/path/to/indata/:/input_data/ \
    -v /host/path/to/outdata/:/output_data/ \
    fischbachlab/ninjamap \
    python ./scripts/ninjaMap.py \
    -bin /input_data/binmap.tsv \
    -bam /input_data/input.sorted.bam \
    -outdir /output_data/summary \
    -prefix mycommunity
```

## Citation

To be added

## License

[GNU GPL] (https://www.gnu.org/licenses/gpl-3.0.html)

## Questions / Concerns

- [Sunit Jain](https://www.sunitjain.com/)
- Xiandong Meng (xdmeng at stanford.edu)
- PI: Michael Fischbach ( fischbach at fischbachgroup.org )
