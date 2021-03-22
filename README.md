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

1. The ninjaIndex pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [ART] - Generate synthetic short reads for each genome
* [Bowtie2] - align reads to all reference genomes
* [ninjaIndex] - generate ninja index for a give synthetic community

2. The ninjaMap pipeline accurately quantify a strain with abundance.


## Quick usage

### Generate a ninjaIndex
```bash
nextflow run  main.nf --genomes 's3://inoput/*.fna' --outdir 's3://output/' -profile aws
```

### Run ninjaMap with an existing ninjaIndex:
```pyhton
ninjaMap.py -bam name_sorted.bam -bin contig_strain_assignments.tsv -prefix mycommunity
```

## Docker
The generic command to run a ninjaMap docker container:

```bash
docker container run \
    -v /host/path/to/indata/:/indata/ \
    -v /host/path/to/outdata/:/outdata/ \
    fischbachlab/ninjamap \
    python ninjaMap.py \
    -bin /indata/binmap.tsv \
    -bam /indata/input.sortedByCoord.bam \
    -outdir /outdata/summary \
    -prefix mycommunity
```

## Paper


## License

[GNU General Public License, version 3] (https://www.gnu.org/licenses/gpl-3.0.html)


## Questions / Concerns

- [Sunit Jain](microbiome.ninja) (dev)
- Xiandong Meng (xdmeng at stanford.edu)
- Michael Fischbach ( fischbach at fischbachgroup.org )
