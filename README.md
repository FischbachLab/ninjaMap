# NinjaMap

ninjaMap is a software tool to calculate strain abundance for a given microbial database.

This tool runs in two steps, ninjaIndex and ninjaMap. It will accept a directory of your reference genomes (one genome per file). It calculate the uniqueness of the genome in the database along with other contigs related metadata, and return a binmap file along with a concatenated fasta file of your references.

## Requirements

1. Nextflow (https://www.nextflow.io/)
2. Python 3.7
3. numpy
4. scipy
5. numba
6. pysam
7. SeqIO
9. samtools
10. bowtie2
11. bbmap


## Pipeline overview

Step 1. The ninjaIndex pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [ART] - Generate synthetic short reads for each genome
* [Bowtie2] - align reads to all reference genomes
* [ninjaIndex] - generate ninja index for a given synthetic community

Step 2. The ninjaMap pipeline accurately quantify a strain with abundance.


## Quick usage

### 1. Generate a ninjaIndex file on aws

The input of the ninjaIndex is a list of genome files in fasta format.

```bash
nextflow run ./nf-core-ninjaindex/main.nf --genomes 's3://bucket/input/*.fna' --outdir 's3://bucket/output/' -profile aws
```
OR

```{bash}
aws batch submit-job \
    --profile maf \
    --job-name nf-ninjaindex \
    --job-queue default-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=s3://nextflow-pipelines/ninjaindex,\
"--genomes","'s3://dev-scratch/ReferenceDBs/NinjaMap/Index/12Com/fasta/*.fna'",\
"--outdir","s3://genomics-workflow-core/Results/NinjaIndex/12Com/db"
```

### 2. Run ninjaMap with an existing ninjaIndex using 16 threads:

The main input of the ninjaMap is a binmap file generated from the 1st step, and a sorted BAM file and its indexed bam.bai file must be present in same directory.

```pyhton
./ninjaMap/scripts/ninjaMap.py -threads 16 -bam sample_sorted.bam -bin binmap.tsv -prefix mycommunity
```
A wrapper script ninjaMap/ninjaMap_index.sh is provided to run ninjaMap on aws via batch.

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

## ninjaMap Full usage

```python
python  ./ninjaMap/scripts/ninjaMap.py --help
Description:
This script will calculate the abundance of a strain in a defined microbial community.
Usage: ninjaMap.py -bam sorted.bam -bin binmap.tsv -prefix my_community

optional arguments:
  -h, --help          show this help message and exit
  -bam BAMFILE        sorted bam file and its indexed bam.bai file must be present in same directory.
  -bin BINMAP         tab-delimited file with Col1= contig name and Col2=Bin/Strain name
  -outdir OUTDIR      output directory
  -prefix PREFIX      output prefix
  -threads THREADS    number of threads available for this job and subprocesses
  -debug              save intermediate false positives bam file
  -truth TRUTH        If using debug, please provide one strain name that you would like to track.
  -mbq MIN_BASE_QUAL  minimum read base quality to consider for coverage calculations.
```
Output files for each sample
====================

The output files are organized into 4 folders.

## bowtie2 folder

The alignment file of all input reads aligned the defined community database in the bam format

## Logs folder

The running logs of various scripts

## ninjaMap folder

1. **\*.ninjaMap.abundance.csv**: this file shows the statistics of the abundance, coverage and depth of each strain in the defined community

+ Strain_Name: strain name<br>
+ Read_Fraction: the abundance in the defined community in percentage<br>
+ Percent_Coverage: the average coverage per strain in percentage<br>
+ Coverage_Depth: the average coverage depth<br>

2. **\*.ninjaMap.read_stats.csv**: this file shows the statistics of input reads

+ File_Name: sample name <br>
+ Reads_Aligned: the number of aligned reads<br>
+ Reads_wPerfect_Aln: the number of perfectly aligned reads<br>
+ Reads_wSingular_Votes: the number of reads voted as singular<br>
+ Reads_wEscrowed_Votes: the number of reads voted as escrow<br>
+ Discarded_Reads_w_Perfect_Aln: the number of discarded perfectly aligned reads

3. **\*.ninjaMap.strain_stats.csv**: this file shows the various statistics of each strains

4. **\*.ninjaMap.votes.csv.gz**: the statistics of reads voting (singular or escrow)

## Stats folder
1. **adapter_trimming_stats_per_ref.txt**: this file shows the statistics of adapter trimming

2. **read_accounting.csv**: this file shows the statistics shows the total number of reads, the number of reads after trimming and the number of aligned reads


Aggregated output files for each study
====================

The aggregated output files are organized into 6 files.

1. **\*.covDepth.csv**: this file shows the average coverage depth per strain by samples
2. **\*.host_contaminants.csv**: this file shows the detected host contaminants (Human or Mouse) by samples if the unalignment rate is over 5%
3. **\*.long.csv**: this is the long format of three files (\*.readFraction.csv, \*.covDepth.csv and \*.percCoverage.csv)
4. **\*.percCoverage.csv**: this file shows the average coverage per strain in percentage by samples
5. **\*.reads_stats.csv**: this file shows the reads statistics in read numbers by samples
6. **\*.readFraction.csv**: this file shows the abundance in the defined community in percentage by samples


## Citation

To be added

## License

[GNU GPL] (https://www.gnu.org/licenses/gpl-3.0.html)

## Questions / Concerns

- [Sunit Jain](https://www.sunitjain.com/)
- Xiandong Meng (xdmeng at stanford.edu)
- PI: Michael Fischbach ( fischbach at fischbachgroup.org )
