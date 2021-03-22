#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

set -e
set -u
set -o pipefail

S3FASTA="${1}"
#S3OUTPUTPATH="${2}"
#S3OUTPUTPATH=${S3OUTPUTPATH%/}
PIGZ_COMPRESSION_THREADS=${CORE_NUM:-4}

FASTA=$(basename -- "$S3FASTA")
PREFIX="${FASTA%.*}"

LOCAL=$(pwd)
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
RAW_FASTQ="${OUTPUTDIR}/interleaved_fastq"
PAIRED_FASTQ="${LOCAL_OUTPUT}/paired_fastq"
LOCAL_FASTA="${OUTPUTDIR}/${FASTA}"

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${PAIRED_FASTQ}"
#aws s3 cp ${S3FASTA} ${LOCAL_FASTA}
cp ${S3FASTA} ${LOCAL_FASTA}

FWD="${PAIRED_FASTQ}/${PREFIX}.R1.fastq.gz"
REV="${PAIRED_FASTQ}/${PREFIX}.R2.fastq.gz"

# NUM_SEQS=$(grep -c ">" ${LOCAL_FASTA})
COV_FOLD=${COV_FOLD:-10}

# Grinder Command
grinder \
    -cf ${COV_FOLD}\
    -rd 140 \
    -id 800 \
    -mo FR \
    -dc '-~*NX' \
    -md poly4 3e-3 3.3e-8 \
    -rs 1712 \
    -am uniform \
    -ql 30 10 \
    -fq 1 \
    -od "${RAW_FASTQ}" \
    -bn "${PREFIX}" \
    -rf "${LOCAL_FASTA}" &> "${LOG_DIR}.log"

# Deinterleave and compress
paste - - - - - - - - < "${RAW_FASTQ}/${PREFIX}-reads.fastq" | sed "s/@/@${PREFIX}_/g" | tee >(cut -f 1-4 | tr "\t" "\n" | pigz --best --processes ${PIGZ_COMPRESSION_THREADS} > $FWD) | cut -f 5-8 | tr "\t" "\n" | pigz --best --processes ${PIGZ_COMPRESSION_THREADS} > $REV

# Sync
#aws s3 sync ${LOCAL_OUTPUT} ${S3OUTPUTPATH}
#aws s3 cp ${LOCAL_OUTPUT} ${S3OUTPUTPATH}
