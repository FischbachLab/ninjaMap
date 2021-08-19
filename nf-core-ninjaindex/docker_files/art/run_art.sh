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
# RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
PAIRED_FASTQ="${LOCAL_OUTPUT}/paired_fastq"
LOCAL_FASTA="${OUTPUTDIR}/${FASTA}"

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${PAIRED_FASTQ}"
#aws s3 cp ${S3FASTA} ${LOCAL_FASTA}
cp ${S3FASTA} ${LOCAL_FASTA}

FWD="${PAIRED_FASTQ}/${PREFIX}.R1.fastq.gz"
REV="${PAIRED_FASTQ}/${PREFIX}.R2.fastq.gz"

COV_FOLD=${2:-10}

# ART Command
/mnt/art_bin_MountRainier/art_illumina \
-na \
-ef \
-ss HS25 \
-i "${LOCAL_FASTA}" \
-p \
-l 150 \
-f ${COV_FOLD} \
-m 500 \
-s 10 \
-rs 1712 \
-o "${OUTPUTDIR}/${PREFIX}" &> "${LOG_DIR}.log"

# Get error free reads
#java -jar $PICARD SamToFastq I=${PREFIX}_errFree.sam F=${PREFIX}1.fq F2=${PREFIX}2.fq
# picard SamToFastq I="${OUTPUTDIR}/${PREFIX}_errFree.sam" F="${OUTPUTDIR}/${PREFIX}_R1.fq" F2="${OUTPUTDIR}/${PREFIX}_R2.fq"
picard SamToFastq -I "${OUTPUTDIR}/${PREFIX}_errFree.sam" -F "${OUTPUTDIR}/${PREFIX}_R1.fq" -F2 "${OUTPUTDIR}/${PREFIX}_R2.fq"

#compress fastq files
gzip < "${OUTPUTDIR}/${PREFIX}_R1.fq" > $FWD
gzip < "${OUTPUTDIR}/${PREFIX}_R2.fq" > $REV

# Sync
#aws s3 sync ${LOCAL_OUTPUT} ${S3OUTPUTPATH}
#aws s3 cp ${LOCAL_OUTPUT} ${S3OUTPUTPATH}
