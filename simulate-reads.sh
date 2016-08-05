#! /bin/bash
# Use ART to Simulate reads
#
GENOMES=$1/*
OUT_DIR=$2
ART=$3
SEED=12345
COVERAGES=( 0.25 0.5 1 5 )
LOGS=${OUT_DIR}/logs

mkdir -p ${LOGS}
for g in ${GENOMES}; do
    FILE=${g##*/}
    GENOME=${FILE%.*}
    echo "Simulating reads for ${GENOME}"
    for cov in "${COVERAGES[@]}"; do
        echo "    ${cov}x coverage..."
        mkdir -p ${OUT_DIR}/${cov}x
        OUTPUT=${OUT_DIR}/${cov}x/${GENOME}
        ${ART} -l 100 -f ${cov} -na -ss HS20 -rs ${SEED} -i ${g} -o ${OUTPUT} 1> ${LOGS}/${GENOME}-${cov}x.log 2>&1
    done
done
