#!/bin/bash

JAVA_BIN_PATH='java'
SEP=$([ $(java -h |& grep ";" | wc -l ) == 1 ] && echo ";" || echo ":")
JAVA_CLASS_PATH="../MSynRec/bin${SEP}../MSynRec/lib/commons-math3-3.3.0.jar"

num=2
spe="Sea Lamprey,Japanese Lamprey"
chr="198,191"
mul="4,4"
K=18
L=0.1

WORK_DIR='./vertebrate'
INPUT_FILE="${WORK_DIR}/probSynDataAll.txt"
OUTPUT_DIR="${WORK_DIR}/probsyn_CVB_${K}_${L}"
mkdir -p ${OUTPUT_DIR}


${JAVA_BIN_PATH} \
-DLOGSUMEXP=true \
-cp ${JAVA_CLASS_PATH} analysis.macrosynteny.MSynRec \
-i "${INPUT_FILE}" \
-o "${OUTPUT_DIR}" \
-T ${num} \
-S "${spe}" \
-C ${chr} \
-D ${mul} \
-A 0.1 \
-B 0.1 \
-K ${K} \
-l ${L} \
-E true,true \
-M CVB0 \
