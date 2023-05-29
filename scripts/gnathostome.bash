#!/bin/bash

JAVA_BIN_PATH='java'
JAVA_CLASS_PATH='../MSynDup/bin:../MSynRec/lib/commons-math3-3.6.1.jar'

K=10
L="0.1"

WORK_DIR='./gnathostome'
SUPP_DATA_DIR="../SuppData1"
SEG_DIR="${SUPP_DATA_DIR}/segments/seg"
PAR_DIR="${SUPP_DATA_DIR}/paralogs"
ORT_DIR="${SUPP_DATA_DIR}/orthologs"
INPUT_SEG_GROUP_FILE="${WORK_DIR}/psm_rec.txt"
OUTPUT_FILE="${WORK_DIR}/psm_rec_ALL.txt"

speciesList=("Human" "Mouse" "Dog" "Opossum" "Chicken" "Turkey" "Zebra Finch" "Spotted Gar" "Elephant Shark")
swcodeList=(HUMAN MOUSE CANFA MONDO CHICK MELGA POEGU LEPOC CALMI)

args="${JAVA_BIN_PATH} -cp ${JAVA_CLASS_PATH} reconstruction.multispecies.MSynDup -T ALL -F 3"
args="${args} -g \"${INPUT_SEG_GROUP_FILE}\" -o \"${OUTPUT_FILE}\""

for x in ${!speciesList[@]}
do
 SPX=${speciesList[$x]}
 SWX=${swcodeList[$x]}
 args="${args} -s \"${SPX},${SEG_DIR}/${SWX}.txt,${PAR_DIR}/${SWX}-${SWX}.txt\""
done

for x in ${!speciesList[@]}
do
 SPX=${speciesList[$x]}
 SWX=${swcodeList[$x]}
 for y in ${!speciesList[@]}
 do
  [[ $x -ge $y ]] && continue
  SPY=${speciesList[$y]}
  SWY=${swcodeList[$y]}
  args="${args} -S \"${SPX},${SPY},${ORT_DIR}/${SWX}-${SWY}.txt\""
 done
done

eval ${args}
