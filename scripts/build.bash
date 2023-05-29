#!/bin/bash

# init
mkdir ./tmp

# build msynrec
rm -r ../analysis/pvrec/bin
mkdir ../analysis/pvrec/bin

find ../analysis/pvrec/src -name "*.java" > tmp/sources.txt
javac @tmp/sources.txt -cp ../analysis/pvrec/lib/commons-math3-3.6.1.jar -d ../analysis/pvrec/bin

# build msyndup
rm -r ../analysis/pgrec/bin
mkdir ../analysis/pgrec/bin

find ../analysis/pgrec/src -name "*.java" > tmp/sources.txt
javac @tmp/sources.txt -cp ../analysis/pgrec/lib/commons-math3-3.6.1.jar -d ../analysis/pgrec/bin

# cleanup
rm -r ./tmp

echo "build finished"
