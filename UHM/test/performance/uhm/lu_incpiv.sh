#!/bin/bash

echo '****** Test for various block size *******'

#for i in 96 128 192 216 256 512 ; do \
#    ../uhmtest 12 4 $i ../../../uhmfile/sphere3.uhm
#done ;

echo '****** Test for various thread size *******'
for i in 1 2 4 6 8 10 12 14 16 18 20 22 24 ; do \
    ../uhmtest $i 4 192 ../../../uhmfile/sphere3.uhm
done ;

for i in 4 8 12 16 20 24 1; do \
    ../uhmtest $i 4 256 ../../../uhmfile/sphere4.uhm
done ;
