#!/bin/bash

# symmetric case 
echo '****** Test for various thread size *******'
for i in 1 2 4 6 8 10 12 14 16 18 20 22 24 ; do \
    ../mumpstest $i 3 1 ../../../uhmfile/sphere3.uhm
done ;

for i in 1 4 8 12 16 20 24 ; do \
    ../mumpstest $i 3 1 ../../../uhmfile/sphere4.uhm
done ;

# un symmetric case
echo '****** Test for various thread size *******'
for i in 1 2 4 6 8 10 12 14 16 18 20 22 24 ; do \
    ../mumpstest $i 3 0 ../../../uhmfile/sphere3.uhm
done ;

for i in 1 4 8 12 16 20 24 ; do \
    ../mumpstest $i 2 0 ../../../uhmfile/sphere4.uhm
done ;


