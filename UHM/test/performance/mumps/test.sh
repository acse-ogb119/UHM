#!/bin/bash

# Sym
echo '****** Test for all decomposition method for one element *******'

../mumpstest 1 3 1 ../../uhmfile/toy_1.uhm
../mumpstest 1 3 1 ../../uhmfile/toy_1.uhm
../mumpstest 1 3 1 ../../uhmfile/toy_1.uhm
../mumpstest 1 3 1 ../../uhmfile/toy_1.uhm

../mumpstest 2 3 1 ../../uhmfile/toy_1.uhm
../mumpstest 2 3 1 ../../uhmfile/toy_1.uhm
../mumpstest 2 3 1 ../../uhmfile/toy_1.uhm
../mumpstest 2 3 1 ../../uhmfile/toy_1.uhm

echo '****** Test for all decomposition method for two element *******'

../mumpstest 1 3 1 ../../uhmfile/toy_2.uhm
../mumpstest 1 3 1 ../../uhmfile/toy_2.uhm
../mumpstest 1 3 1 ../../uhmfile/toy_2.uhm
../mumpstest 1 3 1 ../../uhmfile/toy_2.uhm

../mumpstest 2 3 1 ../../uhmfile/toy_2.uhm
../mumpstest 2 3 1 ../../uhmfile/toy_2.uhm
../mumpstest 2 3 1 ../../uhmfile/toy_2.uhm
../mumpstest 2 3 1 ../../uhmfile/toy_2.uhm

echo '****** Test for all decomposition method for two element *******'

../mumpstest 1 3 1 ../../uhmfile/toy_3.uhm
../mumpstest 1 3 1 ../../uhmfile/toy_3.uhm
../mumpstest 1 3 1 ../../uhmfile/toy_3.uhm
../mumpstest 1 3 1 ../../uhmfile/toy_3.uhm

../mumpstest 2 3 1 ../../uhmfile/toy_3.uhm
../mumpstest 2 3 1 ../../uhmfile/toy_3.uhm
../mumpstest 2 3 1 ../../uhmfile/toy_3.uhm
../mumpstest 2 3 1 ../../uhmfile/toy_3.uhm

echo '****** Test for all decomposition method for toy mesh *******'

../mumpstest 1 3 1 ../../uhmfile/toy_mesh.uhm
../mumpstest 1 3 1 ../../uhmfile/toy_mesh.uhm
../mumpstest 1 3 1 ../../uhmfile/toy_mesh.uhm
../mumpstest 1 3 1 ../../uhmfile/toy_mesh.uhm

../mumpstest 2 3 1 ../../uhmfile/toy_mesh.uhm
../mumpstest 2 3 1 ../../uhmfile/toy_mesh.uhm
../mumpstest 2 3 1 ../../uhmfile/toy_mesh.uhm
../mumpstest 2 3 1 ../../uhmfile/toy_mesh.uhm


# Non-Sym
echo '****** Test for all decomposition method for one element *******'

../mumpstest 1 3 0 ../../uhmfile/toy_1.uhm
../mumpstest 1 3 0 ../../uhmfile/toy_1.uhm
../mumpstest 1 3 0 ../../uhmfile/toy_1.uhm
../mumpstest 1 3 0 ../../uhmfile/toy_1.uhm

../mumpstest 2 3 0 ../../uhmfile/toy_1.uhm
../mumpstest 2 3 0 ../../uhmfile/toy_1.uhm
../mumpstest 2 3 0 ../../uhmfile/toy_1.uhm
../mumpstest 2 3 0 ../../uhmfile/toy_1.uhm

echo '****** Test for all decomposition method for two element *******'

../mumpstest 1 3 0 ../../uhmfile/toy_2.uhm
../mumpstest 1 3 0 ../../uhmfile/toy_2.uhm
../mumpstest 1 3 0 ../../uhmfile/toy_2.uhm
../mumpstest 1 3 0 ../../uhmfile/toy_2.uhm

../mumpstest 2 3 0 ../../uhmfile/toy_2.uhm
../mumpstest 2 3 0 ../../uhmfile/toy_2.uhm
../mumpstest 2 3 0 ../../uhmfile/toy_2.uhm
../mumpstest 2 3 0 ../../uhmfile/toy_2.uhm

echo '****** Test for all decomposition method for two element *******'

../mumpstest 1 3 0 ../../uhmfile/toy_3.uhm
../mumpstest 1 3 0 ../../uhmfile/toy_3.uhm
../mumpstest 1 3 0 ../../uhmfile/toy_3.uhm
../mumpstest 1 3 0 ../../uhmfile/toy_3.uhm

../mumpstest 2 3 0 ../../uhmfile/toy_3.uhm
../mumpstest 2 3 0 ../../uhmfile/toy_3.uhm
../mumpstest 2 3 0 ../../uhmfile/toy_3.uhm
../mumpstest 2 3 0 ../../uhmfile/toy_3.uhm

echo '****** Test for all decomposition method for toy mesh *******'

../mumpstest 1 3 0 ../../uhmfile/toy_mesh.uhm
../mumpstest 1 3 0 ../../uhmfile/toy_mesh.uhm
../mumpstest 1 3 0 ../../uhmfile/toy_mesh.uhm
../mumpstest 1 3 0 ../../uhmfile/toy_mesh.uhm

../mumpstest 2 3 0 ../../uhmfile/toy_mesh.uhm
../mumpstest 2 3 0 ../../uhmfile/toy_mesh.uhm
../mumpstest 2 3 0 ../../uhmfile/toy_mesh.uhm
../mumpstest 2 3 0 ../../uhmfile/toy_mesh.uhm
