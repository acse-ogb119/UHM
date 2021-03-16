#!/bin/bash

# Sym
echo '****** Test for all decomposition method for one element *******'

../pardisotest 1 3 1 ../../uhmfile/toy_1.uhm
../pardisotest 1 3 1 ../../uhmfile/toy_1.uhm
../pardisotest 1 3 1 ../../uhmfile/toy_1.uhm
../pardisotest 1 3 1 ../../uhmfile/toy_1.uhm

../pardisotest 2 3 1 ../../uhmfile/toy_1.uhm
../pardisotest 2 3 1 ../../uhmfile/toy_1.uhm
../pardisotest 2 3 1 ../../uhmfile/toy_1.uhm
../pardisotest 2 3 1 ../../uhmfile/toy_1.uhm

echo '****** Test for all decomposition method for two element *******'

../pardisotest 1 3 1 ../../uhmfile/toy_2.uhm
../pardisotest 1 3 1 ../../uhmfile/toy_2.uhm
../pardisotest 1 3 1 ../../uhmfile/toy_2.uhm
../pardisotest 1 3 1 ../../uhmfile/toy_2.uhm

../pardisotest 2 3 1 ../../uhmfile/toy_2.uhm
../pardisotest 2 3 1 ../../uhmfile/toy_2.uhm
../pardisotest 2 3 1 ../../uhmfile/toy_2.uhm
../pardisotest 2 3 1 ../../uhmfile/toy_2.uhm

echo '****** Test for all decomposition method for two element *******'

../pardisotest 1 3 1 ../../uhmfile/toy_3.uhm
../pardisotest 1 3 1 ../../uhmfile/toy_3.uhm
../pardisotest 1 3 1 ../../uhmfile/toy_3.uhm
../pardisotest 1 3 1 ../../uhmfile/toy_3.uhm

../pardisotest 2 3 1 ../../uhmfile/toy_3.uhm
../pardisotest 2 3 1 ../../uhmfile/toy_3.uhm
../pardisotest 2 3 1 ../../uhmfile/toy_3.uhm
../pardisotest 2 3 1 ../../uhmfile/toy_3.uhm

echo '****** Test for all decomposition method for toy mesh *******'

../pardisotest 1 3 1 ../../uhmfile/toy_mesh.uhm
../pardisotest 1 3 1 ../../uhmfile/toy_mesh.uhm
../pardisotest 1 3 1 ../../uhmfile/toy_mesh.uhm
../pardisotest 1 3 1 ../../uhmfile/toy_mesh.uhm

../pardisotest 2 3 1 ../../uhmfile/toy_mesh.uhm
../pardisotest 2 3 1 ../../uhmfile/toy_mesh.uhm
../pardisotest 2 3 1 ../../uhmfile/toy_mesh.uhm
../pardisotest 2 3 1 ../../uhmfile/toy_mesh.uhm


# Non-Sym
echo '****** Test for all decomposition method for one element *******'

../pardisotest 1 3 0 ../../uhmfile/toy_1.uhm
../pardisotest 1 3 0 ../../uhmfile/toy_1.uhm
../pardisotest 1 3 0 ../../uhmfile/toy_1.uhm
../pardisotest 1 3 0 ../../uhmfile/toy_1.uhm

../pardisotest 2 3 0 ../../uhmfile/toy_1.uhm
../pardisotest 2 3 0 ../../uhmfile/toy_1.uhm
../pardisotest 2 3 0 ../../uhmfile/toy_1.uhm
../pardisotest 2 3 0 ../../uhmfile/toy_1.uhm

echo '****** Test for all decomposition method for two element *******'

../pardisotest 1 3 0 ../../uhmfile/toy_2.uhm
../pardisotest 1 3 0 ../../uhmfile/toy_2.uhm
../pardisotest 1 3 0 ../../uhmfile/toy_2.uhm
../pardisotest 1 3 0 ../../uhmfile/toy_2.uhm

../pardisotest 2 3 0 ../../uhmfile/toy_2.uhm
../pardisotest 2 3 0 ../../uhmfile/toy_2.uhm
../pardisotest 2 3 0 ../../uhmfile/toy_2.uhm
../pardisotest 2 3 0 ../../uhmfile/toy_2.uhm

echo '****** Test for all decomposition method for two element *******'

../pardisotest 1 3 0 ../../uhmfile/toy_3.uhm
../pardisotest 1 3 0 ../../uhmfile/toy_3.uhm
../pardisotest 1 3 0 ../../uhmfile/toy_3.uhm
../pardisotest 1 3 0 ../../uhmfile/toy_3.uhm

../pardisotest 2 3 0 ../../uhmfile/toy_3.uhm
../pardisotest 2 3 0 ../../uhmfile/toy_3.uhm
../pardisotest 2 3 0 ../../uhmfile/toy_3.uhm
../pardisotest 2 3 0 ../../uhmfile/toy_3.uhm

echo '****** Test for all decomposition method for toy mesh *******'

../pardisotest 1 3 0 ../../uhmfile/toy_mesh.uhm
../pardisotest 1 3 0 ../../uhmfile/toy_mesh.uhm
../pardisotest 1 3 0 ../../uhmfile/toy_mesh.uhm
../pardisotest 1 3 0 ../../uhmfile/toy_mesh.uhm

../pardisotest 2 3 0 ../../uhmfile/toy_mesh.uhm
../pardisotest 2 3 0 ../../uhmfile/toy_mesh.uhm
../pardisotest 2 3 0 ../../uhmfile/toy_mesh.uhm
../pardisotest 2 3 0 ../../uhmfile/toy_mesh.uhm
