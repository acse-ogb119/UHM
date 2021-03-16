#!/bin/bash

echo '****** Test for all decomposition method for one element *******'

../uhmtest 1 3 256 ../../uhmfile/toy_1.uhm
../uhmtest 1 3 256 ../../uhmfile/toy_1.uhm
../uhmtest 1 3 256 ../../uhmfile/toy_1.uhm
../uhmtest 1 3 256 ../../uhmfile/toy_1.uhm

../uhmtest 2 3 256 ../../uhmfile/toy_1.uhm
../uhmtest 2 3 256 ../../uhmfile/toy_1.uhm
../uhmtest 2 3 256 ../../uhmfile/toy_1.uhm
../uhmtest 2 3 256 ../../uhmfile/toy_1.uhm

echo '****** Test for all decomposition method for two element *******'

../uhmtest 1 3 256 ../../uhmfile/toy_2.uhm
../uhmtest 1 3 256 ../../uhmfile/toy_2.uhm
../uhmtest 1 3 256 ../../uhmfile/toy_2.uhm
../uhmtest 1 3 256 ../../uhmfile/toy_2.uhm

../uhmtest 2 3 256 ../../uhmfile/toy_2.uhm
../uhmtest 2 3 256 ../../uhmfile/toy_2.uhm
../uhmtest 2 3 256 ../../uhmfile/toy_2.uhm
../uhmtest 2 3 256 ../../uhmfile/toy_2.uhm

echo '****** Test for all decomposition method for two element *******'

../uhmtest 1 3 256 ../../uhmfile/toy_3.uhm
../uhmtest 1 3 256 ../../uhmfile/toy_3.uhm
../uhmtest 1 3 256 ../../uhmfile/toy_3.uhm
../uhmtest 1 3 256 ../../uhmfile/toy_3.uhm

../uhmtest 2 3 256 ../../uhmfile/toy_3.uhm
../uhmtest 2 3 256 ../../uhmfile/toy_3.uhm
../uhmtest 2 3 256 ../../uhmfile/toy_3.uhm
../uhmtest 2 3 256 ../../uhmfile/toy_3.uhm

echo '****** Test for all decomposition method for toy mesh *******'

../uhmtest 1 3 256 ../../uhmfile/toy_mesh.uhm
../uhmtest 1 3 256 ../../uhmfile/toy_mesh.uhm
../uhmtest 1 3 256 ../../uhmfile/toy_mesh.uhm
../uhmtest 1 3 256 ../../uhmfile/toy_mesh.uhm

../uhmtest 2 3 256 ../../uhmfile/toy_mesh.uhm
../uhmtest 2 3 256 ../../uhmfile/toy_mesh.uhm
../uhmtest 2 3 256 ../../uhmfile/toy_mesh.uhm
../uhmtest 2 3 256 ../../uhmfile/toy_mesh.uhm
