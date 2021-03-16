#!/bin/bash

#
#   Copyright © 2011, Kyungjoo Kim
#   All rights reserved.
#  
#   This file is part of LINAL.
#  
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are met:
#
#   1. Redistributions of source code must retain the above copyright notice,
#     this list of conditions and the following disclaimer.
#
#   2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#
#   3. Neither the name of the owner nor the names of its contributors
#     may be used to endorse or promote products derived from this software
#     without specific prior written permission.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
#   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#   POSSIBILITY OF SUCH DAMAGE.
#


# Record the top path
TOP_PATH=$(pwd)
echo "**** TOP_PATH :: $TOP_PATH"

# Collect all testlist
TESTS=$(find . -name "testlist" -print)

# Make lib
if [ -f "../lib/liblinal.a" ]; then
    echo "**** liblinal.a exist"
else
    echo "**** liblinal.a NOT exist"
    cd ..
    make clean; make lib;
    cd $TOP_PATH
fi

# Make clean and all
echo "**** Build all TESTS"
for TEST in $TESTS
do 
    cd ${TEST%/*} 
    make all >> make.log
    cd $TOP_PATH 
done

# Counting
count=0;

for TEST in $TESTS
do
    while read LINE
    do
	count=$(( count + 1 ))
    done < $TEST
done
echo "**** Number of TESTS :: $count"

# Execute test code
n_test=0; n_fail=0;
for TEST in $TESTS
do
    cd ${TEST%/*}
    while read LINE
    do
  	./$LINE >> test.log
  	n_test=$(( n_test + 1 ))
 	if [ $? -ne 0 ]; then
 	    n_fail=$(( n_fail + 1 ))
 	    echo "[$n_test/$count] FAILED :: ${TEST%/*}/$LINE " 
 	else
 	    echo "[$n_test/$count] PASSED :: ${TEST%/*}/$LINE "
 	fi
     done < testlist
    cd $TOP_PATH 
done

echo "----------------------------------------------------"
echo "                    REPORT"
echo "----------------------------------------------------"
echo " Number of tests :: $count "
echo " FAILURES        :: $n_fail "
echo "----------------------------------------------------"

exit 0
