/*
  Copyright © 2011, Kyungjoo Kim
  All rights reserved.
  
  This file is part of LINAL.
  
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

  3. Neither the name of the owner nor the names of its contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.
*/
#include "linal/common.hxx"
#include "linal/const.hxx"
#include "linal/util.hxx"
#include "linal/matrix.hxx"
#include "linal/flat.hxx"
#include "linal/hier.hxx"

#ifdef LINAL_GPU_ENABLE

#include "cublas.h"
#include "linal/gpu.hxx"

namespace linal {

  static int number_of_thread = 0, number_of_gpu = 0, number_of_flat_gpu = 0;
  static std::vector< Flat_GPU_ > bin;

  // --------------------------------------------------------------
  // ** GPU init
  void init_gpu(int n_thread, int n_gpu) {
    init_gpu(n_thread, n_gpu, 3);
  }
  void init_gpu(int n_thread, int n_gpu, int n_flat_gpu) {

    // ** cublas init
    cublasInit();

    // ** initialization of variables
    bin.clear(); 
    number_of_thread   = n_thread;
    number_of_gpu      = n_gpu;
    number_of_flat_gpu = n_flat_gpu;

    int n_size_of_bin  = n_thread * n_flat_gpu;

    // ** add gpu container to bin
    for (int j=0;j<n_thread;++j) {
      for (int i=0;i<n_flat_gpu;++i) {
        Flat_GPU_ in;
        in.set_device( get_device_gpu( j ) );
        bin.push_back( in );
      }
    }
  }

  void finalize_gpu() {
    // ** clear all variables
    bin.clear();    
    number_of_gpu      = 0;
    number_of_flat_gpu = 0;

    // ** cublas shut down
    cublasShutdown();
  }

  int  get_device_gpu( int thread ) {
    return (thread%get_n_gpu());
  }

  void set_device_gpu( int device ) {
    cudaSetDevice( device );
  }

  // this thread number does not have to be same as actual thread number
  int  get_n_thread()   { return number_of_thread; }
  int  get_n_gpu()      { return number_of_gpu; }
  int  get_n_flat_gpu() { return number_of_flat_gpu; }
  
  Flat_GPU_& get_flat_gpu( int thread, int index ) {
    return bin.at(thread*number_of_flat_gpu + index);
  } 
  
  void create_flat_gpu( int type, int m, int n ) {
    for (int i=0;i<bin.size();++i) 
      bin.at(i).create( type, m, n );
  }

  void free_flat_gpu() {
    for (int i=0;i<bin.size();++i) 
      bin.at(i).free();
  }
}

#endif
