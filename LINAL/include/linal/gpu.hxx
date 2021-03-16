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
#ifndef LINAL_FLAT_GPU_HXX
#define LINAL_FLAT_GPU_HXX

namespace linal {
  typedef class  Matrix_*     Matrix;
  typedef class  Flat_GPU_*   Flat_GPU;
  typedef struct FLA_Obj_view FLA_Obj;

  // -------------------------------------------------------------- 
  // ** Global variable handling
  extern void init_gpu          ( int n_thread, int n_gpu );
  extern void init_gpu          ( int n_thread, int n_gpu, int n_flat_gpu );
  extern void finalize_gpu      ();
  extern void set_device_gpu    ( int device );
  extern int  get_device_gpu    ( int thread );
  extern int  get_n_thread      (); 
  extern int  get_n_gpu         ();
  extern int  get_n_flat_gpu    ();
  extern Flat_GPU_& get_flat_gpu( int thread, int index );
  extern void create_flat_gpu   ( int type, int m, int n );
  extern void free_flat_gpu     ();
  
  // -------------------------------------------------------------- 
  // ** Flat GPU Matrix
  /*!
    Flat_GPU_ class only stored in the GPU memory
  */
  class Flat_GPU_ : public Matrix_ {
  private:
  protected:
    int  device, buffer_gpu_created;
    __device__ void *buffer_gpu;

  public:

    /*!
      Default constructor.
    */
    Flat_GPU_() : Matrix_(),
                  device(0), 
                  buffer_gpu_created(0), 
                  buffer_gpu(NULL) {  }

    /*!
      Wrapping constructor - It is thought not being created.
      \param obj [in] Input FLA_Obj
    */
    Flat_GPU_(FLA_Obj obj) : Matrix_(),
                             device(0), 
                             buffer_gpu_created(0), 
                             buffer_gpu(NULL) { 
      LINAL_ERROR(obj.base != NULL &&
                  obj.base->elemtype == FLA_SCALAR,
                  ">> FLA_Obj is null or MATRIX type");
      // if linal wrap the FLA_Obj, linal does not control
      // FLA_Obj buffer allocation. 
      this->fla = obj; 
      this->created = 0;
    }

    /*! 
      Copy constructor.
      \param obj [in] Input Flat_GPU_
    */
    Flat_GPU_(const Flat_GPU_& obj) : Matrix_(obj),
                                      device(0), 
                                      buffer_gpu_created(0), 
                                      buffer_gpu(NULL)  { 
      this->buffer_gpu = obj.buffer_gpu;
      this->created =0;
    }

    /*!
      Destructor.
    */
    virtual ~Flat_GPU_() { this->free(); }

    /*! Set device ID which gpu is bound
     */
    inline void set_device(int device) { 
      LINAL_ERROR(device < get_n_gpu(),
                  ">> Device ID is greater than the number of devices");
      this->device = device; 
    }
    
    /*! Get device ID
     */
    inline int  get_device() { 
      return this->device; 
    }

    /*!
      Create m by n Flat_ matrix.
      \param type [in] LINAL_INT, LINAL_REAL, or LINAL_COMPLEX
      \param m    [in] number of rows
      \param n    [in] number of columns
    */
    inline void create(int type, int m, int n) {
      this->create_without_buffer( type, m, n );
      this->create_buffer();
    }

    /*!
      Create m by n Flat_ matrix without bufer.
      \param type [in] LINAL_INT, LINAL_REAL, or LINAL_COMPLEX
      \param m    [in] number of rows
      \param n    [in] number of columns
    */
    inline void create_without_buffer(int type, int m, int n) {
      // either m or n is 0, do not create
      if (!m || !n) return;

      LINAL_ERROR( check_all_scalar_type(type),
                  ">> Datatype is wrong (should be INT, REAL, or COMPLEX");
      LINAL_ERROR(!this->is_created(),
                  ">> This is already created");

      FLA_Obj_create_without_buffer(type, m, n, &this->fla);
      this->created = 1;
    }
    
    /*! 
      Free the matrix. If buffer was allocated, free it.
    */
    inline void free() {
      if (this->is_created()) {
        FLA_Obj_free_without_buffer(&this->fla);
        this->created = 0;
      }
      if (this->buffer_gpu_created)
        this->free_buffer();
    }

    /*!
      Create buffer on the gpu
    */
    inline void create_buffer() {
      LINAL_ERROR( !is_base_null(), 
                   ">> Base object is null");

      cudaSetDevice( this->get_device() );
      
      cublasStatus status;

      // create m by n matrix
      status = cublasAlloc( this->get_m()*this->get_n(),
                            this->get_data_size(),
                            (void**)&this->buffer_gpu );

      LINAL_ERROR( status == CUBLAS_STATUS_SUCCESS,
                   ">> CUBLAS fail to alloc" );
      this->buffer_gpu_created = 1;
    }

    /*!
      Free buffer on the gpu
    */
    inline void free_buffer() {

      cublasFree( this->buffer_gpu );

      this->buffer_gpu         = NULL;
      this->buffer_gpu_created = 0;
    }

    /*!
      Execute cublasGetMatrix ( memory moves from GPU to CPU )
      \param obj [in] FLA_Obj will read GPU buffer
    */
    inline void get_matrix_to(FLA_Obj obj) {

      cublasStatus status;

      // Read buffer on GPU to CPU
      status = cublasGetMatrix( this->get_m(),
                                this->get_n(),
                                this->get_data_size(),
                                this->buffer_gpu,
                                this->get_m(),
                                FLA_Obj_buffer_at_view( obj ),
                                FLA_Obj_col_stride( obj ) );

      LINAL_ERROR( status == CUBLAS_STATUS_SUCCESS,
                   ">> CUBLAS fail GetMatrix" );
    }

    /*!
      Execute cublasSetMatrix ( memory moves from CPU to GPU )
      \param obj [in] FLA_Obj whose contenst will be written in GPU
    */
    inline void set_matrix_from(FLA_Obj obj) {

      cublasStatus status;

      // Write the CPU memory into GPU
      status = cublasSetMatrix( FLA_Obj_length( obj ),
                                FLA_Obj_width( obj ),
                                FLA_Obj_datatype_size( FLA_Obj_datatype( obj ) ),
                                FLA_Obj_buffer_at_view( obj ), 
                                FLA_Obj_col_stride( obj ),
                                this->buffer_gpu,
                                this->get_m() );

      LINAL_ERROR( status == CUBLAS_STATUS_SUCCESS,
                   ">> CUBLAS fail SetMatrix" );
    }

    /*!
      Partition the matrix on the GPU.
     */
    inline void extract(Flat_GPU_ &A, int m, int n) {
      this->extract(A, m, n, 0, 0);
    }

    inline void extract(Flat_GPU_ &A, int m, int n, int offm, int offn) {
      if (!(offm+m <= this->get_m() &&
	    offn+n <= this->get_n())) {

	std::ostringstream msg;
	msg << ">> Out of range" << std::endl
	    << "  range (" << this->get_m() << "," << this->get_n() << ")" 
	    << std::endl
	    << "  offset(" << offm << "," << offn << ")" << std::endl
	    << "  dim   (" << m    << "," << n    << ")" << std::endl;
	
	LINAL_ERROR( false, msg.str().c_str() );
      }

      A.fla.m       = m;
      A.fla.n       = n;
      A.fla.offm    = this->get_offm() + offm;
      A.fla.offn    = this->get_offn() + offn;
      A.fla.base    = this->get_fla().base;

      // share the gpu buffer
      A.buffer_gpu         = this->get_buffer();
      A.buffer_gpu_created = 0;
      return;
    }

    inline void  set_buffer(void *buf_ptr) { this->buffer_gpu = buf_ptr; }
    inline void* get_buffer()              { return (double*)this->buffer_gpu; }
    
    inline bool operator==(Flat_GPU_ &b) {
      return (this->get_m()         == b.get_m() &&
              this->get_n()         == b.get_n() &&
              this->get_data_type() == b.get_data_type());
    }
    inline bool operator!=(Flat_GPU_ &b) { return !(*this == b); }

    inline void disp() { this->disp(stdout, (char*)"- Flat GPU -"); }
    inline void disp(char *title) { this->disp(stdout, title); }
    inline void disp(FILE* stream) { this->disp(stream, (char*)"- Flat GPU -"); }
    inline void disp(FILE* stream, char *title){ 
      if (this->get_m() && this->get_n()) {
        fprintf(stream, "%s", title);
	switch (this->get_data_type()) {
	case LINAL_INT:            fprintf(stream, "INT     TYPE\n");  break;
	case LINAL_SINGLE_REAL:   
	case LINAL_DOUBLE_REAL:    fprintf(stream, "REAL    TYPE\n");  break;
	case LINAL_SINGLE_COMPLEX: 
	case LINAL_DOUBLE_COMPLEX: fprintf(stream, "COMPLEX TYPE\n");  break;
	}
        fprintf(stream, "(%4d, %4d )::( %4d, %4d )\n",
                this->get_offm(), this->get_offn(),
                this->get_m(), this->get_n());

      }
    }
  };
}

#endif
