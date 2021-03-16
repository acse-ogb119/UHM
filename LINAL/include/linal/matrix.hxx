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
#ifndef LINAL_MATRIX_HXX
#define LINAL_MATRIX_HXX

namespace linal {
  typedef class  Matrix_*     Matrix;
  typedef struct FLA_Obj_view FLA_Obj;

  // -------------------------------------------------------------- 
  // ** Matrix
  /*! 
    Matrix_ class is a simple wrapper of FLA_Obj. This class has
    some extentional member functions to safe free or malloc.
   */
  class Matrix_ {
  private:
  protected:
    int  created; //!< Flag indicate matrix is created
    FLA_Obj fla;  //!< FLA_Obj
  public:
    /*!
      Default constructor initialize the FLA_Obj. 
      Initially FLA_Obj has NULL base object.
    */
    Matrix_() { 
      this->created  = 0; 
      this->fla.m    = 0; 
      this->fla.n    = 0;
      this->fla.offm = 0;
      this->fla.offn = 0;
      this->fla.base = NULL;
    }

    /*!
      Copy constructor does not create FLA_Obj but share the FLA_Obj with
      an input argument.
      \param b [in] Input matrix
    */
    Matrix_(const Matrix_&b) { 
      this->fla     = b.fla; 
      this->created = 0; 
    }

    /*! 
      Destructor calls free().
      \sa free
     */
    virtual ~Matrix_() { this->free(); }
    
    /*! 
      Free the FLA_Obj only if it is created
      \sa ~Matrix_
     */
    inline void free() {  
      if (this->is_created()) {
          FLA_Obj_free(&this->fla);
          this->created = 0;
      }
    }

    /*!
      Create buffer based on the datatype and dimension of FLA_Obj.
    */
    inline void create_buffer() { 
      FLA_Obj_create_buffer( 0, 0, &this->fla ); 
      FLA_Obj_set_to_scalar( FLA_ZERO, this->fla );
    }
    
    /*!
      Free buffer if it has buffer allocated.
    */
    inline void free_buffer() {
      FLA_Obj_free_buffer( &this->fla );
    }

    // ** query functions
    inline bool is_created() { return this->created; }
    inline bool is_buffer_null() { 
      return (is_base_null() || FLA_Obj_buffer_is_null(this->fla)); 
    } 
    inline bool is_base_null() { return (this->fla.base == NULL); }
    inline bool is_hier() { 
      return (FLA_Obj_elemtype(this->fla) == LINAL_MATRIX); 
    }
    inline bool is_real_datatype() {
      return (FLA_Obj_datatype(this->fla) == LINAL_REAL);
    }
    inline bool is_complex_datatype() {
      return (FLA_Obj_datatype(this->fla) == LINAL_COMPLEX);
    }

    // ** these does not have set interface
    inline int get_data_type()  { 
      return FLA_Obj_datatype(this->fla);
    }
    inline int get_data_size() {
      return FLA_Obj_datatype_size(this->get_data_type());
    }
    inline int get_buffer_size() {
      return (this->get_data_size() * this->get_m() * this->get_n());
    }
    /*!
      Get column stride
    */
    inline int get_cs()    { return this->fla.base->cs; }
    /*!
      Get row stride
    */
    inline int get_rs()    { return this->fla.base->rs; }

    // ** these can be used for extration of matrices 
    //    so it has set and get interfaces
    inline void set_m(int m)          { this->fla.m    = m; }
    inline void set_n(int n)          { this->fla.n    = n; }
    inline void set_offm(int offm)    { this->fla.offm = offm; }
    inline void set_offn(int offn)    { this->fla.offn = offn; }
    inline void set_fla(FLA_Obj obj)  { this->fla      = obj; }
    
    inline int  get_m()               { return this->fla.m; }
    inline int  get_n()               { return this->fla.n; }
    inline int  get_offm()            { return this->fla.offm; }
    inline int  get_offn()            { return this->fla.offn; }

    /*!
      Return the FLA_Obj to use FLAME
    */
    inline FLA_Obj& get_fla()         { return this->fla; }

    /*! 
      Operator returns the FLA_Obj.
    */
    inline FLA_Obj& operator ~()      { return this->get_fla(); }
    /*!
      Copy the FLA_Obj in input matrix b. This share the buffer with b.
      \param b [in] Input matrix
    */
    inline Matrix_& operator=(const Matrix_ &b) { 
      this->fla = b.fla;
      this->created= 0;
      return *this;
    }
  };
}
#endif
