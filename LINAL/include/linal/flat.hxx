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
#ifndef LINAL_FLAT_HXX
#define LINAL_FLAT_HXX

namespace linal {
  typedef class  Matrix_*     Matrix;
  typedef class  Flat_*       Flat;
  typedef struct FLA_Obj_view FLA_Obj;
  
  // -------------------------------------------------------------- 
  // ** Flat Matrix
  /*!
    Flat_ class contains scalars operation. Flat_ can be thought as wrapper 
    for the FLA_Obj.
  */
  class Flat_ : public Matrix_ {
  private:
  protected:
  public:
    /*!
      Default constructor.
    */
    Flat_() : Matrix_() { }
    /*!
      Wrapping constructor.
      \param obj [in] Input FLA_Obj
    */
    Flat_(FLA_Obj obj) : Matrix_() { 
      LINAL_ERROR(obj.base != NULL &&
                  obj.base->elemtype == FLA_SCALAR,
                  ">> FLA_Obj is null or MATRIX type");

      // if linal wrap the FLA_Obj, linal does not control
      // FLA_Obj buffer allocation. 
      this->fla = obj; 
      this->created = 0;
    }
    /*!
      Copy constrcutor. 
      \param obj [in] Input Flat_ matrix object.
    */
    Flat_(const Flat_& obj) : Matrix_(obj) { }
    /*!
      Destructor.
    */
    virtual ~Flat_() { this->free(); }

    /*!
      Create m by n Flat_ matrix.
      \param type [in] LINAL_INT, LINAL_REAL, or LINAL_COMPLEX
      \param m    [in] number of rows
      \param n    [in] number of columns
     */
    inline void create(int type, int m, int n) {
      // either m or n is 0, do not create
      if (!m || !n) return;
      
      LINAL_ERROR( check_all_scalar_type(type),
                  ">> Datatype is wrong (should be INT, REAL, or COMPLEX");
      LINAL_ERROR(!this->is_created(),
                  ">> This is already created");

      FLA_Obj_create(type, m, n, 0, 0, &this->fla);
      // FLA_Obj_set_to_scalar(FLA_ZERO, this->fla);
      this->created = 1;
    }

    /*!
      Create Flat_ matrix which has same size of input matrix Flat_ b.
      \param trans [in] LINAL_NO_TRANSPOSE, LINAL_TRANSPOSE
      \param b     [in] Input matrix
     */
    inline void create(int trans, Flat_ b) {
      // if b is empty, do not create the matrix
      if (!b.get_m() || !b.get_n()) return;

      LINAL_ERROR( check_all_scalar_type( b.get_data_type() ),
                  ">> Datatype is wrong (should be INT, REAL, or COMPLEX");
      LINAL_ERROR(!this->is_created(),
                  ">> This is already created");

      // switch m and n according to trans
      switch (trans) {
      case LINAL_NO_TRANSPOSE:
	FLA_Obj_create(b.get_data_type(), b.get_m(), b.get_n(), 0,0, &this->fla);
	break;
      case LINAL_TRANSPOSE:
	FLA_Obj_create(b.get_data_type(), b.get_n(), b.get_m(), 0,0, &this->fla);
	break;
      }

      // FLA_Obj_set_to_scalar(FLA_ZERO, this->fla);
      this->created = 1;
    }

    /*!
      Symbolic matrix is created without buffer. Even if users try to free 
      this object, it does not double free the buffer.
      \param type [in] LINAL_INT, LINAL_REAL, or LINAL_COMPLEX
      \param m    [in] number of rows
      \param n    [in] number of columns        
     */
    inline void create_without_buffer(int type, int m, int n) {
      // either m or n is 0, do not create
      if (!m || !n) return;

      LINAL_ERROR( check_all_scalar_type(type) ,
                  ">> Datatype is wrong (should be INT, REAL, or COMPLEX");
      LINAL_ERROR(!this->is_created(),
                  ">> This is already created");

      FLA_Obj_create_without_buffer(type, m, n, &this->fla);
      this->created = 1;
    }

    /*!
      Wrapping. Flat_ matrix is not created. It supports the access by Flat_
      to FLA_Obj.
      \param obj [in] Input FLA_Obj
    */
    inline void wrap(FLA_Obj &obj)  {
      LINAL_ERROR(this->created == 0,
                  ">> Flat_ is already created");
      LINAL_ERROR(obj.base != NULL &&
                  obj.base->elemtype == FLA_SCALAR,
                  ">> FLA_Obj is null or MATRIX type");

      // if linal wrap the FLA_Obj, linal does not control
      // FLA_Obj buffer allocation.
      this->fla = obj;
      this->created = 0;
    }

    /*!
      Extract block matrix of m by n from this to Flat_ A.
      
      \param A    [out] Matrix object which has block matrix
      \param m    [in]  number of rows
      \param n    [in]  number of columns     
      \param offm [in]  number of offset in rows
      \param offn [in]  number of offset in columns
     */
    inline void extract(Flat_ &A, int m, int n, int offm, int offn) {
      if (!(offm+m <= this->get_m() &&
	    offn+n <= this->get_n())) {

	std::ostringstream msg;
	msg << ">> Out of range" << std::endl
	    << "  range (" << this->get_m() << "," << this->get_n() 
	    << std::endl
	    << "  offset(" << offm << "," << offn << ")" << std::endl
	    << "  dim   (" << m    << "," << n    << ")" << std::endl;
	
	LINAL_ERROR( false, msg.str().c_str() );
      }

      A.fla.m    = m;
      A.fla.n    = n;
      A.fla.offm = this->get_offm() + offm;
      A.fla.offn = this->get_offn() + offn;
      A.fla.base = this->get_fla().base;
      return;
    }

    /*!
      Linal only consider double type for the buffer.
     */
    inline double* get_buffer() {
      return (double*)FLA_Obj_buffer_at_view(this->get_fla());
    }
    
    /*!
      Overloaded operator to see this matrix equal to input matrix. 
      It compare the type and dimension.
      \param b [in] Input matrix
    */
    inline bool operator==(Flat_ &b) {
      return (this->get_m()         == b.get_m() &&
	      this->get_n()         == b.get_n() &&
	      this->get_data_type() == b.get_data_type());
    }

    /*!
      Overloaded operator to see this matrix is NOT equal to input matrix.
      Simply negate the operator "=="
      \param b [in] Input matrix
    */
    inline bool operator!=(Flat_ &b) { return !(*this == b); }
    

    /*!
      Return the reference of member indicated by row and col indices
      in double type matrix.
      \param row [in] row index
      \param col [in] column index
     */
    inline double& operator() (unsigned int row, unsigned int col) {
      // for real
      if (!(row<this->get_m() && col<this->get_n())) {
	std::ostringstream msg;
        msg << ">> Out of range" << std::endl
            << "  range (" << this->get_m() << "," << this->get_n() << ")" 
            << std::endl
            << "  index (" << row  << "," << col    << ")" 
            << std::endl;

        LINAL_ERROR( false, msg.str().c_str() );
      }
      return (this->get_buffer())[col*this->get_cs()+row];
    }
    
    /*!
      Return the reference of member indicated by row and col indices
      in complex type matrix.
      \param row [in] row index
      \param col [in] column index
      \param idx [in] 0 - real, 1 - imag
    */
    inline double& operator() (unsigned int row, unsigned int col,
			       unsigned int idx) {
      // for complex idx 0 : real , idx 1 : imag
      if (!(row<this->get_m() && col<this->get_n())) {
	std::ostringstream msg;
        msg << ">> Out of range" << std::endl
            << "  range (" << this->get_m() << "," << this->get_n() << ")"
            << std::endl
            << "  index (" << row    << "," << col    << ")" << std::endl;

        LINAL_ERROR( false, msg.str().c_str() );
      }
      return (this->get_buffer())[2*col*this->get_cs()+2*row+idx];
    }
    inline void disp() { this->disp(stdout, (char*)"- Flat -"); }
    inline void disp(char *title) { this->disp(stdout, title); }
    inline void disp(FILE* stream) { this->disp(stream, (char*)"- Flat -"); }
    inline void disp(FILE* stream, char *title){ 
      if (this->get_m() && this->get_n()) {
	switch (this->get_data_type()) {
	case LINAL_INT:
	  FLA_Obj_fshow(stream, title, this->get_fla(), (char*)"% d",
                        (char*)"----------");
	  break;
	case LINAL_SINGLE_REAL:
	case LINAL_DOUBLE_REAL:
	  FLA_Obj_fshow(stream, title, this->get_fla(), (char*)"% E",
			(char*)"----------");
          break;
	case LINAL_SINGLE_COMPLEX:
	case LINAL_DOUBLE_COMPLEX:
	  FLA_Obj_fshow(stream, title, this->get_fla(), (char*)"( % E , % E )",
			(char*)"----------");
          break;
	}
      }
    }
  };
}

#endif
