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
#ifndef LINAL_HIER_HXX
#define LINAL_HIER_HXX

namespace linal {
  typedef class  Matrix_*     Matrix;
  typedef class  Flat_*       Flat;
  typedef class  Hier_*       Hier;
  typedef struct FLA_Obj_view FLA_Obj;
  
  // -------------------------------------------------------------- 
  // ** Hier Matrix
  /*!
    Hier_ class contains hierarhical matrix operations.
    This class only uses 1-level of blocking.
  */
  class Hier_ : public Matrix_ {
  private:
  protected:
  public:
    /*!
      Default constructor.
    */
    Hier_() : Matrix_() { }

    /*!
      Wrapping constructor.
      \param obj [in] Input FLA_Obj
    */
    Hier_(FLA_Obj obj) : Matrix_() {
      LINAL_ERROR(obj.base != NULL &&
                  obj.base->elemtype == FLA_MATRIX,
                  ">> FLA_Obj is null or SCALAR type");

      // if linal wrap the FLA_Obj, linal does not control 
      // FLA_Obj buffer allocation.
      this->fla = obj;
      this->created = 0;
    }

    /*!
      Copy constructor.
      \param obj [in] Input Hier_ matrix object.
    */
    Hier_(const Hier_& obj) : Matrix_(obj) { }

    /*!
      Destructor.
    */
    virtual ~Hier_() { this->free(); }

    /*!
      Create m by n Hier_ matrix. Its buffer contains matrices. 
      \param type [in] LINAL_REAL, or LINAL_COMPLEX
      \param m    [in] number of rows
      \param n    [in] number of columns   
     */
    inline void create(int type, int m, int n) {
      // either m or n is 0, do not create
      if (!m || !n) return;

      // Hier_ matrix does not support integer type.
      LINAL_ERROR( check_all_scalar_type(type) ,
                  ">> Datatype is wrong (should be REAL, or COMPLEX");
      LINAL_ERROR(!this->is_created(),
                  ">> This is already created");
      
      FLA_Obj_create_ext(type, FLA_MATRIX, m, n, 0, 0, 0, 0, &this->fla);
      this->created = 1;
    }

    /*!
      Create Hier_ matrix with block size (mb, nb) for the matrix size (m,n)
      \param type [in] LINAL_REAL, or LINAL_COMPLEX
      \param m    [in] number of rows
      \param n    [in] number of columns
      \param mb   [in] rowise blocksize
      \param nb   [in] columnwise blocksize 
    */
    inline void create(int type, int m, int n, int mb, int nb) {
      // either m or n is 0, do not create
      if (!m || !n) return;

      LINAL_ERROR( check_all_scalar_type(type), 
                  ">> Datatype is wrong (should be REAL, or COMPLEX");
      LINAL_ERROR(!this->is_created(),
                  ">> This is already created");
      
      std::vector< int > offm, offn;
      int offs, mm,nn;
      FLA_Obj *buf;

      // calculate the size of Hier_ matrix
      mm = m/mb + (m%mb > 0);
      nn = n/nb + (n%nb > 0);

      // create                                                     
      FLA_Obj_create_ext(type, FLA_MATRIX,
			 mm, nn, 0, 0, 0, 0,
			 &this->fla);

      buf = (FLA_Obj*)FLA_Obj_buffer_at_view(this->fla);

      // calculate offset values
      offs = 0;
      for (int i=0;i<mm;++i) { offm.push_back(offs); offs += mb; }
      offm.push_back(m);

      offs = 0;
      for (int j=0;j<nn;++j) { offn.push_back(offs); offs += nb; }
      offn.push_back(n);

      // assign block matrices to the matrix objects in the Hier_ buffer
      for (int j=0;j<nn;j++) {
	for (int i=0;i<mm;i++) {
	  FLA_Obj *part = &buf[ j*mm + i ];

	  part->offm = offm.at(i);
	  part->offn = offn.at(j);

	  part->m    = offm.at(i+1) - offm.at(i);
	  part->n    = offn.at(j+1) - offn.at(j);

	  part->base = NULL;
	}
      }
      
      // this is created
      this->created = 1;
    }

    /*!
      Create Hier_ matrix with block size (mb, nb) for the Flat_ matrix.
      \param F  [in] Input matrix
      \param mb [in] rowise blocksize
      \param nb [in] columnwise blocksize
    */
    inline void create(Flat_ F, int mb, int nb) {
      if (!F.get_m() || !F.get_n()) return;
      this->create(F.get_data_type(), F.get_m(), F.get_n(), mb, nb);
      for (int j=0;j<this->get_n();++j) 
	for (int i=0;i<this->get_m();++i) 
	  (*this)(i,j).base = (~F).base;
    }

    /*!
      Create Hier_ matrix with merging four 2d partitioned Hier_ submatrices 
      \param ATL [in] Top Left matrix
      \param ABL [in] Bottom Left matrx
      \param ATR [in] Top Right matrix
      \param ABR [in] Bottom Right matrix
    */
    inline void create_h4(Hier_ ATL, Hier_ ATR,
                          Hier_ ABL, Hier_ ABR) {

      this->create(ATL.get_data_type(), 
                   ATL.get_m()+ABL.get_m(), ATL.get_n()+ATR.get_n());
      
      int offm = ATL.get_m(), offn = ATL.get_n();

      // ATL
      for (int j=0;j<ATL.get_n();j++) 
	for (int i=0;i<ATL.get_m();i++) 
	  (*this)(i,j) = ATL(i,j);

      // ATR
      for (int j=0;j<ATR.get_n();j++) 
	for (int i=0;i<ATR.get_m();i++) 
	  (*this)(i,offn+j) = ATR(i,j);

      // ABL
      for (int j=0;j<ABL.get_n();j++) 
	for (int i=0;i<ABL.get_m();i++) 
	  (*this)(offm+i,j) = ABL(i,j);
      
      // ABR
      for (int j=0;j<ABR.get_n();j++) 
	for (int i=0;i<ABR.get_m();i++) 
	  (*this)(offm+i,offn+j) = ABR(i,j);
    }

    /*!
      Create Hier_ matrix with merging vertically partitioned Hier_ submatrices 
      \param ATL [in] Top Left matrix
      \param ABL [in] Bottom Left matrx
    */
    inline void create_h2_ver(Hier_ ATL, 
                              Hier_ ABL) {

      this->create(ATL.get_data_type(), 
                   ATL.get_m()+ABL.get_m(), ATL.get_n());
      
      int offm = ATL.get_m();

      // ATL
      for (int j=0;j<ATL.get_n();j++) 
	for (int i=0;i<ATL.get_m();i++) 
	  (*this)(i,j) = ATL(i,j);

      // ABL
      for (int j=0;j<ABL.get_n();j++) 
	for (int i=0;i<ABL.get_m();i++) 
	  (*this)(offm+i,j) = ABL(i,j);
    }


    /*!
      Create Hier_ matrix with merging horizontally partitioned Hier_ submatrices 
      \param ATL [in] Top Left matrix
      \param ATR [in] Top Right matrix
    */
    inline void create_h2_hor(Hier_ ATL, Hier_ ATR) {
                              
      this->create(ATL.get_data_type(), 
                   ATL.get_m(), ATL.get_n()+ATR.get_n());
      
      int offn = ATL.get_n();

      // ATL
      for (int j=0;j<ATL.get_n();j++) 
	for (int i=0;i<ATL.get_m();i++) 
	  (*this)(i,j) = ATL(i,j);

      // ATR
      for (int j=0;j<ATR.get_n();j++) 
	for (int i=0;i<ATR.get_m();i++) 
	  (*this)(i,offn+j) = ATR(i,j);
    }

    /*!
      Wrapping. Hier_ matrix is not created. 
      It supports the access by Hier_ to FLA_Obj.
      \param obj [in] Input FLA_Obj
    */
    inline void wrap(FLA_Obj &obj)  {
      LINAL_ERROR(this->created == 0,
                  ">> Hier_ is already created");
      LINAL_ERROR(obj.base != NULL &&
                  obj.base->elemtype == FLA_MATRIX,
                  ">> FLA_Obj is null or MATRIX type");

      // if linal wrap the FLA_Obj, linal does not control
      // FLA_Obj buffer allocation.
      this->fla = obj;
      this->created = 0;
    }
    /*!
      Extract block matrix of m by n from this to Hier_ A.
      
      Remark - This corresponds to partitioning in FLAME.
     */
    inline void extract(Hier_ &A, int m, int n, int offm, int offn) {
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

      A.fla.m    = m; 
      A.fla.n    = n;
      A.fla.offm = this->get_offm() + offm;
      A.fla.offn = this->get_offn() + offn;
      A.fla.base = this->get_fla().base;
      return;
    }

    /*!
      Convert this Hier_ matrix into Flat_ matrix object.
     */
    inline Flat_ flat() {
      FLA_Obj obj;
      int m=0, n=0, offm=0, offn=0;

      if (this->get_n()>0) 
        for (int i=0;i<this->get_m();++i) 
          m+= Flat_((*this)(i,0)).get_m();

      if (this->get_m()>0)
        for (int j=0;j<this->get_n();++j)
          n+= Flat_((*this)(0,j)).get_n();

      if (m && n) { 
        obj = (*this)(0,0);
        obj.m    = m; obj.n    = n; 
      } else {
        obj.offm = 0; obj.offn = 0;
        obj.m    = m; obj.n    = n; 
      }
      return Flat_(obj);
    }

    /*!
      Return the buffer of Hier_ matrix.
     */
    inline FLA_Obj* get_buffer() {
      return (FLA_Obj*)FLA_Obj_buffer_at_view(this->get_fla());
    }
    /*!
      Display the matrix to stdout
    */
    inline void disp() { disp(stdout); }
    /*!
      Display the matrix to stream
    */
    inline void disp(FILE* stream) {
      fprintf(stream, " - Hier -\n");
      for (int i=0;i<this->get_m();++i) {
	for (int j=0;j<this->get_n();++j) {
          if ( (*this)(i,j).base == NULL ) {
            fprintf( stream, "  (%4d, %4d )::( %4d, %4d ) ",
                     (*this)(i,j).offm,
                     (*this)(i,j).offn,
                     FLA_Obj_length((*this)(i,j)),
                     FLA_Obj_width((*this)(i,j)) );
          } else {
            fprintf( stream, "x (%4d, %4d )::( %4d, %4d ) ",
                     (*this)(i,j).offm,
                     (*this)(i,j).offn,
                     FLA_Obj_length((*this)(i,j)),
                     FLA_Obj_width((*this)(i,j)) );
          }
	}
	fprintf(stream, "\n");
      }
    }

    /*!
      Overloaded operator to see this matrix equal to input matrix.
      It compare the type and dimension.
      \param b [in] Input matrix
    */
    inline bool operator==(Hier_ &b) {
      return (this->get_m()         == b.get_m() &&
	      this->get_n()         == b.get_n() &&
	      this->get_data_type() == b.get_data_type());
    }

    /*!
      Overloaded operator to see this matri
      Overloaded operator to see this matrix is NOT equal to input matrix.
      Simply negate the operator "=="
      \param b [in] Input matrix  
     */
    inline bool operator!=(Hier_ &b) { return !(*this == b); }

    /*!
      Return the reference of block matrix indicated by row and col indicies.
      \param row [in] row index
      \param col [in] column index
     */
    inline FLA_Obj& operator() (unsigned int row, unsigned int col) {
      if (!(row<this->get_m() && col<this->get_n())) {
	std::ostringstream msg;
        msg << ">> Out of range" << std::endl
            << "  range (" << this->get_m() << "," << this->get_n() << ")"
            << std::endl
            << "  index (" << row    << "," << col    << ")" << std::endl;

        LINAL_ERROR( false, msg.str().c_str() );
      }
      return (this->get_buffer())[col*this->get_cs()+row];
    }
  };
}

#endif
