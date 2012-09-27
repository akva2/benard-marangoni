/***************************************************************************
 *   Copyright (C) 2005 by Arne Morten Kvarving                            *
 *   spiff@mspiggy                                                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef MNL_MATRIX_H_
#define MNL_MATRIX_H_

#include "config.h"
#include "vector.h"
#include "range.h"
#include "memtracker.h"

#include <string>
#include <iostream>
#include <memory.h>

/** muu numerics library **
- this is the matrix class, include operator overloads
**/

namespace mnl {
  namespace basics {
    /* forwards  */
    class complexMatrix;

    /* basic matrix class */
    template<class T, class V>
      class basicMatrix {
        public:
          class matrixDotter {
            public:
              T operator()(const basicMatrix<T,V>& A,
                  const basicMatrix<T,V>& B) const
              {
                return A.dot(B);
              }
          };

          basicMatrix(const std::string& n_name, int n_rows, int n_cols, bool clear=true, T** n_data=NULL) :
            m_name(n_name), 
            external_data(false), 
            m_rows(n_rows), 
            m_cols(n_cols)
        {
          if( !n_data ) {
            ALLOCATELOG(m_name,n_rows*n_cols);
            m_data = allocate(n_rows,n_cols,clear);
#ifdef MNL_MEMORY_VERBOSE_
            m_memId = utilities::g_tracker.addChunk(utilities::memTracker::Matrix,length()*sizeof(T));
#endif
          } else {
            external_data = true;
            m_data = n_data;
          }

          setupVectors();
        }

          basicMatrix(const basicMatrix<T,V>& input) :
            external_data(false),
            m_name(input.name()),
            m_rows(input.rows()),
            m_cols(input.cols())
        {
          ALLOCATELOG(name()+" (copy)");
          m_data = allocate(m_rows,m_cols,false);
#ifdef MNL_MEMORY_VERBOSE_
          m_memId = utilities::g_tracker.addChunk(utilities::memTracker::Matrix,m_rows*m_cols*sizeof(Real));
#endif
          setupVectors();
        }

          virtual ~basicMatrix()
          {
            if( !external_data ) {
              DEALLOCATELOG(name());
              deallocate(m_data);
#ifdef MNL_MEMORY_VERBOSE_
              utilities::g_tracker.removeChunk(m_memId);
#endif
            }
            destroyVectors();
          }

          enum dataMode {
            Copy = 1,
            Reference = 2
          };

          inline T** data()
          {
            return m_data;
          }

          T** data(dataMode mode)
          {
            //                T** result;
            //                if( mode == Copy ) {
            //                    std::cout << "wtf why are we here so often?" << std::endl;
            //                    ALLOCATELOG(name()+" (copy of data)");
            //                    result = allocate(rows(),cols(),false);
            //#ifdef MNL_MEMORY_VERBOSE_
            //                    m_memId = utilities::g_tracker.addChunk(utilities::memTracker::Matrix,rows()*cols()*sizeof(T));
            //#endif
            //                    memcpy(m_data[0],m_data[0],rows()*cols()*sizeof(T));
            //                } else 
            //                    result = m_data;

            //                return( result );
          }

          void data(T** n_data)
          {
            if( !external_data ) {
              DEALLOCATELOG(name()+" (assigning external data)");
              deallocate(m_data);
#ifdef MNL_MEMORY_VERBOSE_
              utilities::g_tracker.removeChunk(m_memId);
#endif
            }
            m_data = n_data;
            external_data = true;
            setupVectors();
          }

          inline const T** data() const
          {
            return( const_cast<const T**>(m_data) );
          }

          T max() const
          {
            return *std::max_element(m_data()[0],m_data()[0]+length());
          }

          V row(int i) const
          {
            assert( i < rows() );
            V result("Row from "+name(),cols(),false);
            for( int j=0;j<cols();++j )
              result[j] = m_data[j][i];

            return( result );
          }

          void clear()
          {
            assert( m_data );

            memset(m_data[0],0,rows()*cols()*sizeof(T));
          }

          void fill(const T& in)
          {
            assert( m_data );
            for( int i=0;i<length();++i )
              m_data[0][i] = in;
          }

          void clearRow(int k)
          {
            assert( k < rows() );

            for( int j=0;j<cols();++j )
              m_data[j][k] = 0;
          }

          const std::string& name() const
          {
            return( m_name );
          }

          void setName(const std::string& name)
          {
            m_name = name;
          }

          inline const int& rows() const
          {
            return( m_rows );
          }

          inline const int& cols() const
          {
            return( m_cols );
          }

          inline int length() const
          {
            return rows()*cols();
          }

          void info() const
          {
            std::cout << m_name << ": " << rows() << "x" << cols() << std::endl;
          }

          T sum() const
          {
            T result=0;
            for( int i=0;i<length();++i )
              result += m_data[0][i];

            return( result );
          }

          virtual T dot(const basicMatrix<T,V>& in) const = 0;

          const matrixDotter& getDotter() const
          {
            return m_dotter;
          }

          const V& operator[](const int index) const // returns column vector !
          {
            assert( index < cols() );

            return( *m_vector[index] );
          }

          inline V& operator[](const int index)
          {
            assert( index < cols() );

            return( *m_vector[index] );
          }
        protected:
          std::string m_name;
          std::vector<V*> m_vector;
          bool external_data;
          int m_rows, m_cols;
          T** m_data;
          matrixDotter m_dotter;
#ifdef MNL_MEMORY_VERBOSE_
          int m_memId;
#endif
          static T** allocate(int n_rows, int n_cols, bool clear)
          {
            T** result = new T*[n_cols];
            result[0] = new T[n_rows*n_cols];
            for( int i=1;i<n_cols;i++ )
              result[i] = (T*)(result[0]+i*n_rows); // manually ensure nice linear memory alignment in fortran order!

            if( clear )
              memset(result[0],0,n_rows*n_cols*sizeof(T));

            return( result );
          }

          static void deallocate(T** array)
          {
            delete[] array[0];
            delete[] array;
          }

          void setupVectors()
          {
            destroyVectors();

            std::stringstream s;
            for( int i=0;i<cols();++i ) {
              //                    s << "Column " << i << " from " << name();
              m_vector.push_back(new V(s.str(),m_rows,false,m_data[i]));
            }
          }

          void destroyVectors()
          {
            for( int i=0;i<m_vector.size();++i )
              delete m_vector[i];

            m_vector.clear();
          }
        private:
          void operator=(const basicMatrix<T,V>& input)
          {
          }
      };

    class Matrix : public basicMatrix<Real,Vector> {
      public:
        enum fileFormat {
          ASCII = 1,
          HDF5  = 2
        };

        using basicMatrix<Real,Vector>::operator[];
        using basicMatrix<Real,Vector>::row;

        Matrix(const std::string& n_name, int n_rows, int n_cols, bool clear=true, Real** n_data=NULL);
        Matrix(const Matrix&);
        Matrix(const std::string& n_name, const Matrix& matrix);
        virtual ~Matrix();

        void print(int precision=5) const;
        void print(std::string& str, int precision=5) const;
        void save(const std::string& filename, fileFormat format=ASCII) const;
        void load(const std::string& filename);

        void row(int i, const Vector& n_row);

        static Matrix identity(int M, int N=-1);
        static Matrix diag(const Vector& diag);
        static Matrix ones(int M, int N=-1);
        Matrix transposed() const;
        Matrix submatrix(const utilities::Range& r1, const utilities::Range& r2) const;

        complexVector eigenValues(bool destroy=false) const;
        complexVector eigenValues(complexMatrix &matrix) const;
        Vector eigenValuesSym(bool destroy=false) const;
        void eigenValuesSym(Matrix& evec, Vector& eigs) const;
        void generalizedEigenvaluesSPD(const Matrix& B, Matrix& Q, Vector& eigs) const;
        void LUFactorize();
        void LUSolve(Vector& b);
        void LLSolve(Vector& b);
        void LUSolve(Matrix& B);
        void LLSolve(Matrix& B);
        bool invert();
        Matrix choleskyFactorize();

        void axpy(const Real alpha, const Matrix& x);
        void axpyRow(int k, const Real alpha, const Vector& x);
        void axpyRow(int k, const Real alpha, const Matrix& x);
        void scaleRow(int k, const Real alpha);
        void copyRow(int k, const Vector& x);
        void copyRow(int k, const Matrix& x);

        Real dot(const basicMatrix<Real,Vector>& x) const;

        Matrix& operator =(const Matrix&);
        Matrix& operator =(const Real** const);
        Matrix& operator =(const Real);

        void operator *=(const Matrix& matrix);
        void operator *=(Real factor);
        void operator *=(const Vector& vector); // vector == diag(vector)

        void operator +=(const Matrix& matrix);
        void operator +=(const Real);
        void operator -=(const Matrix& matrix);
        void operator -=(const Real);
        void operator +=(const Vector& vector); // vector == diag(vector)
      protected:
        int* pivots;
    }; 

    class complexMatrix : public basicMatrix<std::complex<Real>,complexVector> {
      public:
        using basicMatrix<std::complex<Real>,complexVector>::row;
        using basicMatrix<std::complex<Real>,complexVector>::operator[];

        complexMatrix(const std::string& n_name, int n_rows, int n_cols, bool clear=true, std::complex<Real>** n_data = NULL);
        complexMatrix(const complexMatrix& matrix);
        complexMatrix(const Matrix& matrix);
        complexMatrix(const std::string& n_name, const complexMatrix& matrix);
        virtual ~complexMatrix();

        void print(int precision=5) const;
        void save(const std::string name, Matrix::fileFormat format=Matrix::ASCII) const;

        Matrix real() const;
        Matrix imag() const;
        complexMatrix transposed() const;
        complexMatrix submatrix(const utilities::Range& r1, const utilities::Range& r2) const;

        void row(int i, const complexVector& n_row);

        void axpy(const std::complex<Real> alpha, const complexMatrix& x);
        void axpy(const Real alpha, const complexMatrix& x);
        void axpy(const Real alpha, const Matrix& x);
        void axpyRow(int k, const Real alpha, const Vector& x);
        void axpyRow(int k, const Real alpha, const complexVector& x);
        void axpyRow(int k, const std::complex<Real>& alpha, const complexVector& x);
        void axpyRow(int k, const Real alpha, const Matrix& x);
        void axpyRow(int k, const Real alpha, const complexMatrix& x);
        void axpyRow(int k, const std::complex<Real>& alpha, const complexMatrix& x);
        void scaleRow(int k, const Real alpha);
        void scaleRow(int k, const std::complex<Real>& alpha);
        void copyRow(int k, const complexMatrix& x);
        void copyRow(int k, const complexVector& x);

        std::complex<Real> dot(const basicMatrix<std::complex<Real>,complexVector>& x) const;

        complexMatrix& operator=(const complexMatrix&);
        complexMatrix& operator=(const Matrix& matrix);

        void operator +=(const Matrix&);
        void operator +=(const complexMatrix&);

        void operator -=(const Matrix&);
        void operator -=(const complexMatrix&);
        void operator -=(const Real);
        void operator -=(const std::complex<Real>);

        void operator *=(std::complex<Real> factor);
        void operator *=(Real factor);
    };

    /* matrix-matrix operators */
    const Matrix operator*(const Matrix&, const Matrix&);
    const Matrix multTranspose(const Matrix&, const Matrix&, char transA, char transB, Real scaleA=1.f);
    const Vector multTranspose(const Matrix&, const Vector&, char transA, Real alpha=1.f);
    void multTranspose(Matrix&, const Matrix&, const Matrix&, char transA, char transB, Real scaleA=1.f, Real scaleB=0.f, int skipcol=0);
    const complexMatrix multTranspose(const complexMatrix&, const complexMatrix&, char transA, char transB, const std::complex<Real>& scaleA=mnlComplexOne);
    void multTranspose(complexMatrix&, const complexMatrix&, const complexMatrix&, char transA, char transB, const std::complex<Real>& scaleA=mnlComplexOne, const std::complex<Real>& scaleB=mnlComplexZero);
    const complexMatrix operator*(const complexMatrix&, const complexMatrix&);
    void multPointwise(complexMatrix& result, const complexMatrix& phi, const complexMatrix& u);
    void multPointwise(Matrix& result, const Matrix& phi, const Matrix& u);
    void multPointwise(Matrix& result, const Matrix& u);
    void multPointwise(Matrix& result, const Matrix& phi, const Matrix& u, Real alpha, Real beta);
    void multPointwiseReal(complexMatrix& result, const complexMatrix& phi, const complexMatrix& u);
    void multPointwiseReal(complexMatrix& phi, const complexMatrix& u);
    const Matrix kron(const Matrix& A, const Matrix& B);

    const Matrix operator -(const Matrix&,const Matrix&);
    const Matrix operator +(const Matrix&,const Matrix&);
    const Matrix operator *(Real, const Matrix&);

    /* matrix-vector products */
    const Vector operator *(const Matrix&, const Vector&);
    const complexVector operator *(const complexMatrix&, const complexVector&);
    const complexVector operator *(const Matrix&, const complexVector&);
    const complexVector operator *(const complexMatrix&, const Vector&);

    /* matrix-scalar products  */
    const complexMatrix operator *(const complexMatrix&, const std::complex<Real>);

    const Matrix loadMatrix(const std::string& filename);
  } // namespace basics
} // namespace mnl
 
#endif

