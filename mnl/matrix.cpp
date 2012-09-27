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

#include "matrix.h"
#ifdef MNL_MEMORY_VERBOSE_
#include "memtracker.h"
#endif

#include <assert.h>
#include <sstream>
#include <fstream>
#include <functional>
#include <algorithm>

using namespace std;

namespace mnl {
  namespace basics {
    Matrix::Matrix(const string& n_name, int n_rows, int n_cols, bool clear, Real** n_data) : 
      basicMatrix<Real,Vector>(n_name,n_rows,n_cols,clear,n_data)
    {
      pivots = NULL;
    }

    Matrix::Matrix(const Matrix& matrix) :
      basicMatrix<Real,Vector>(matrix)
    {
      pivots = NULL;

      *this = matrix;
    }

    Matrix::Matrix(const string& n_name, const Matrix& matrix) : 
      basicMatrix<Real,Vector>(n_name,matrix.rows(),matrix.cols(),false,const_cast<Real**>(matrix.data()))
    {
      pivots = NULL;
    }

    Matrix::~Matrix()
    {
      if( pivots )
        delete[] pivots;
    }

    void Matrix::print(int precision) const
    {
      cout << name() << "(" << rows() << "x" << cols() << ") = " << endl;
      int oldp = cout.precision();
      cout.precision(precision);
      for( int i=0;i<rows();i++ ) {
        cout << "[ ";
        for( int j=0;j<cols();j++ ) {
          if( fabs(data()[j][i]) < 1.e-13f )
            cout << 0 << " ";
          else	
            cout << (Real)data()[j][i] << " ";
        }
        cout << "]" << endl;
      }
      cout.precision(oldp);
    }

    void Matrix::print(string& str, int precision) const
    {
      stringstream s;
      s.precision(precision);
      s << name() << "(" << rows() << "x" << cols() << ") = " << endl;
      for( int i=0;i<rows();i++ ) {
        s << "[ ";
        for( int j=0;j<cols();j++ ) {
          if( fabs(data()[j][i]) < 1.e-14f )
            s << 0 << " ";
          else	
            s << (Real)data()[j][i] << " ";
        }
        s << "]" << endl;
      }
      str = s.str();
    }

    void Matrix::save(const string& filename, fileFormat format) const
    {
      if( format == ASCII ) {
        ofstream f;
        if( filename.find(".asc") != string::npos )
          f.open(filename.c_str());
        else
          f.open((filename+".asc").c_str());
        f.precision(16);
        f << fixed;
        f << "% " << rows() << '\t' << cols() << '\t' << name() << endl;
        for( size_t i=0;i<rows();++i ) {
          for( size_t j=0;j<cols();++j )
            f << (m_data[j])[i] << '\t';
          f << endl;
        }
      }
    }

    void Matrix::load(const string& filename)
    {
      int n_rows;
      int n_cols;

      stringstream s;
      string line;
      ifstream f;
      f.open(filename.c_str());
      getline(f, line); // read dimensions
      line.erase(0,2);
      int pos = line.rfind('\t');
      if( pos != string::npos )
        line.erase(pos,line.size()-pos);
      s << line;
      s >> n_rows;
      s >> n_cols;

      cout << "rows " << rows() << " cols " << cols() << endl;
      assert( n_rows == rows() && n_cols == cols() );

      for (size_t i=0;i<n_rows;++i) {  // read one line
        getline(f, line); // read actual data
        s.clear();
        s << line;
        Real t;
        for (size_t j=0;j<n_cols;++j) {
          s >> t;
          m_data[j][i]= t;
        }
      }
    }

    void Matrix::row(int i, const Vector& n_row)
    {
      assert( n_row.length() == cols() );

      BLASRPFX(copy,BLASINT(m_cols),BLASPREAL(n_row.data()),BLASINT(mnlIntOne),BLASPREAL(m_data[0]+i),BLASINT(m_rows));
    }

    Matrix Matrix::diag(const Vector& diag)
    {
      string name2("Diagonal matrix");
      Matrix result(name2,diag.length(),diag.length());

      int inc=result.rows()+1;
      BLASRPFX(copy,BLASINT(diag.length()),BLASPREAL(diag.data()),BLASINT(mnlIntOne),BLASPREAL(result.data()[0]),BLASINT(inc));

      return( result );
    }

    Matrix Matrix::ones(int M, int N)
    {
      int len;
      if( N == -1 )
        N=M;
      if( M < N )
        len = M;
      else
        len = N;
      Matrix result("One matrix",M,N);
      for( int i=0;i<result.length();++i )
        result[0][i] = 1;

      return( result );
    }

    Matrix Matrix::identity(int M, int N)
    {
      int len;
      if( N == -1 )
        N=M;
      if( M < N )
        len = M;
      else
        len = N;
      Matrix result("Identity matrix",M,N);
      for( int i=0;i<len;++i )
        result[i][i] = 1;

      return( result );
    }

    Matrix Matrix::transposed() const
    {
      Matrix result("transposed of "+name(),cols(),rows());
      for( int j=0;j<rows();++j )
        result[j] = row(j);

      return result;
    }

    Matrix Matrix::submatrix(const utilities::Range& r1, const utilities::Range& r2) const
    {
      Matrix result("Submatrix of "+name(),r1.size(),r2.size(),false);

      for( int j=0;j<r2.size();++j ) 
        for( int i=0;i<r1.size();++i ) {
          assert( (r1[i] < rows()) && (r2[j] < cols()) );
          result[j][i] = data()[r2[j]][r1[i]];
        }

      return( result );
    }

    complexVector Matrix::eigenValues(bool destroy) const
    {
      Real **d2;
      if( !destroy ) {
        d2 = allocate(m_rows,m_cols,false);
        int N=m_rows*m_cols;
        BLASRPFX(copy,BLASINT(N),m_data[0],BLASINT(mnlIntOne),d2[0],BLASINT(mnlIntOne));
      } else
        d2 = m_data;
      Real *wr = new Real[m_rows];
      Real *wi= new Real[m_rows];
      int ilo=1;
      int ihi=m_rows;
      Real *work = new Real[m_rows]; // INVESTIGATE BLOCK PROPERTIES!
      Real *tau = new Real[m_rows-1];

      int info;

      /* reduce to upper Hessenberg form */
      LAPACKRPFX(gehrd,&ihi,&ilo,&ihi,d2[0],&ihi,tau,work,&ihi,&info);
      assert( info == 0 );

      /* find eigenvalues */
      char job = 'E';
      char compz = 'N';
      LAPACKRPFX(hseqr,&job,&compz,&ihi,&ilo,&ihi,d2[0],&ihi,wr,wi,NULL,&ihi,work,&ihi,&info);
      assert( info == 0 );

      complexVector result("Eigenvalues of "+name(),m_cols,wr,wi);

      if( !destroy )
        deallocate(d2);
      delete[] tau;
      delete[] work;
      delete[] wr;
      delete[] wi;

      return( result );
    }

    complexVector Matrix::eigenValues(complexMatrix &matrix) const
    {
      LOG("Finding eigenvalues and eigenvectors of "+name());
      Real **d2 = allocate(m_rows,m_cols,false);
      Real *wr = new Real[m_rows];
      Real *wi= new Real[m_rows];
      int N=m_rows*m_cols;
      BLASRPFX(copy,BLASINT(N),m_data[0],BLASINT(mnlIntOne),d2[0],BLASINT(mnlIntOne));
      int ilo=1;
      int ihi=m_rows;
      Real *work = new Real[m_rows]; // INVESTIGATE BLOCK PROPERTIES!
      Real *tau = new Real[m_rows-1];

      int info;

      /* reduce to upper Hessenberg form */
      LAPACKRPFX(gehrd,&ihi,&ilo,&ihi,d2[0],&ihi,tau,work,&ihi,&info);
      assert( info == 0 );

      /* find eigenvalues */
      char job = 'E';
      char compz = 'N';

      LAPACKRPFX(hseqr,&job,&compz,&ihi,&ilo,&ihi,d2[0],&ihi,wr,wi,NULL,&ihi,work,&ihi,&info);
      assert( info == 0 );

      complexVector result("Eigenvalues of "+name(),m_cols,wr,wi);

      deallocate(d2);
      delete[] tau;
      delete[] work;
      delete[] wr;
      delete[] wi;

      return( result );
    }

    Vector Matrix::eigenValuesSym(bool destroy) const
    {
      LOG("Finding eigenvalues of "+name()+" (symmetric)");
      assert( m_rows == m_cols ); // required for symmetry
      Real** d2;
      if( !destroy ) {
        ALLOCATELOG("buffers for eigenValuesSym");
        d2 = allocate(m_rows,m_cols,false);
        int N=m_rows*m_cols;
        BLASRPFX(copy,BLASINT(N),m_data[0],BLASINT(mnlIntOne),d2[0],BLASINT(mnlIntOne));
      } else
        d2 = m_data;
      int N=m_rows;

      Real* D = new Real[N];
      Real* E = new Real[N-1];
      Real* tau = new Real[N-1];
      Real* work = new Real[2*N];

      int info;
      char compz = 'N';
      char uplo = 'U';

      /* reduce to tridiagonal form */
      LAPACKRPFX(sytrd,&uplo,&N,d2[0],&N,D,E,tau,work,&N,&info);
      assert( info == 0 );

      /* find eigenvalues */
      LAPACKRPFX(steqr,&compz,&N,D,E,NULL,&N,work,&info);
      assert( info == 0 );

      Vector result("Eigenvalues of "+name(),m_cols);
      result = D;

      if( !destroy ) {
        DEALLOCATELOG("buffers for eigenValuesSym");
        deallocate(d2);
      }
      delete[] tau;
      delete[] work;
      delete[] D;
      delete[] E;

      return( result );
    }

    void Matrix::eigenValuesSym(Matrix& evec, Vector& eigs) const
    {
      LOG("Finding eigenvalues and eigenvectors of "+name()+" (symmetric)");
      assert( m_rows == m_cols ); // required for symmetry
      evec = *this;
      Real** d2 = evec.data();
      int N=m_rows;

      ALLOCATELOG("buffers for eigenValuesSym");
      Real* D = new Real[N];
      Real* E = new Real[N-1];
      Real* tau = new Real[N-1];
      Real* work = new Real[2*N-2];
      int info;
      char uplo = 'U';

      /* reduce to tridiagonal form */
      int N2 = 2*N-2;
      LAPACKRPFX(sytrd,&uplo,&N,d2[0],&N,D,E,tau,work,&N2,&info);
      assert( info == 0 );

      /* reconstruct orthogonal matrix */
      LAPACKRPFX(orgtr,&uplo,&N,d2[0],&N,tau,work,&N2,&info);
      assert( info == 0 );

      /* find eigenvalues and eigenvectors */
      char compz = 'V';
      LAPACKRPFX(steqr,&compz,&N,D,E,d2[0],&N,work,&info);
      assert( info == 0 );

      eigs = D;

      DEALLOCATELOG("buffers for eigenValuesSym");
      delete[] tau;
      delete[] work;
      delete[] D;
      delete[] E;
    }

    void Matrix::generalizedEigenvaluesSPD(const Matrix& B, Matrix& Q, Vector& eigs) const
    {
      LOG("Finding generalized eigenvalues and eigenvectors of "+name()+" (SPD)");
      Matrix B2(B);
      Matrix C(*this);
      char Uplo = 'U';
      int N = B.rows();
      int info;

      /* We first find the Cholesky factorization of B = U*U^T */
#ifdef _AIX
#undef dpotrf
#endif
      LAPACKRPFX(potrf,&Uplo,&N,B2.data()[0],&N,&info);
      assert( info == 0 );

      /* Moving on, let Lapack construct the C matrix we then pass to eigenvalues. CQ=LQ */
      int itype = 1;
      Real** cData = C.data();
      LAPACKRPFX(sygst,&itype,&Uplo,&N,cData[0],&N,B2.data()[0],&N,&info);
      assert( info == 0 );

      /* Find eigenvalues and eigenvectors of C */
      C.eigenValuesSym(Q,eigs);
      eigs.setName("Generalized eigenvalues of "+C.name());

      /* Obtain correct eigenvectors Q = U^-1*Q */
      Real** qData = Q.data();
      for( int i=0;i<Q.cols();i++ )
        BLASRPFX(trsv,BLASCHAR(Uplo),BLASCHAR(mnlNoTrans),BLASCHAR(mnlNoTrans),BLASINT(Q.rows()),B2.data()[0],BLASINT(Q.rows()),qData[i],BLASINT(mnlIntOne));
    }

    void Matrix::LUFactorize()
    {
      if( !pivots ) {
        ALLOCATELOG("allocating pivot array");
        int dim = rows();
        if( cols() < dim )
          dim = cols();
        pivots = new int[dim];
      }
      int info;
#ifdef _AIX
#undef dgetrf
#endif
      LAPACKRPFX(getrf,&m_rows,&m_cols,data()[0],&m_rows,pivots,&info);

      assert( info == 0 );
    }

    Matrix Matrix::choleskyFactorize()
    {
      assert( rows() == cols() );

      LOG("Finding Cholesky factorization of "+name());

      Matrix result(*this);
      for( int i=0;i<result.rows();++i )
        for( int j=0; j < i;++j )
          result[i][j] = 0;

      result.setName("Cholesky factorization of "+name());
      char Uplo = 'L';
      int N = rows();
      int info;

      /* find the Cholesky factorization of B = U*U^T */
      LAPACKRPFX(potrf,&Uplo,&N,result.data()[0],&N,&info);
      assert( info == 0 );

      return result;
    }

    void Matrix::LUSolve(Vector& b)
    {
      if( !pivots )
        LUFactorize();

      char trans='N';
      int info;
      int N=1;
      int M=b.length();
#ifdef _AIX
#undef dgetrs
#endif
      LAPACKRPFX(getrs,&trans,&m_cols,&N,data()[0],&m_cols,pivots,b.data(),&M,&info);

      assert( info == 0 );
    }

    void Matrix::LUSolve(Matrix& B) // use with caution!
    {
      if( !pivots )
        LUFactorize();

      char trans='N';
      int info;
      int N=1;
      int M=B.length();
#ifdef _AIX
#undef dgetrs
#endif
      LAPACKRPFX(getrs,&trans,&m_cols,&N,data()[0],&m_cols,pivots,B.data()[0],&M,&info);

      assert( info == 0 );
    }

    void Matrix::LLSolve(Vector& b)
    {
      char trans='N';
      int info;
      int N=1;
      int M=b.length();
#ifdef _AIX
#undef dpotrs
#endif
      char Uplo='L';
      LAPACKRPFX(potrs,&Uplo,&m_cols,&N,data()[0],&m_cols,b.data(),&M,&info);

      assert( info == 0 );
    }


    void Matrix::LLSolve(Matrix& B) // use with caution!
    {
      char trans='N';
      int info;
      int N=1;
      int M=B.length();
#ifdef _AIX
#undef dpotrs
#endif
      char Uplo='L';
      LAPACKRPFX(potrs,&Uplo,&m_cols,&N,data()[0],&m_cols,B.data()[0],&M,&info);

      assert( info == 0 );
    }

    bool Matrix::invert()
    {
      if( !pivots )
        LUFactorize();

      int info;
#ifdef _AIX
#undef dgetri
#endif
      Real* work = new Real[2*m_cols];
      int lwork = 2*m_cols;
      LAPACKRPFX(getri,&m_cols,data()[0],&m_rows,pivots,work,&lwork,&info);

      assert( info == 0 );
      delete[] work;

      return true;
    }

    void Matrix::axpy(const Real alpha, const Matrix& x)
    {
      assert( (x.cols() == cols() ) && (x.rows() == rows()) );

      int N=rows()*cols();
      BLASRPFX(axpy,BLASINT(N),BLASREAL(alpha),const_cast<Real*>(x.data()[0]),BLASINT(mnlIntOne),data()[0],BLASINT(mnlIntOne));
    }

    void Matrix::axpyRow(int k, const Real alpha, const Vector& x)
    {
      assert( k < rows() && x.length() == cols() );

      BLASRPFX(axpy,BLASINT(cols()),BLASREAL(alpha),BLASPREAL(x.data()),BLASINT(mnlIntOne),BLASPREAL(data()[0]+k),BLASINT(rows()));
    }

    void Matrix::axpyRow(int k, const Real alpha, const Matrix& x)
    {
      assert( k < rows() && x.cols() == cols() && x.rows() == rows() );

      BLASRPFX(axpy,BLASINT(cols()),BLASREAL(alpha),BLASPREAL(x.data()[0]+k),BLASINT(rows()),BLASPREAL(data()[0]+k),BLASINT(rows()));
    }

    void Matrix::scaleRow(int k, const Real alpha)
    {
      assert( k < rows() );

      BLASRPFX(scal,BLASINT(cols()),BLASREAL(alpha),BLASPREAL(data()[0]+k),BLASINT(rows()));
    }

    void Matrix::copyRow(int k, const Matrix& x)
    {
      assert( k < rows() && x.cols() == cols() && rows() == x.rows() );

      BLASRPFX(copy,BLASINT(cols()),BLASPREAL(x.data()[0]+k),BLASINT(rows()),BLASPREAL(data()[0]+k),BLASINT(rows()));
    }

    void Matrix::copyRow(int k, const Vector& x)
    {
      assert( k < rows() && x.length() == cols() );

      BLASRPFX(copy,BLASINT(cols()),BLASPREAL(x.data()),BLASINT(mnlIntOne),BLASPREAL(data()[0]+k),BLASINT(rows()));
    }

    Real Matrix::dot(const basicMatrix<Real,Vector>& x) const
    {
      assert( cols() == x.cols() && rows() == x.rows() );

      int N=cols()*rows();
      return BLASRPFX(dot,BLASINT(N),BLASPREAL(x.data()[0]),BLASINT(mnlIntOne),BLASPREAL(data()[0]),BLASINT(mnlIntOne)); 
    }

    Matrix& Matrix::operator =(const Matrix& matrix)
    {
      assert( (cols() == matrix.cols()) && (rows() == matrix.rows()) );

      int N=m_rows*m_cols;
      BLASRPFX(copy,BLASINT(N),BLASPREAL(matrix.data()[0]),BLASINT(mnlIntOne),m_data[0],BLASINT(mnlIntOne));

      return( *this );
    }

    Matrix& Matrix::operator=(const Real** const data)
    {
      int N=m_rows*m_cols;
      BLASRPFX(copy,BLASINT(N),const_cast<Real*>(data[0]),BLASINT(mnlIntOne),m_data[0],BLASINT(mnlIntOne));

      return( *this );
    }

    Matrix& Matrix::operator =(const Real alpha)
    {
      for( int i=0;i<length();++i )
        m_data[0][i] = alpha;

      return( *this );
    }

    void Matrix::operator *=(const Matrix& matrix)
    {
      assert( cols() == matrix.rows() );

      Real** m_data2 = allocate(cols(),matrix.rows(),false);

      BLASRPFX(gemm,BLASCHAR(mnlNoTrans),BLASCHAR(mnlNoTrans),BLASINT(m_rows),BLASINT(m_cols),BLASINT(matrix.cols()),BLASREAL(mnlRealOne),data()[0],BLASINT(m_rows),BLASPREAL(matrix.data()[0]),BLASINT(matrix.rows()),BLASREAL(mnlRealZero),BLASPREAL(m_data2[0]),BLASINT(m_rows));
      DEALLOCATELOG(name() + " (MxM)");
      deallocate(m_data);
#ifdef MNL_MEMORY_VERBOSE_
      utilities::g_tracker.removeChunk(m_memId);
      m_memId = utilities::g_tracker.addChunk(utilities::memTracker::Matrix,cols()*matrix.rows()*sizeof(Real));
#endif
      m_data = m_data2;
    }

    void Matrix::operator *=(const Vector& vector) // vector == diag(vector)
    {
#if defined(DEBUG)
      int len;
      if( rows() < cols() )
        len = rows();
      else
        len = cols();
      assert( vector.length() == len );
#endif

      for( int i=0;i<vector.length();i++ )
        for( int j=0;j<cols();j++ )
          m_data[i][j] *= vector[i];
    }

    void Matrix::operator *=(Real factor)
    {
      int N=rows()*cols();
      BLASRPFX(scal,BLASINT(N),BLASREAL(factor),BLASPREAL(m_data[0]),BLASINT(mnlIntOne));
    }

    void Matrix::operator +=(const Matrix& matrix)
    {
      assert( (rows() == matrix.rows()) && (cols() == matrix.cols()) );

      int N=rows()*cols();
      BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealOne),BLASPREAL(matrix.data()[0]),BLASINT(mnlIntOne),BLASPREAL(data()[0]),BLASINT(mnlIntOne));
    }

    void Matrix::operator +=(const Vector& vector)
    {
#if defined(DEBUG)
      int len;
      if( rows() < cols() )
        len = rows();
      else
        len = cols();
      assert( vector.length() == len );
#endif

      int N=rows()+1;
      BLASRPFX(axpy,BLASINT(vector.length()),BLASREAL(mnlRealOne),const_cast<Real*>(vector.data()),BLASINT(mnlIntOne),m_data[0],BLASINT(N));
    }

    void Matrix::operator +=(const Real alpha)
    {
      for( int i=0;i<length();++i )
        m_data[0][i] += alpha;
    }

    void Matrix::operator -=(const Matrix& matrix)
    {
      assert( (rows() == matrix.rows()) && (cols() == matrix.cols()) );

      int N=rows()*cols();
      BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealMinusOne),BLASPREAL(matrix.data()[0]),BLASINT(mnlIntOne),BLASPREAL(data()[0]),BLASINT(mnlIntOne));
    }

    void Matrix::operator -=(const Real alpha)
    {
      for( int i=0;i<length();++i )
        m_data[0][i] -= alpha;
    }

    complexMatrix::complexMatrix(const string& n_name, int n_rows, int n_cols, bool clear, std::complex<Real>** n_data) : 
      basicMatrix<std::complex<Real>,complexVector>(n_name,n_rows,n_cols,clear,n_data)
    {
    }

    complexMatrix::complexMatrix(const complexMatrix& matrix) : 
      basicMatrix<std::complex<Real>,complexVector>(matrix.name()+" (copy)",matrix.rows(),matrix.cols(),false)
    {
      *this = matrix;
    }

    complexMatrix::complexMatrix(const Matrix& matrix) : 
      basicMatrix<std::complex<Real>,complexVector>(matrix.name(),matrix.rows(),matrix.cols())
    {
      *this = matrix;
    }

    complexMatrix::complexMatrix(const string& name, const complexMatrix& matrix) : 
      basicMatrix<std::complex<Real>,complexVector>(name,matrix.rows(),matrix.cols(),false,const_cast<std::complex<Real>**>(matrix.data()))
    {
    }	

    complexMatrix::~complexMatrix()
    {
    }

    void complexMatrix::print(int precision) const
    {
      cout << name() << "(" << rows() << "x" << cols() << ") = " << endl;
      int oldp = cout.precision();
      cout.precision(precision);
      for( int i=0;i<rows();i++ ) {
        cout << "[ ";
        for( int j=0;j<cols();j++ ) {
          Real re,im;
          if( fabs(data()[j][i].real()) < 1.e-18f )
            re = 0.f;
          else
            re = data()[j][i].real();
          if( fabs(data()[j][i].imag()) < 1.e-18f )
            im = 0.f;
          else
            im = data()[j][i].imag();
          if( im >= 0.f )
            cout << re << "+" << im << "i ";
          else 
            cout << re << "" << im << "i ";
        }
        cout << "]" << endl;
      }
      cout.precision(oldp);
    }

    void complexMatrix::save(const string name, Matrix::fileFormat format) const
    {
      if( format == Matrix::ASCII ) {
        FILE* f = fopen(name.c_str(),"rw");
        for( int i=0;i<m_rows;i++ ) {
          fprintf(f,"\t");
          for( int j=0;j<m_cols-1;j++ )
            fprintf(f,"%f\t",m_data[j][i].real());
          fprintf(f,"%f\n",m_data[m_cols-1][i].real());
        }
        fclose(f);
      }
    }

    void complexMatrix::row(int i, const complexVector& n_row)
    {
      assert( n_row.length() == cols() );

      BLASCPFX(copy,BLASINT(m_cols),BLASPCOMPLEX(n_row.data()),BLASINT(mnlIntOne),BLASPCOMPLEX(m_data[0]+i),BLASINT(m_rows));
    }

    Matrix complexMatrix::real() const
    {
      Matrix result("Real part of "+name(),rows(),cols(),false);
      Real** data = result.data();

      int N=cols()*rows();
      BLASRPFX(copy,BLASINT(N),BLASPREAL(m_data[0]),BLASINT(mnlIntTwo),data[0],BLASINT(mnlIntOne));

      return( result );
    }

    Matrix complexMatrix::imag() const
    {
      Matrix result("Real part of "+name(),rows(),cols(),false);
      Real** data = result.data();

      int N=cols()*rows();
      BLASRPFX(copy,BLASINT(N),BLASPREAL(m_data[0])+1,BLASINT(mnlIntTwo),data[0],BLASINT(mnlIntOne));

      return( result );
    }

    complexMatrix complexMatrix::transposed() const
    {
      complexMatrix result("transposed of "+name(),cols(),rows());
      for( int j=0;j<rows();++j )
        result[j] = row(j);

      return result;
    }

    complexMatrix complexMatrix::submatrix(const utilities::Range& r1, const utilities::Range& r2) const
    {
      complexMatrix result("Submatrix of "+name(),r1.size(),r2.size(),false);

      for( int j=0;j<r2.size();++j ) 
        for( int i=0;i<r1.size();++i ) {
          assert( (r2[j] < cols()) && (r1[i] < rows()) );
          result[j][i] = data()[r2[j]][r1[i]];
        }

      return( result );
    }

    void complexMatrix::axpy(const std::complex<Real> alpha, const complexMatrix& x)
    {
      assert( (rows() == x.rows()) && (cols() == x.cols()) );

      int N=rows()*cols();
      BLASCPFX(axpy,BLASINT(N),BLASCOMPLEX(alpha),BLASPCOMPLEX(x.data()[0]),BLASINT(mnlIntOne),BLASPCOMPLEX(data()[0]),BLASINT(mnlIntOne));
    }

    void complexMatrix::axpy(const Real alpha, const complexMatrix& x)
    {
      assert( (rows() == x.rows()) && (cols() == x.cols()) );

      int N=rows()*cols()*2;
      BLASRPFX(axpy,BLASINT(N),BLASREAL(alpha),BLASPREAL(x.data()[0]),BLASINT(mnlIntOne),BLASPREAL(data()[0]),BLASINT(mnlIntOne));
    }

    void complexMatrix::axpy(const Real alpha, const Matrix& x)
    {
      assert( (rows() == x.rows()) && (cols() == x.cols()) );

      int N=rows()*cols();
      BLASRPFX(axpy,BLASINT(N),BLASREAL(alpha),const_cast<Real*>(x.data()[0]),BLASINT(mnlIntOne),BLASPREAL(data()[0]),BLASINT(mnlIntTwo));
    }

    void complexMatrix::axpyRow(int k, const Real alpha, const Vector& x)
    {
      assert( k < rows() && x.length() == cols() );

      int N=2*rows();
      BLASRPFX(axpy,BLASINT(cols()),BLASREAL(alpha),BLASPREAL(x.data()),BLASINT(mnlIntOne),BLASPREAL(data()[0]+k),BLASINT(N));
    }

    void complexMatrix::axpyRow(int k, const Real alpha, const Matrix& x)
    {
      assert( k < rows() && x.cols() == cols() && rows() == x.rows() );

      int N=2*rows();
      BLASRPFX(axpy,BLASINT(cols()),BLASREAL(alpha),BLASPREAL(x.data()[0]+k),BLASINT(x.rows()),BLASPREAL(data()[0]+k),BLASINT(N));
    }

    void complexMatrix::axpyRow(int k, const Real alpha, const complexVector& x)
    {
      assert( k < rows() && x.length() == cols() );

      int N=2*rows();
      BLASRPFX(axpy,BLASINT(cols()),BLASREAL(alpha),BLASPREAL(x.data()),BLASINT(mnlIntTwo),BLASPREAL(data()[0]+k),BLASINT(N));
      BLASRPFX(axpy,BLASINT(cols()),BLASREAL(alpha),BLASPREAL(x.data())+1,BLASINT(mnlIntTwo),BLASPREAL(data()[0]+k)+1,BLASINT(N));
    }

    void complexMatrix::axpyRow(int k, const Real alpha, const complexMatrix& x)
    {
      assert( k < rows() && x.cols() == cols() && rows() == x.rows() );

      int N=2*rows();
      BLASRPFX(axpy,BLASINT(cols()),BLASREAL(alpha),BLASPREAL(x.data()[0]+k),BLASINT(N),BLASPREAL(data()[0]+k),BLASINT(N));
      BLASRPFX(axpy,BLASINT(cols()),BLASREAL(alpha),BLASPREAL(x.data()[0]+k)+1,BLASINT(N),BLASPREAL(data()[0]+k)+1,BLASINT(N));
    }

    void complexMatrix::axpyRow(int k, const std::complex<Real>& alpha, const complexVector& x)
    {
      assert( k < rows() && x.length() == cols() );

      BLASCPFX(axpy,BLASINT(cols()),BLASCOMPLEX(alpha),BLASPCOMPLEX(x.data()),BLASINT(mnlIntOne),BLASPCOMPLEX(data()[0]+k),BLASINT(rows()));
    }

    void complexMatrix::axpyRow(int k, const std::complex<Real>& alpha, const complexMatrix& x)
    {
      assert( k < rows() && x.cols() == cols() && rows() == x.rows() );

      BLASCPFX(axpy,BLASINT(cols()),BLASCOMPLEX(alpha),BLASPCOMPLEX(x.data()[0]+k),BLASINT(rows()),BLASPCOMPLEX(data()[0]+k),BLASINT(rows()));
    }

    void complexMatrix::scaleRow(int k, const Real alpha)
    {
      assert( k < rows() );

      int N=2*rows();
      BLASRPFX(scal,BLASINT(cols()),BLASREAL(alpha),BLASPREAL(data()[0]+k),BLASINT(N));
      BLASRPFX(scal,BLASINT(cols()),BLASREAL(alpha),BLASPREAL(data()[0]+k)+1,BLASINT(N));
    }

    void complexMatrix::scaleRow(int k, const std::complex<Real>& alpha)
    {
      assert( k < rows() );

      BLASCPFX(scal,BLASINT(cols()),BLASCOMPLEX(alpha),BLASPCOMPLEX(data()[0]+k),BLASINT(rows()));
    }

    void complexMatrix::copyRow(int k, const complexMatrix& x)
    {
      assert( k < rows() && x.cols() == cols() && rows() == x.rows() );

      BLASCPFX(copy,BLASINT(cols()),BLASPCOMPLEX(x.data()[0]+k),BLASINT(rows()),BLASPCOMPLEX(data()[0]+k),BLASINT(rows()));
    }

    void complexMatrix::copyRow(int k, const complexVector& x)
    {
      assert( k < rows() && x.length() == cols() );

      BLASCPFX(copy,BLASINT(cols()),BLASPCOMPLEX(x.data()),BLASINT(mnlIntOne),BLASPCOMPLEX(data()[0]+k),BLASINT(rows()));
    }

    std::complex<Real> complexMatrix::dot(const basicMatrix<std::complex<Real>,complexVector>& x) const
    {
      assert( cols() == x.cols() && rows() == x.rows() );

      int N=cols()*rows();

      std::complex<Real> result;
      //#ifndef __freebsd__
      //            ZDOTC(result,BLASINT(N),BLASPCOMPLEX(x.data()[0]),BLASINT(mnlIntOne),BLASPCOMPLEX(data()[0]),BLASINT(mnlIntOne));
      //#else
      for( int j=0;j<N;++j )
        result += conj(x.data()[0][j])*data()[0][j];
      //#endif
      return result;
    }

    complexMatrix& complexMatrix::operator=(const complexMatrix& matrix)
    {
      assert( (cols() == matrix.cols()) && (rows() == matrix.rows()) );

      int N=m_rows*m_cols;
      BLASCPFX(copy,BLASINT(N),BLASPCOMPLEX(matrix.data()[0]),BLASINT(mnlIntOne),BLASPCOMPLEX(m_data[0]),BLASINT(mnlIntOne));

      return( *this );
    }

    complexMatrix& complexMatrix::operator=(const Matrix& matrix)
    {
      assert( (cols() == matrix.cols()) && (rows() == matrix.rows()) );

      clear();
      int N=m_rows*m_cols;
      BLASRPFX(copy,BLASINT(N),BLASPREAL(matrix.data()[0]),BLASINT(mnlIntOne),BLASPREAL(m_data[0]),BLASINT(mnlIntTwo));

      return( *this );
    }

    void complexMatrix::operator -=(const Matrix& matrix)
    {
      assert( (rows() == matrix.rows()) && (cols() == matrix.cols()) );

      int N=rows()*cols();
      BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealMinusOne),const_cast<Real*>(matrix.data()[0]),BLASINT(mnlIntOne),BLASPREAL(data()[0]),BLASINT(mnlIntTwo));
    }

    void complexMatrix::operator -=(const complexMatrix& matrix)
    {
      assert( (rows() == matrix.rows()) && (cols() == matrix.cols()) );

      int N=rows()*cols();
      BLASCPFX(axpy,BLASINT(N),BLASCOMPLEX(mnlComplexMinusOne),BLASPCOMPLEX(matrix.data()[0]),BLASINT(mnlIntOne),BLASPCOMPLEX(data()[0]),BLASINT(mnlIntOne));
    }

    void complexMatrix::operator -=(const Real alpha)
    {
      for( int i=0;i<length();++i )
        m_data[0][i] -= alpha;
    }

    void complexMatrix::operator -=(const std::complex<Real> alpha)
    {
      for( int i=0;i<length();++i )
        m_data[0][i] -= alpha;
    }

    void complexMatrix::operator +=(const Matrix& matrix)
    {
      assert( (rows() == matrix.rows()) && (cols() == matrix.cols()) );

      int N=rows()*cols();
      BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealOne),BLASPREAL(matrix.data()[0]),BLASINT(mnlIntOne),BLASPREAL(data()[0]),BLASINT(mnlIntTwo));
    }

    void complexMatrix::operator +=(const complexMatrix& matrix)
    {
      assert( (rows() == matrix.rows()) && (cols() == matrix.cols()) );

      int N=rows()*cols()*2;
      BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealOne),BLASPREAL(matrix.data()[0]),BLASINT(mnlIntOne),BLASPREAL(data()[0]),BLASINT(mnlIntOne));
    }

    void complexMatrix::operator *=(std::complex<Real> factor)
    {
      int N=rows()*cols();
      BLASCPFX(scal,BLASINT(N),BLASCOMPLEX(factor),BLASPCOMPLEX(m_data[0]),BLASINT(mnlIntOne));
    }

    void complexMatrix::operator *=(Real factor)
    {
      int N=rows()*cols()*2;
      BLASRPFX(scal,BLASINT(N),BLASREAL(factor),BLASPREAL(m_data[0]),BLASINT(mnlIntOne));
    }

    const Vector operator*(const Matrix& lhs, const Vector& rhs)
    {
      assert( lhs.cols() == rhs.length() );

      Vector result(lhs.name()+"*"+rhs.name(),lhs.rows(),false);
      BLASRPFX(gemv,BLASCHAR(mnlNoTrans),BLASINT(lhs.rows()),BLASINT(lhs.cols()),BLASREAL(mnlRealOne),const_cast<Real*>(lhs.data()[0]),BLASINT(lhs.rows()),const_cast<Real*>(rhs.data()),BLASINT(mnlIntOne),BLASREAL(mnlRealZero),result.data(),BLASINT(mnlIntOne));

      return( result );
    }

    const complexVector operator*(const Matrix& lhs, const complexVector& rhs)
    {
      assert( (lhs.cols() == rhs.length()) );

      complexVector result(lhs.name()+"*"+rhs.name(),rhs.length(),false);
      BLASRPFX(gemv,BLASCHAR(mnlNoTrans),BLASINT(lhs.rows()),BLASINT(lhs.cols()),BLASREAL(mnlRealOne),BLASPREAL(lhs.data()[0]),BLASINT(lhs.rows()),BLASPREAL(rhs.data()),BLASINT(mnlIntTwo),BLASREAL(mnlRealZero),BLASPREAL(result.data()),BLASINT(mnlIntTwo));
      BLASRPFX(gemv,BLASCHAR(mnlNoTrans),BLASINT(lhs.rows()),BLASINT(lhs.cols()),BLASREAL(mnlRealOne),BLASPREAL(lhs.data()[0]),BLASINT(lhs.rows()),BLASPREAL(rhs.data())+1,BLASINT(mnlIntTwo),BLASREAL(mnlRealZero),BLASPREAL(result.data())+1,BLASINT(mnlIntTwo));

      return( result );
    }

    const complexVector operator*(const complexMatrix& lhs, const complexVector& rhs)
    {
      assert( lhs.cols() == rhs.length() );

      complexVector result(lhs.name()+"*"+rhs.name(),lhs.rows(),false);
      BLASCPFX(gemv,BLASCHAR(mnlNoTrans),BLASINT(lhs.rows()),BLASINT(lhs.cols()),BLASCOMPLEX(mnlComplexOne),BLASPCOMPLEX(lhs.data()[0]),BLASINT(lhs.rows()),BLASPCOMPLEX(rhs.data()),BLASINT(mnlIntOne),BLASCOMPLEX(mnlComplexZero),BLASPCOMPLEX(result.data()),BLASINT(mnlIntOne));

      return( result );
    }

    const complexVector operator*(const complexMatrix& lhs, const Vector& rhs)
    {
      assert( lhs.cols() == rhs.length() );

      complexVector fool(rhs);
      complexVector result(lhs.name()+"*"+rhs.name(),lhs.rows(),false);
      BLASCPFX(gemv,BLASCHAR(mnlNoTrans),BLASINT(lhs.rows()),BLASINT(lhs.cols()),BLASCOMPLEX(mnlComplexOne),BLASPCOMPLEX(lhs.data()[0]),BLASINT(lhs.rows()),BLASPCOMPLEX(fool.data()),BLASINT(mnlIntOne),BLASCOMPLEX(mnlComplexZero),BLASPCOMPLEX(result.data()),BLASINT(mnlIntOne));

      return( result );
    }

    const Matrix operator*(const Matrix& one, const Matrix& two)
    {
      assert( one.cols() == two.rows() );

      Matrix result(one.name()+"*"+two.name(),one.rows(),two.cols(),false);
      BLASRPFX(gemm,BLASCHAR(mnlNoTrans),BLASCHAR(mnlNoTrans),BLASINT(one.rows()),BLASINT(two.cols()),BLASINT(two.rows()),BLASREAL(mnlRealOne),BLASPREAL(one.data()[0]),BLASINT(one.rows()),BLASPREAL(two.data()[0]),BLASINT(two.rows()),BLASREAL(mnlRealZero),BLASPREAL(result.data()[0]),BLASINT(one.rows()));

      return( result );
    }

    const Matrix multTranspose(const Matrix& one, const Matrix& two, char trans1, char trans2, Real scaleA)
    {
      int M,N,K,lda,ldb;
      if( trans1 == 'N' ) {
        M = one.rows();
        lda = M;
      }
      else {
        M = one.cols();
      }
      if( trans2 == 'N' ) {
        N = two.cols();
        K = two.rows();
        ldb = K;
      }
      else {
        N = two.rows();
        K = two.cols();
        ldb = N;
      }
      if( trans1 != 'N')
        lda = K;

      string nam = one.name();
      if( trans1 == 'T')
        nam += "^T";
      nam += "*"+two.name();
      if( trans2 == 'T' )
        nam += "^T";
      Matrix result(nam,M,N,false);

      BLASRPFX(gemm,BLASCHAR(trans1),BLASCHAR(trans2),BLASINT(M),BLASINT(N),BLASINT(K),BLASREAL(scaleA),BLASPREAL(one.data()[0]),BLASINT(lda),BLASPREAL(two.data()[0]),BLASINT(ldb),BLASREAL(mnlRealZero),BLASPREAL(result.data()[0]),BLASINT(M));

      return( result );
    }

    const Vector multTranspose(const Matrix& lhs, const Vector& rhs, char trans, Real alpha)
    {
      int N=lhs.rows();
      if( trans == 'T')
        N = lhs.cols();
      Vector result(lhs.name()+"*"+rhs.name(),N,false);
      BLASRPFX(gemv,BLASCHAR(trans),BLASINT(lhs.rows()),BLASINT(lhs.cols()),BLASREAL(mnlRealOne),const_cast<Real*>(lhs.data()[0]),BLASINT(lhs.rows()),const_cast<Real*>(rhs.data()),BLASINT(mnlIntOne),BLASREAL(mnlRealZero),result.data(),BLASINT(mnlIntOne));

      return( result );
    }

    void multTranspose(Matrix& result, const Matrix& one, const Matrix& two, char trans1, char trans2, Real scaleA, Real scaleB, int skipcol)
    {
      int M,N,K,lda,ldb;
      if( trans1 == 'N' ) {
        M = one.rows();
        lda = M;
      }
      else {
        M = one.cols()-skipcol;
      }
      if( trans2 == 'N' ) {
        N = two.cols();
        K = two.rows();
        ldb = K;
      }
      else {
        N = two.rows();
        K = two.cols();
        ldb = N;
      }
      if( trans1 != 'N')
        lda = K;

      assert( result.rows()-skipcol == M && result.cols() == N );

      BLASRPFX(gemm,BLASCHAR(trans1),BLASCHAR(trans2),BLASINT(M),BLASINT(N),BLASINT(K),BLASREAL(scaleA),BLASPREAL(one.data()[skipcol]),BLASINT(lda),BLASPREAL(two.data()[0]),BLASINT(ldb),BLASREAL(scaleB),BLASPREAL(result.data()[skipcol]),BLASINT(M));
    }

    const complexMatrix multTranspose(const complexMatrix& one, const complexMatrix& two, char trans1, char trans2, const std::complex<Real>& scaleA)
    {
      int M,N,K,lda,ldb;
      if( trans1 == 'N' ) {
        M = one.rows();
        lda = M;
      }
      else {
        M = one.cols();
      }
      if( trans2 == 'N' ) {
        N = two.cols();
        K = two.rows();
        ldb = K;
      }
      else {
        N = two.rows();
        K = two.cols();
        ldb = N;
      }
      if( trans1 != 'N')
        lda = K;

      string nam = one.name();
      if( trans1 == 'T')
        nam += "^T";
      nam += "*"+two.name();
      if( trans2 == 'T' )
        nam += "^T";

      complexMatrix result(nam,M,N,false);
      BLASCPFX(gemm,BLASCHAR(trans1),BLASCHAR(trans2),BLASINT(M),BLASINT(N),BLASINT(K),BLASCOMPLEX(scaleA),BLASPCOMPLEX(one.data()[0]),BLASINT(lda),BLASPCOMPLEX(two.data()[0]),BLASINT(ldb),BLASCOMPLEX(mnlComplexZero),BLASPCOMPLEX(result.data()[0]),BLASINT(M));

      return( result );
    }

    void multTranspose(complexMatrix& result, const complexMatrix& one, const complexMatrix& two, char trans1, char trans2, const std::complex<Real>& scaleA, const std::complex<Real>& scaleB)
    {
      int M,N,K,lda,ldb;
      if( trans1 == 'N' ) {
        M = one.rows();
        lda = M;
      }
      else {
        M = one.cols();
      }
      if( trans2 == 'N' ) {
        N = two.cols();
        K = two.rows();
        ldb = K;
      }
      else {
        N = two.rows();
        K = two.cols();
        ldb = N;
      }
      if( trans1 != 'N')
        lda = K;

      assert( result.rows() == M && result.cols() == N );

      BLASCPFX(gemm,BLASCHAR(trans1),BLASCHAR(trans2),BLASINT(M),BLASINT(N),BLASINT(K),BLASCOMPLEX(scaleA),BLASPCOMPLEX(one.data()[0]),BLASINT(lda),BLASPCOMPLEX(two.data()[0]),BLASINT(ldb),BLASCOMPLEX(scaleB),BLASPCOMPLEX(result.data()[0]),BLASINT(M));
    }

    const Matrix operator -(const Matrix& one, const Matrix& two)
    {
      assert( one.rows() == two.rows() && one.cols() == two.cols() );

      Matrix result(one);
      int N=one.rows()*one.cols();
      BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealMinusOne),const_cast<Real*>(two.data()[0]),BLASINT(mnlIntOne),result.data()[0],BLASINT(mnlIntOne));

      return result;
    }

    const Matrix operator +(const Matrix& one, const Matrix& two)
    {
      assert( one.rows() == two.rows() && one.cols() == two.cols() );

      Matrix result(one);
      int N=one.rows()*one.cols();
      BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealOne),const_cast<Real*>(two.data()[0]),BLASINT(mnlIntOne),result.data()[0],BLASINT(mnlIntOne));

      return result;
    }

    const Matrix operator *(Real alpha, const Matrix& one)
    {
      Matrix result(one);
      result *= alpha;

      return result;
    }

    const complexMatrix operator*(const complexMatrix& one, const complexMatrix& two)
    {
      assert( one.cols() == two.rows() );
      complexMatrix result(one.name()+"*"+two.name(),one.rows(),two.cols(),false);

      BLASCPFX(gemm,BLASCHAR(mnlNoTrans),BLASCHAR(mnlNoTrans),BLASINT(one.rows()),BLASINT(two.cols()),BLASINT(two.rows()),BLASCOMPLEX(mnlComplexOne),BLASPCOMPLEX(one.data()[0]),BLASINT(one.rows()),BLASPCOMPLEX(two.data()[0]),BLASINT(two.rows()),BLASCOMPLEX(mnlComplexZero),BLASPCOMPLEX(result.data()[0]),BLASINT(one.rows()));

      return( result );
    }

    const complexMatrix operator *(const complexMatrix& one, const std::complex<Real> factor)
    {
      complexMatrix result(one);
      result *= factor;

      return( result );
    }

    const Matrix kron(const Matrix& A, const Matrix& B)
    {
      Matrix result("kron("+A.name()+","+B.name()+")",A.rows()*B.rows(),A.cols()*B.cols());

      for( int i=0;i<A.cols();++i )
        for( int j=0;j<A.rows();++j )
          for( int k=0;k<B.cols();++k )
            for( int l=0;l<B.rows();++l )
              result[j*B.rows()+k][i*B.cols()+l] = A[i][j]*B[k][l];

      return result;
    }

    void multPointwise(complexMatrix& result, const complexMatrix& phi, const complexMatrix& u)
    {
      assert( (result.rows() == phi.rows()) && (result.cols() == phi.cols()) && (result.rows() == u.rows()) && (result.cols() == u.cols()) );

      transform(u.data()[0],u.data()[0]+u.rows()*u.cols(),phi.data()[0],result.data()[0],multiplies< std::complex<Real> >());
    }

    void multPointwise(Matrix& result, const Matrix& phi, const Matrix& u)
    {
      assert( result.rows() == phi.rows() && phi.rows() == u.rows() && result.cols() == phi.cols() && phi.cols() == u.cols() );

      int M = u.length();
      for( int i=0;i<M;++i )
        result.data()[0][i] = u.data()[0][i]*phi.data()[0][i];
    }

    void multPointwise(Matrix& result, const Matrix& u)
    {
      int M = u.length();
      for( int i=0;i<M;++i )
        result.data()[0][i] *= u.data()[0][i];
    }

    void multPointwise(Matrix& result, const Matrix& phi, const Matrix& u, Real alpha, Real beta)
    {
      assert( result.rows() == phi.rows() && phi.rows() == u.rows() && result.cols() == phi.cols() && phi.cols() == u.cols() );

      int M = u.length();
      for( int i=0;i<M;++i )
        result.data()[0][i] = alpha*result.data()[0][i]+beta*u.data()[0][i]*phi.data()[0][i];
    }

    void multPointwiseReal(complexMatrix& result, const complexMatrix& phi, const complexMatrix& u)
    { 
      assert( (result.rows() == phi.rows()) && (result.cols() == phi.cols()) && (result.rows() == u.rows()) && (result.cols() == u.cols()) );

      const Real* uData = BLASPREAL(u.data()[0]);
      const Real* pData = BLASPREAL(phi.data()[0]);
      std::complex<Real>* rData = result.data()[0];
      for( int i=0;i<result.rows()*result.cols();++i ) {
        rData[i] = (*uData)*(*pData);
        uData +=2; pData +=2;
      }
    }

    void multPointwiseReal(complexMatrix& phi, const complexMatrix& u)
    {
      assert( (phi.rows() == u.rows()) && (phi.cols() == u.cols()) );

      const Real* uData = BLASPREAL(u.data()[0]);
      Real* pData = BLASPREAL(phi.data()[0]);
      for( int i=0;i<phi.rows()*phi.cols();++i ) {
        *pData *= (*uData);
        uData += 2; pData += 2;
      }
    }			 

    const Matrix loadMatrix(const string& filename)
    {
      int n_rows;
      int n_cols;

      stringstream s;
      string line;
      ifstream f;
      f.open(filename.c_str());
      getline(f, line); // read dimensions
      line.erase(0,2);
      int pos = line.rfind('\t');
      string name;
      if( pos != string::npos ) {
        name = line.substr(pos+1);
        line.erase(pos,line.size()-pos);
      }
      s << line;
      s >> n_rows;
      s >> n_cols;
      f.close();

      Matrix result(name,n_rows,n_cols,false);
      result.load(filename);

      return result;
    }

  } // namespace basics
} // namespace mnl

