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

#include "matrices.h"
#ifdef MNL_MEMORY_VERBOSE_
#include "memtracker.h"
#endif

#include <sstream>
#include <algorithm>
#include <functional>

#ifdef __INTEL_COMPILER
#include <xmmintrin.h>
#endif

using namespace std;

namespace mnl {
  namespace basics {
    Matrices::Matrices(const string& n_name, int n_rows, int n_cols, int n_k, bool doClear, Real* n_data) : 
      basicMatrices<Real,Matrix>(n_name,n_rows,n_cols,n_k,doClear,n_data)
    {
    }

    Matrices::Matrices(const Matrices& matrices) : 
      basicMatrices<Real,Matrix>(matrices)
    {
      *this = matrices;
    }

    Matrices::Matrices(const string& n_name, const Matrices& matrices) : 
      basicMatrices<Real,Matrix>(n_name,matrices.rows(),matrices.cols(),matrices.matrices(),false,const_cast<Real*>(matrices.data()))
    {
    }

    Matrices::~Matrices()
    {
    }

    Matrices& Matrices::operator=(const Matrices& matrices)
    {
      assert( (rows() == matrices.rows()) && (cols() == matrices.cols()) && (this->matrices() == matrices.matrices()) );

      int N=rows()*cols()*this->matrices();
      BLASRPFX(copy,BLASINT(N),BLASPREAL(matrices.data()),BLASINT(mnlIntOne),BLASPREAL(the_data),BLASINT(mnlIntOne));

      return( *this );
    }

    Matrices& Matrices::operator=(const Real in)
    {
      for( int i=0;i<length();++i )
        the_data[i] = in;	
      return( *this );
    }

    void Matrices::operator +=(const Matrices& matrices)
    {
      assert( (m_k == matrices.matrices()) && (m_rows == matrices.rows()) && (m_cols == matrices.cols()) );

      int N=m_k*m_rows*m_cols;
      BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealOne),BLASPREAL(matrices.data()),BLASINT(mnlIntOne),BLASPREAL(the_data),BLASINT(mnlIntOne));
    }

    void Matrices::operator +=(const Real alpha)
    {
      for( int i=0;i<length();++i )
        the_data[i] += alpha;
    }

    void Matrices::operator -=(const Matrices& matrices)
    {
      assert( (m_k == matrices.matrices()) && (m_rows == matrices.rows()) && (m_cols == matrices.cols()) );

      int N=m_k*m_rows*m_cols;
      BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealMinusOne),BLASPREAL(matrices.data()),BLASINT(mnlIntOne),BLASPREAL(the_data),BLASINT(mnlIntOne));
    }

    void Matrices::operator -=(const Real alpha)
    {
      for( int i=0;i<length();++i )
        the_data[i] -= alpha;
    }

    void Matrices::operator*=(const Real alpha)
    {
      int N=m_k*m_rows*m_cols;
      BLASRPFX(scal,BLASINT(N),BLASREAL(alpha),the_data,BLASINT(mnlIntOne));
    }

    void Matrices::axpy(const Real alpha, const Matrices& x)
    {
      assert( (rows() == x.rows()) && (cols() == x.cols()) && (matrices() == x.matrices()) );

      int N=m_k*m_rows*m_cols;
      BLASRPFX(axpy,BLASINT(N),BLASREAL(alpha),BLASPREAL(x.data()),BLASINT(mnlIntOne),BLASPREAL(the_data),BLASINT(mnlIntOne));
    }

    Real Matrices::dot(const basicMatrices<Real,Matrix>& x) const
    {
      assert( (rows() == x.rows()) && (cols() == x.cols()) && (matrices() == x.matrices()) );
      int N=cols()*rows()*matrices();
      return BLASRPFX(dot,BLASINT(N),BLASPREAL(x.data()),BLASINT(mnlIntOne),BLASPREAL(data()),BLASINT(mnlIntOne)); 

    }

    complexMatrices::complexMatrices(const string& n_name, int n_rows, int n_cols, int n_k, bool clear) : 
      basicMatrices<std::complex<Real>,complexMatrix>(n_name,n_rows,n_cols,n_k,clear)
    {
    }

    complexMatrices::complexMatrices(const complexMatrices& matrices) : 
      basicMatrices<std::complex<Real>,complexMatrix>(matrices)
    {
      *this = matrices;
    }

    complexMatrices::complexMatrices(const string& name, const complexMatrices& matrices) : 
      basicMatrices<std::complex<Real>,complexMatrix>(name,matrices.rows(),matrices.cols(),matrices.matrices(),false,const_cast<std::complex<Real>* >(matrices.data()))
    {
    }

    complexMatrices::~complexMatrices()
    {
    }

    Matrices complexMatrices::real() const
    { 
      Matrices result("Real part of "+name(),rows(),cols(),matrices(),false);

      for( int i=0;i<matrices();++i )
        result[i] = (*this)[i].real();

      return( result );
    }

    Matrices complexMatrices::imag() const
    {
      Matrices result("Imaginary part of "+name(),rows(),cols(),matrices(),false);
      for( int i=0;i<matrices();++i )
        result[i] = (*this)[i].imag();

      return( result );
    } 

    void complexMatrices::axpy(const std::complex<Real> alpha, const complexMatrices& x)
    {
      assert( (rows() == x.rows()) && (cols() == x.cols()) && (matrices() == x.matrices()) );

      int N=m_k*m_rows*m_cols;
      BLASCPFX(axpy,BLASINT(N),BLASCOMPLEX(alpha),BLASPCOMPLEX(x.data()),BLASINT(mnlIntOne),BLASPCOMPLEX(the_data),BLASINT(mnlIntOne));
    }

    void complexMatrices::axpy(const Real alpha, const complexMatrices& x)
    {
      assert( (rows() == x.rows()) && (cols() == x.cols()) && (matrices() == x.matrices()) );

      int N=m_k*m_rows*m_cols*2;
      BLASRPFX(axpy,BLASINT(N),BLASREAL(alpha),BLASPREAL(x.data()),BLASINT(mnlIntOne),BLASPREAL(the_data),BLASINT(mnlIntOne));
    }

    std::complex<Real> complexMatrices::dot(const basicMatrices<std::complex<Real>,complexMatrix>& x) const
    {
      assert( cols() == x.cols() && rows() == x.rows() && matrices() == x.matrices() );

      int N=cols()*rows()*matrices();

      std::complex<Real> result;
#ifndef __freebsd__
      ZDOTC(result,BLASINT(N),BLASPCOMPLEX(x.data()),BLASINT(mnlIntOne),BLASPCOMPLEX(data()),BLASINT(mnlIntOne));
#else
      for( int j=0;j<N;++j )
        result += conj(*(x.data()+j))*(*(data()+j));
#endif
      return result;
    }

    complexMatrices& complexMatrices::operator=(const complexMatrices& matrices)
    {
      assert( (rows() == matrices.rows()) && (cols() == matrices.cols()) && (this->matrices() == matrices.matrices()) );

      int N=rows()*cols()*this->matrices();
      BLASCPFX(copy,BLASINT(N),BLASPCOMPLEX(matrices.data()),BLASINT(mnlIntOne),BLASPCOMPLEX(the_data),BLASINT(mnlIntOne));

      return( *this );
    }

    complexMatrices& complexMatrices::operator=(const Matrices& matrices)
    {
      assert( (rows() == matrices.rows()) && (cols() == matrices.cols()) && (this->matrices() == matrices.matrices()) );

      int N=rows()*cols()*this->matrices();
      memset(the_data,0,N*sizeof(std::complex<Real>));
      BLASRPFX(copy,BLASINT(N),BLASPREAL(matrices.data()),BLASINT(mnlIntOne),BLASPREAL(the_data),BLASINT(mnlIntTwo));

      return( *this );
    }

    void complexMatrices::operator -=(const Real mean)
    {
      int size = rows()*cols();
      for( int l=0;l<matrices();++l ) {
        for( int i=0;i<size;++i )
          *(data()+l*size+i) -= mean;
      }
    }

    void complexMatrices::operator -=(const complexMatrices& matrices)
    {
      assert( (m_k == matrices.matrices()) && (m_rows == matrices.rows()) && (m_cols == matrices.cols()) ); // maybe throw exception instead.. dunno. in release this SHOULD not happen anyway.

      int N=m_k*m_rows*m_cols*2;
      BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealMinusOne),BLASPREAL(matrices.data()),BLASINT(mnlIntOne),BLASPREAL(the_data),BLASINT(mnlIntOne));
    }

    void complexMatrices::operator -=(const Matrices& matrices)
    {
      assert( (m_k == matrices.matrices()) && (m_rows == matrices.rows()) && (m_cols == matrices.cols()) ); // maybe throw exception instead.. dunno. in release this SHOULD not happen anyway.

      int N=m_k*m_rows*m_cols;
      BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealMinusOne),BLASPREAL(matrices.data()),BLASINT(mnlIntOne),BLASPREAL(the_data),BLASINT(mnlIntTwo));
    }

    void complexMatrices::operator +=(const complexMatrices& matrices)
    {
      assert( (m_k == matrices.matrices()) && (m_rows == matrices.rows()) && (m_cols == matrices.cols()) ); // maybe throw exception instead.. dunno. in release this SHOULD not happen anyway.

      int N=m_k*m_rows*m_cols*2;
      BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealOne),BLASPREAL(matrices.data()),BLASINT(mnlIntOne),BLASPREAL(the_data),BLASINT(mnlIntOne));
    }

    void complexMatrices::operator +=(const Matrices& matrices)
    {
      assert( (m_k == matrices.matrices()) && (m_rows == matrices.rows()) && (m_cols == matrices.cols()) ); // maybe throw exception instead.. dunno. in release this SHOULD not happen anyway.

      int N=m_k*m_rows*m_cols;
      BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealOne),BLASPREAL(matrices.data()),BLASINT(mnlIntOne),BLASPREAL(the_data),BLASINT(mnlIntTwo));
    }

    void complexMatrices::operator *=(const std::complex<Real> alpha)
    {
      int N=m_k*m_rows*m_cols;
      BLASCPFX(scal,BLASINT(N),BLASCOMPLEX(alpha),BLASPCOMPLEX(the_data),BLASINT(mnlIntOne));
    }

    void complexMatrices::operator *=(const Real alpha)
    {
      int N=2*m_k*m_rows*m_cols;
      BLASRPFX(scal,BLASINT(N),BLASREAL(alpha),BLASPREAL(the_data),BLASINT(mnlIntOne));
    }

    complexMatrices complexMatrices::operator -(const complexMatrices& matrices) const
    {
      string n_name=name()+"-"+matrices.name();
      complexMatrices result(*this);
      result -= matrices;
      result.setName(n_name);

      return( result );
    }

    void multPointwise(complexMatrices& result, const complexMatrices& phi, const complexMatrices& u)
    {
      assert( (result.rows() == phi.rows()) && (result.cols() == phi.cols()) && (result.matrices() == phi.matrices()) && (result.rows() == u.rows()) && (result.cols() == u.cols()) && (result.matrices() == u.matrices()) );

      transform(u.data(),u.data()+u.rows()*u.cols()*u.matrices(),phi.data(),result.data(),multiplies<std::complex<Real> >());
    } 

    void multPointwise(Matrices& result, const Matrices& phi, const Matrices& u)
    {
      assert( (result.rows() == phi.rows()) && (result.cols() == phi.cols()) && (result.matrices() == phi.matrices()) && (result.rows() == u.rows()) && (result.cols() == u.cols()) && (result.matrices() == u.matrices()) );

      transform(u.data(),u.data()+u.rows()*u.cols()*u.matrices(),phi.data(),result.data(),std::multiplies<Real>());
    } 

    void multPointwise(Matrices& result, const Matrices& phi)
    {
      transform(result.data(),result.data()+result.rows()*result.cols()*result.matrices(),phi.data(),result.data(),std::multiplies<Real>());
    } 

    void multPointwiseReal(complexMatrices& result, const complexMatrices& phi, const complexMatrices& u)
    {
      assert( (result.rows() == phi.rows()) && (result.cols() == phi.cols()) && (result.matrices() == phi.matrices()) && (result.rows() == u.rows()) && (result.cols() == u.cols()) && (result.matrices() == u.matrices()) );

      for( int l=0;l<phi.matrices();++l ) {
        const Real* uData = BLASPREAL(u[l].data()[1]);
        const Real* pData = BLASPREAL(phi[l].data()[1]);
        Real *rData = BLASPREAL(result[l].data()[1]);
        for( int i=0;i<phi.rows()*(phi.cols()-2);++i ) {
          *rData = (*uData)*(*pData);
          *(rData+1) = 0;
          uData += 2; pData += 2; rData += 2;
        }
      }
    }			 

    void multPointwise(complexMatrices& phi, const complexMatrices& u)
    {
      assert( (phi.rows() == u.rows()) && (phi.cols() == u.cols()) && (phi.matrices() == u.matrices()) );

      transform(u.data(),u.data()+u.rows()*u.cols()*u.matrices(),phi.data(),phi.data(),std::multiplies<std::complex<Real> >());
    } 

    void multPointwiseReal(complexMatrices& phi, const complexMatrices& u)
    {
      assert( (phi.rows() == u.rows()) && (phi.cols() == u.cols()) && (phi.matrices() == u.matrices()) );

#pragma omp parallel for schedule(static)
      for( int l=0;l<phi.matrices();++l ) {
        const Real* uData = BLASPREAL(u[l].data()[1]);
        Real* pData = BLASPREAL(phi[l].data()[1]);
        for( int i=0;i<phi.rows()*(phi.cols()-2);++i ) {
          *pData *= (*uData);
          *(pData+1) = 0;
          uData += 2; pData += 2;
        }
      }
    }

    void applyLocalGlobal(basics::Matrices& result, const basics::Matrices& victim, const basics::Matrix& op, 
        char transA, char transB, int ofs, Real scaleA, Real scaleB)
    {
      int planes = result.matrices();
      if( ofs == 3 )
        planes -= 1;
      else
        planes -= ofs;
      int skip=0;
      if( ofs == 1 || ofs == 2 )
        skip = 1;

      int M,N,K,lda,ldb;
      if( transA == 'N' )
        lda = M = result.rows()*result.cols();
      else
        M = planes;

      if( transB == 'N' ) {
        N = op.cols();
        K = op.rows();
        ldb = K;
      }
      else {
        N = op.rows();
        K = op.cols();
        ldb = N;
      }
      if( transA != 'N')
        lda = K;

      BLASRPFX(gemm,BLASCHAR(transA),BLASCHAR(transB),BLASINT(M),BLASINT(N),BLASINT(K),BLASREAL(scaleA),BLASPREAL(victim[skip].data()[0]),BLASINT(lda),BLASPREAL(op.data()[0]),BLASINT(ldb),BLASREAL(scaleB),BLASPREAL(result[skip].data()[0]),BLASINT(M));
    }
  } // namespace basics
} // namespace mnl

