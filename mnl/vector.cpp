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

#include "vector.h"
#include "matrix.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

namespace mnl {
  namespace basics {
    Vector::Vector(const string& n_name, int n_length, bool clear, Real* n_data) : 
      basicVector<Real>(n_name, n_length, clear, n_data)
    {
    }

    Vector::Vector(const Vector& vector) : 
      basicVector<Real>(vector)
    {
      *this = vector;
    }

    Vector::Vector(const string& n_name, const Vector& vector) :
      basicVector<Real>(n_name,vector.length(),false,const_cast<Real*>(vector.data()))
    {
    }

    Vector::~Vector()
    {
    }

    void Vector::print() const
    {
      cout << v_name << "(" << length() << ") = " << "[";
      for( int i=0;i<length()-1;++i ) {
        if( fabs(data()[i]) < 1.e-14f )
          cout << data()[i] << " ";
        else
          cout << data()[i] << " ";
      }
      if( fabs(data()[length()-1]) < 1.e-14f )
        cout << 0 << "]" << endl;
      else
        cout << data()[length()-1] << "]" << endl;
    }

    void Vector::save(const string& filename) const
    {
      ofstream f;
      f.open(filename.c_str());
      f.precision(16);
      f << fixed;
      f << "% " << length() << '\t' << name() << endl;
      for( size_t i=0;i<length();++i )
        f << v_data[i] << endl;
      f.close();
    }

    void Vector::load(const string& filename)
    {
      int n;
      stringstream s;
      string line;
      ifstream f;
      f.open(filename.c_str());
      getline(f, line); // read dimensions
      line.erase(0,2);
      int pos=line.find('\t');
      if( pos != string::npos )
        line.erase(pos,line.size()-pos);

      s << line;
      s >> n;

      assert( n == length() );

      for( size_t i=0;i<n;++i ) { // read one line
        getline(f, line); // read actual data
        s.clear();
        s << line;
        Real t;
        s >> t;
        v_data[i] = t;
      }
      f.close();
    }

    void Vector::axpy(const Real alpha, const Vector& x)
    {
      assert( x.length() == length() );

      BLASRPFX(axpy,BLASINT(m_N),BLASREAL(alpha),BLASPREAL(x.data()),BLASINT(mnlIntOne),v_data,BLASINT(mnlIntOne));
    }

    void Vector::gemv(const Vector& x, const Matrix& A, const char trans, const Real alpha, const Real beta)
    {
      BLASRPFX(gemv,BLASCHAR(trans),BLASINT(A.rows()),BLASINT(A.cols()),BLASREAL(alpha),BLASPREAL(A.data()[0]),BLASINT(A.rows()),BLASPREAL(x.data()),BLASINT(mnlIntOne),BLASREAL(beta),BLASPREAL(data()),BLASINT(mnlIntOne));
    }

    Real Vector::dot(const basicVector<Real>& v2) const
    {
      assert( length() == v2.length() );

      Real result;
      result = BLASRPFX(dot,BLASINT(m_N),v_data,BLASINT(mnlIntOne),BLASPREAL(v2.data()),BLASINT(mnlIntOne));

      return( result );	
    } 

    Vector Vector::invert() const
    {
      Vector result("Inverse of "+name(),length(),false);
      for( int i=0;i<length();++i )
        result[i] = Real(1)/v_data[i];

      return( result );
    } 

    Vector& Vector::operator =(const Vector& v2)
    {
      assert( length() == v2.length() );

      BLASRPFX(copy,BLASINT(m_N),BLASPREAL(v2.data()),BLASINT(mnlIntOne),v_data,BLASINT(mnlIntOne));

      return *this;
    } 

    Vector& Vector::operator =(const Real *v2)
    {
      BLASRPFX(copy,BLASINT(m_N),BLASPREAL(v2),BLASINT(mnlIntOne),BLASPREAL(v_data),BLASINT(mnlIntOne));

      return *this;
    } 

    Vector& Vector::operator =(const Real val)
    {
      for( int i=0;i<length();++i )
        v_data[i] = val;

      return *this;
    }

    void Vector::operator +=(const Vector& v2)
    {
      assert( v2.length()==length() ); 

      BLASRPFX(axpy,BLASINT(m_N),BLASREAL(mnlRealOne),BLASPREAL(v2.data()),BLASINT(mnlIntOne),v_data,BLASINT(mnlIntOne));
    } 

    void Vector::operator +=(const Real* v2)
    {
      BLASRPFX(axpy,BLASINT(m_N),BLASREAL(mnlRealOne),BLASPREAL(v2),BLASINT(mnlIntOne),BLASPREAL(v_data),BLASINT(mnlIntOne));
    }

    void Vector::operator -=(const Vector& v2)
    {
      assert( v2.length()==length() ); 

      BLASRPFX(axpy,BLASINT(m_N),BLASREAL(mnlRealMinusOne),BLASPREAL(v2.data()),BLASINT(mnlIntOne),v_data,BLASINT(mnlIntOne));
    }

    void Vector::operator -=(const Real* v2)
    {
      BLASRPFX(axpy,BLASINT(m_N),BLASREAL(mnlRealMinusOne),BLASPREAL(v2),BLASINT(mnlIntOne),BLASPREAL(v_data),BLASINT(mnlIntOne));
    }


    void Vector::operator *=(const Real factor)
    {
      BLASRPFX(scal,BLASINT(m_N),BLASREAL(factor),v_data,BLASINT(mnlIntOne));
    } 


    void Vector::operator /=(const Real factor)
    {
      *this *= Real(1)/factor;
    } 

    Vector Vector::operator[](const utilities::Range& range) const
    {
      string s("Range from ");
      s += v_name;
      Vector result(s,range.size(),false);
      for( int i=0;i<range.size();++i )
        result[i] = v_data[range[i]];

      return( result );
    } 

    complexVector::complexVector(const string& n_name, int n_length, bool clear, std::complex<Real>* n_data) : 
      basicVector< std::complex<Real> >(n_name, n_length, clear, n_data)
    {
    }

    complexVector::complexVector(const complexVector& vector) : 
      basicVector< std::complex<Real> >(vector)
    {
      *this = vector;
    } 

    complexVector::complexVector(const Vector& vector) : 
      basicVector< std::complex<Real> >(vector.name(),vector.length())
    {
      m_N=vector.length();
#ifdef MNL_MEMORY_VERBOSE_
      m_memId = utilities::g_tracker.addChunk(utilities::memTracker::Vector,m_N*sizeof(std::complex<Real>));
#endif
      *this = vector;
    } 

    complexVector::complexVector(const string& n_name, int n_length, const Real* re, const Real* im) : 
      basicVector< std::complex<Real> >(n_name, n_length)
    {			
#ifdef MNL_MEMORY_VERBOSE_
      m_memId = utilities::g_tracker.addChunk(utilities::memTracker::Vector,m_N*sizeof(std::complex<Real>));
#endif
      for( int i=0;i<m_N;++i )
        v_data[i] = std::complex<Real>(re[i],im[i]);
    }

    complexVector::complexVector(const string& n_name, const complexVector& vector) :
      basicVector< std::complex<Real> >(n_name,vector.length(),false,const_cast<std::complex<Real>*>(vector.data()))
    {
    }

    complexVector::~complexVector()
    { 
    }

    void complexVector::print() const
    {
      Real re,im;
      cout << v_name << " = " << "[";
      cout.precision(5);
      for( int i=0;i<length()-1;++i ) {
        if( fabs(data()[i].real()) < 1.e-14f )
          re = 0;
        else
          re = data()[i].real();
        if( fabs(data()[i].imag()) < 1.e-14f )
          im = 0;
        else
          im = data()[i].imag();
        if( im >= 0 )
          cout << re << "+" << im << "i ";
        else
          cout << re << im << "i ";
      }
      if( fabs(data()[length()-1].real()) < 1.e-14f )
        re = 0;
      else
        re = data()[length()-1].real();
      if( fabs(data()[length()-1].imag()) < 1.e-14f )
        im = 0;
      else
        im = data()[length()-1].imag();
      if( im >= 0 )
        cout << re << "+" << im << "i]" << endl;
      else
        cout << re << im << "i]" << endl;
    } 

    void complexVector::axpy(const std::complex<Real> alpha, const complexVector& x)
    {
      assert( x.length() == length() );

      BLASCPFX(axpy,BLASINT(m_N),BLASCOMPLEX(alpha),BLASPCOMPLEX(x.data()),BLASINT(mnlIntOne),BLASPCOMPLEX(v_data),BLASINT(mnlIntOne));
    }

    void complexVector::axpy(const Real alpha, const complexVector& x)
    {
      assert( x.length() == length() );

      int N=2*length();
      BLASRPFX(axpy,BLASINT(N),BLASREAL(alpha),BLASPREAL(x.data()),BLASINT(mnlIntOne),BLASPREAL(v_data),BLASINT(mnlIntOne));
    }

    void complexVector::axpy(const Real alpha, const Vector& x)
    {
      assert( x.length() == length() );

      BLASRPFX(axpy,BLASINT(m_N),BLASREAL(alpha),BLASPREAL(x.data()),BLASINT(mnlIntOne),BLASPREAL(v_data),BLASINT(mnlIntTwo));
    }

    void complexVector::gemv(const complexVector& x, const Matrix& A, const char trans, const Real alpha, const Real beta)
    {
      BLASRPFX(gemv,BLASCHAR(trans),BLASINT(A.rows()),BLASINT(A.cols()),BLASREAL(alpha),BLASPREAL(A.data()[0]),BLASINT(A.rows()),BLASPREAL(x.data()),BLASINT(mnlIntTwo),BLASREAL(beta),BLASPREAL(data()),BLASINT(mnlIntTwo));
      BLASRPFX(gemv,BLASCHAR(trans),BLASINT(A.rows()),BLASINT(A.cols()),BLASREAL(alpha),BLASPREAL(A.data()[0]),BLASINT(A.rows()),BLASPREAL(x.data())+1,BLASINT(mnlIntTwo),BLASREAL(beta),BLASPREAL(data())+1,BLASINT(mnlIntTwo));
    }

    void complexVector::gemv(const complexVector& x, const complexMatrix& A, const char trans, const std::complex<Real> alpha, const std::complex<Real> beta)
    {
      BLASCPFX(gemv,BLASCHAR(trans),BLASINT(A.rows()),BLASINT(A.cols()),BLASCOMPLEX(alpha),BLASPCOMPLEX(A.data()[0]),BLASINT(A.rows()),BLASPCOMPLEX(x.data()),BLASINT(mnlIntOne),BLASCOMPLEX(beta),BLASPCOMPLEX(data()),BLASINT(mnlIntOne));
    }

    void complexVector::gemv(const complexMatrix& x, const int row, const complexMatrix& A, const char trans, const std::complex<Real> alpha, const std::complex<Real> beta)
    {
      BLASCPFX(gemv,BLASCHAR(trans),BLASINT(A.rows()),BLASINT(A.cols()),BLASCOMPLEX(alpha),BLASPCOMPLEX(A.data()[0]),BLASINT(A.rows()),BLASPCOMPLEX(x.data()[0]+row),BLASINT(x.rows()),BLASCOMPLEX(beta),BLASPCOMPLEX(data()),BLASINT(mnlIntOne));
    }

    Vector complexVector::real() const
    {
      Vector result("Real part of "+name(),length(),false);

      BLASRPFX(copy,BLASINT(m_N),BLASPREAL(v_data),BLASINT(mnlIntTwo),BLASPREAL(result.data()),BLASINT(mnlIntOne));

      return( result );
    } 

    Vector complexVector::imag() const
    {
      Vector result("Imaginary part of "+name(),length(),false);

      BLASRPFX(copy,BLASINT(m_N),BLASPREAL(v_data)+1,BLASINT(mnlIntTwo),BLASPREAL(result.data()),BLASINT(mnlIntOne));

      return( result );
    } 

    std::complex<Real> complexVector::dot(const basicVector<std::complex<Real> >& v2) const
    {
      assert( length() == v2.length() );

      std::complex<Real> result = 0;
      Complex res;

      ZDOTC(res,BLASINT(m_N),BLASPCOMPLEX(v_data),BLASINT(mnlIntOne),BLASPCOMPLEX(v2.data()),BLASINT(mnlIntOne));
      memcpy(&result,&res,2*sizeof(Real));

      return( result );	
    } 

    complexVector& complexVector::operator =(const complexVector& v2)
    {
      assert( length() == v2.length() );

      BLASCPFX(copy,BLASINT(m_N),BLASPCOMPLEX(v2.data()),BLASINT(mnlIntOne),BLASPCOMPLEX(v_data),BLASINT(mnlIntOne));

      return *this;
    } 

    complexVector& complexVector::operator =(const Vector& v2)
    {
      assert( length() == v2.length() );

      int inc=1,inc2=2;
      memset(v_data,0,length()*sizeof(std::complex<Real>));
      BLASRPFX(copy,BLASINT(m_N),BLASPREAL(v2.data()),BLASINT(mnlIntOne),BLASPREAL(v_data),BLASINT(mnlIntTwo));

      return *this;
    } 

    complexVector& complexVector::operator =(const std::complex<Real>* v2)
    {
      BLASCPFX(copy,BLASINT(m_N),BLASPCOMPLEX(v2),BLASINT(mnlIntOne),BLASPCOMPLEX(v_data),BLASINT(mnlIntOne));

      return *this;
    } 

    void complexVector::operator +=(const complexVector& v2)
    {
      assert( v2.length()==length() ); 

      BLASCPFX(axpy,BLASINT(m_N),BLASCOMPLEX(mnlComplexOne),BLASPCOMPLEX(v2.data()),BLASINT(mnlIntOne),BLASPCOMPLEX(v_data),BLASINT(mnlIntOne));
    } 

    void complexVector::operator +=(const Vector& v2)
    {
      assert( length() == v2.length() );

      BLASRPFX(axpy,BLASINT(m_N),BLASREAL(mnlRealOne),BLASPREAL(v2.data()),BLASINT(mnlIntOne),BLASPREAL(v_data),BLASINT(mnlIntTwo));
    } 

    void complexVector::operator +=(const std::complex<Real>* v2)
    {
      BLASCPFX(axpy,BLASINT(m_N),BLASCOMPLEX(mnlComplexOne),BLASPCOMPLEX(v2),BLASINT(mnlIntOne),BLASPCOMPLEX(v_data),BLASINT(mnlIntOne));
    } 

    void complexVector::operator -=(const complexVector& v2)
    {
      assert( v2.length()==length() ); 

      BLASCPFX(axpy,BLASINT(m_N),BLASCOMPLEX(mnlComplexMinusOne),BLASPCOMPLEX(v2.data()),BLASINT(mnlIntOne),BLASPCOMPLEX(v_data),BLASINT(mnlIntOne));
    } 

    void complexVector::operator -=(const Vector& v2)
    { 
      assert( v2.length()==length() ); 

      BLASRPFX(axpy,BLASINT(m_N),BLASREAL(mnlRealMinusOne),BLASPREAL(v2.data()),BLASINT(mnlIntOne),BLASPREAL(v_data),BLASINT(mnlIntTwo));
    }

    void complexVector::operator -=(const std::complex<Real>* v2)
    {
      BLASCPFX(axpy,BLASINT(m_N),BLASCOMPLEX(mnlComplexMinusOne),BLASPCOMPLEX(v2),BLASINT(mnlIntOne),BLASPCOMPLEX(v_data),BLASINT(mnlIntOne));
    }

    void complexVector::operator *=(const Vector& v2)
    {
      assert( v2.length()==length() ); 

      for( int i=0;i<length();++i ) // BLASALIZE
        v_data[i] *= v2[i];
    } 

    void complexVector::operator *=(const std::complex<Real> factor)
    {
      BLASCPFX(scal,BLASINT(m_N),BLASCOMPLEX(factor),BLASPCOMPLEX(v_data),BLASINT(mnlIntOne));
    } 

    void complexVector::operator *=(const Real factor)
    {
      int N2 = 2*m_N;
      BLASRPFX(scal,BLASINT(N2),BLASREAL(factor),BLASPREAL(v_data),BLASINT(mnlIntOne));
    }

    void complexVector::operator /=(const Real& factor)
    {
      *this *= Real(1)/factor;
    } 

    const Vector loadVector(const string& filename)
    {
      int n;
      stringstream s;
      string line;
      ifstream f;
      f.open(filename.c_str());
      getline(f, line); // read dimensions
      line.erase(0,2);
      int pos=line.find('\t');
      string name;
      if( pos != string::npos ) {
        name = line.substr(pos+1);
        line.erase(pos,line.size()-pos);
      }
      s << line;
      s >> n;

      Vector result(name,n);
      result.load(filename);

      return result;
    }

    complexVector complexVector::operator[](const utilities::Range& range) const
    {
      string s("Range from ");
      s += v_name;
      complexVector result(s,range.size(),false);
      for( int i=0;i<range.size();++i )
        result[i] =  v_data[range[i]];

      return( result );
    }

    const Vector operator +(const Vector& one, Real two)
    {
      Vector result(one);
      result += two;

      return( result );
    }

    const Vector operator -(const Vector& one, Real two)
    {
      Vector result(one);
      result -= two;

      return( result );
    }

    const Vector operator*(const Vector& one, Real two)
    {
      Vector result(one);
      result *= two;

      return( result );
    }

    const Vector operator/(const Vector& one, Real two)
    {
      Vector result(one);
      result /= two;

      return( result );
    }
  } // namespace basics
} // namespace mnl

