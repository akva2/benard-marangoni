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
#ifndef MNL_VECTOR_H_
#define MNL_VECTOR_H_

/** -- muu numerics library --
* - this is the vector class, including operator overloads
*/
#include "config.h"
#include "range.h"
#include "memtracker.h"

#include <complex>
#include <vector>
#include <memory>
#include <string>
#include <assert.h>
#include <algorithm>
#include <functional>

namespace mnl {
  namespace basics {
    /* forwards */
    class Matrix;
    class complexMatrix;

    /* basic vector class */
    template<class T>
      class basicVector {
        public:
          class vectorDotter {
            public:
              Real operator()(const basicVector<T>& A, 
                  const basicVector<T>& B) const
              {
                return A.dot(B);
              }
          };

          basicVector(const std::string& n_name, int n_length, bool clear=true, T* n_data=NULL) :
            v_name(n_name), 
            external_data(false), 
            m_N(n_length)
        {
          if( !n_data ) {
            ALLOCATELOG(name());
            v_data = new T[length()];
#ifdef MNL_MEMORY_VERBOSE_
            m_memId = utilities::g_tracker.addChunk(utilities::memTracker::Vector,length()*sizeof(T));
#endif
            if( clear )
              memset(v_data,0,length()*sizeof(T));
          } else {
            external_data = true;
            v_data = n_data;
          }
        }

          basicVector(const basicVector<T>& input) :
            v_name(input.name()),
            external_data(false),
            m_N(input.length())
        {
          ALLOCATELOG(name() + " (from other vector)");
          v_data = new T[length()];
#ifdef MNL_MEMORY_VERBOSE_
          m_memId = utilities::g_tracker.addChunk(utilities::memTracker::Vector,length()*sizeof(T));
#endif
        }

          virtual ~basicVector()
          {
            if( !external_data ) {
              //                        DEALLOCATELOG(name());
              delete[] v_data;
#ifdef MNL_MEMORY_VERBOSE_
              utilities::g_tracker.removeChunk(m_memId);
#endif
            }
          }

          const int& length() const
          {
            return( m_N );
          }

          utilities::Range range() const
          {
            return( utilities::Range::colon(0,length()-1) );
          }

          const std::string& name() const
          {
            return( v_name );
          }

          void setName(const std::string& n_name)
          {
            v_name = n_name;
          }

          void info() const
          {
            std::cout << v_name << ": " << m_N << std::endl;
          }

          T sum() const
          {
            T result=0;
            for( int i=0;i<length();++i )
              result += v_data[i];

            return( result );
          }

          void clear()
          {
            memset(v_data,0,length()*sizeof(T));
          }

          void fill(const T alpha)
          {
            for( int i=0;i<length();++i )
              v_data[i] = alpha;
          }

          inline const T* data() const
          {
            return( v_data );
          }

          inline T* data()
          {
            return( v_data );
          }

          void data(T* n_data) // DANGEROUS!
          {
            if( !external_data ) {
              DEALLOCATELOG(name() + " (assigning external data)");
              delete[] v_data;
            }
            v_data = n_data;
            external_data = true;
          } 

          virtual T dot(const basicVector<T>& in) const = 0;

          const vectorDotter& getDotter() const
          {
            return m_dotter;
          }

          inline const T& operator[](const int index) const
          {
            assert( index < length() );

            return( v_data[index] );
          }

          inline T& operator[](const int index)
          {
            assert( index < length() );

            return( v_data[index] );
          }

          virtual void operator +=(const T data)
          {
            for( int i=0;i<length();++i)
              v_data[i] += data;
          }

          virtual void operator -=(const T data)
          {
            for( int i=0;i<length();++i )
              v_data[i] -= data;
          }

          virtual void operator *=(const basicVector<T>& v2)
          {
            assert( v2.length()==length() ); 

            std::transform(data(),data()+length(),v2.data(),data(),std::multiplies< T >());
          }

          virtual void operator *=(const T* v2)
          {
            std::transform(data(),data()+length(),v2,data(),std::multiplies< T >());
          }

          virtual void operator /=(const basicVector<T>& v2)
          {
            assert( v2.length()==length() ); 

            std::transform(data(),data()+length(),v2.data(),data(),std::divides< T >());
          }

          virtual void operator /=(const T* v2)
          {
            std::transform(data(),data()+length(),v2,data(),std::divides< T >());
          }

          void assign(const utilities::Range& drange, const T data)
          {
            for( int i=0;i<drange.size();++i ) {
              assert( drange[i] < length() );

              v_data[drange[i]] = data;
            }
          }

          void assign(const utilities::Range& drange, const T* data)
          {
            for( int i=0;i<drange.size();++i ) {
              assert( drange[i] < length() );

              v_data[drange[i]] = data[i];
            }
          }

          void assign(const utilities::Range& drange, const basicVector<T>& data)
          {
            for( int i=0;i<drange.size();++i ) {
              assert( drange[i] < length() );

              v_data[drange[i]] = data[i];
            }
          }

          void assign(const utilities::Range& srange, const utilities::Range& drange, const basicVector<T>& data)
          {
            assert( srange.size() == drange.size() );

            for( int i=0;i<srange.size();++i ) {
              assert( srange[i] < data.length() );
              assert( drange[i] < length() );

              v_data[drange[i]] = data[srange[i]];
            }
          } 
        protected:
          T* v_data;
          std::string v_name;
          bool external_data;
          int m_N;
          vectorDotter m_dotter;
#ifdef MNL_MEMORY_VERBOSE_
          int m_memId;
#endif
        private:
          void operator =(const basicVector<T>& input)
          {
          }
      };

    /* real vector class */
    class Vector : public basicVector<Real> {
      public:
        using basicVector<Real>::operator+=;
        using basicVector<Real>::operator-=;
        using basicVector<Real>::operator*=;
        using basicVector<Real>::operator/=;
        using basicVector<Real>::operator[];

        Vector(const std::string& n_name, int size, bool clear=true, Real* n_data=NULL);
        Vector(const Vector& vector);
        Vector(const std::string& n_name, const Vector& vector);
        virtual ~Vector();

        virtual void print() const;

        virtual void load(const std::string& filename);
        virtual void save(const std::string& filename) const;

        void axpy(const Real alpha, const Vector& x);
        void gemv(const Vector& x, const Matrix& A, const char trans=mnlNoTrans, const Real alpha=Real(1), const Real beta=Real(1));

        Real dot(const basicVector<Real>& v2) const;
        Vector invert() const;

        /* operator overloads */
        void operator +=(const Vector&);
        void operator +=(const Real*);

        void operator -=(const Vector&);
        void operator -=(const Real*);

        void operator *=(const Real);

        void operator /=(const Real);

        Vector& operator =(const Vector&);
        Vector& operator =(const Real*);
        Vector& operator =(const Real);

        Vector operator[](const utilities::Range& range) const;
    };

    /* complex vector class */
    class complexVector : public basicVector< std::complex<Real> > {
      public:
        using basicVector< std::complex<Real> >::operator+=;
        using basicVector< std::complex<Real> >::operator-=;
        using basicVector< std::complex<Real> >::operator*=;
        using basicVector< std::complex<Real> >::operator/=;
        using basicVector< std::complex<Real> >::operator[];

        complexVector(const std::string& n_name, int n_length, bool clear=true, std::complex<Real>* n_data=NULL);
        complexVector(const complexVector& vector);
        complexVector(const Vector& vector);
        complexVector(const std::string& n_name, int size, const Real* re, const Real *im);
        complexVector(const std::string& n_name, const complexVector& vector);
        virtual ~complexVector();

        /* read data */
        virtual void print() const;

        void axpy(const Real alpha, const complexVector& x);
        void axpy(const std::complex<Real> alpha, const complexVector& x);
        void axpy(const Real alpha, const Vector& x);
        void gemv(const complexVector& x, const Matrix& A, const char trans=mnlNoTrans, const Real alpha=1.f, const Real beta=1.f);
        void gemv(const complexVector& x, const complexMatrix& A, const char trans=mnlNoTrans, const std::complex<Real> alpha=mnlComplexOne, const std::complex<Real> beta=mnlComplexOne);
        void gemv(const complexMatrix& x, const int row, const complexMatrix& A, const char trans, const std::complex<Real> alpha, const std::complex<Real> beta);
        std::complex<Real> dot(const basicVector<std::complex<Real> >& v2) const;

        Vector real() const;
        Vector imag() const;

        complexVector& operator =(const complexVector&);
        complexVector& operator =(const std::complex<Real>*);
        complexVector& operator =(const Vector&);

        void operator +=(const complexVector&);
        void operator +=(const Vector&);
        void operator +=(const std::complex<Real>*);

        void operator -=(const complexVector&);
        void operator -=(const Vector&);
        void operator -=(const std::complex<Real>*);

        void operator *=(const Vector&);
        void operator *=(const std::complex<Real>);
        void operator *=(const Real);

        void operator /=(const Vector&);
        void operator /=(const Real&);

        complexVector operator[](const utilities::Range& range) const;
    };

    const Vector operator +(const Vector& one, Real two);
    const Vector operator -(const Vector& one, Real two);
    const Vector operator*(const Vector& one, Real two);
    const Vector operator/(const Vector& one, Real two);

    //    template<class T>
    //    const T operator *(Real alpha, const T& one)
    //    {
    //      T result(one);
    //      result *= alpha;

    //      return( result );
    //    }
    //    
    //    template<class T>
    //    const T operator +(const T& one, const T& two)
    //    {
    //      T result(one);
    //      result += two;

    //      return( result );
    //    }
    //    
    //    template<class T>
    //    const T operator -(const T& one, const T& two)
    //    {
    //      T result(one);
    //      result -= two;

    //      return( result );
    //    }

    const Vector loadVector(const std::string& filename);

  } // namespace basics
} // namespace mnl

#endif

