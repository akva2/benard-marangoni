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
#ifndef MNL_MATRICES_H_
#define MNL_MATRICES_H_

#include "config.h"
#include "matrix.h"

#include <vector>
#include <string>
#include <cassert>

#include <sstream>

namespace mnl {
  namespace basics {
    template<class T, class M>
      class basicMatrices {
        public:
          class matricesDotter {
            public:
              Real operator()(const basicMatrices<T,M>& A,
                  const basicMatrices<T,M>& B) const
              {
                return A.dot(B);
              }
          };

          basicMatrices(const std::string& n_name, int n_rows, int n_cols, int n_k, bool clear, T* n_data=NULL) :
            m_name(n_name), 
            m_k(n_k), 
            m_rows(n_rows), 
            m_cols(n_cols)
        {
          external_data = n_data?true:false;
          if( !n_data ) {
            ALLOCATELOG(name(),length());
            the_data = new T[length()]; // containing the ACTUAL data!
#ifdef MNL_MEMORY_VERBOSE_
            m_memId = utilities::g_tracker.addChunk(utilities::memTracker::Matrices,length()*sizeof(T));
#endif
            if( clear )
              memset(the_data,0,length()*sizeof(T));
          } else
            the_data = n_data;

          for( int i=0;i<matrices();++i ) {
            m_data.push_back((T**)malloc(cols()*sizeof(T*)));
            for( int j=0;j<cols();++j )
              m_data[i][j] = the_data+i*(rows()*cols())+j*rows();
          }
          setupMatrices();
        }

          basicMatrices(const basicMatrices<T,M>& matrices) :
            m_name(matrices.name()),
            m_rows(matrices.rows()),
            m_cols(matrices.cols()),
            m_k(matrices.matrices())
        {
          ALLOCATELOG(name()+" (copy constructor)",length());
          the_data = new T[length()]; // containing the ACTUAL data!
#ifdef MNL_MEMORY_VERBOSE_
          m_memId = utilities::g_tracker.addChunk(utilities::memTracker::Matrices,length()*sizeof(T));
#endif
          for( int i=0;i<m_k;++i ) {
            m_data.push_back((T**)malloc(m_cols*sizeof(T*)));
            for( int j=0;j<m_cols;++j )
              m_data[i][j] = the_data+i*(m_rows*m_cols)+j*m_rows;
          }
          external_data = false;
          setupMatrices();
        }

          virtual ~basicMatrices()
          {
            DEALLOCATELOG(name());
            if( !external_data ) {
              delete[] the_data;
#ifdef MNL_MEMORY_VERBOSE_
              utilities::g_tracker.removeChunk(m_memId);
#endif
            }
            destroyMatrices();
            for( int i=0;i<matrices();++i )
              free(m_data[i]);
          }

          inline const int& matrices() const
          {
            return( m_k );
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
            return rows()*cols()*matrices();
          }

          const T* data() const
          {
            return( the_data );
          }

          T* data()
          {
            return( the_data );
          }

          T max() const
          {
            return *std::max_element(the_data,the_data+length());
          }

          void clear()
          {
            assert( the_data );

            memset(the_data,0,length()*sizeof(T));
          }

          void fill(const T& in)
          {
            assert( the_data );
            for( int i=0;i<length();++i )
              the_data[i] = in;
          }

          void save(const std::string& name, const Matrix::fileFormat format=Matrix::ASCII) const
          {
            if( format == Matrix::ASCII ) {
              for( int i=0;i<matrices();++i ) {
                std::stringstream s;
                s << name << "-" << i << ".asc";
                m_matrix[i]->save(s.str());
              }
            }
          }

          virtual void info() const
          {
            std::cout << m_name << ": " << rows() << "x" << cols() << "x" << matrices() << std::endl;
          }

          Real sum() const
          {
            Real result=0;
            int max=length();
#pragma omp parallel for schedule(static) reduction(+:result)
            for( int i=0;i<max;++i ) 
              result += the_data[i];

            return( result );
          }

          virtual T dot(const basicMatrices<T,M>& in) const = 0;
          const matricesDotter& getDotter() const
          {
            return m_dotter;
          }

          const std::string& name() const
          {
            return( m_name );
          }

          void setName(const std::string& name)
          {
            m_name = name;
          }

          void print(int precision=5) const
          {
            for( int i=0;i<m_k;++i )
              m_matrix[i]->print(precision);
          }

          inline const M& operator[](int index) const
          {
            assert( index < m_k );
            return( *m_matrix[index] );
          }

          inline M& operator[](const int index)
          {
            assert( index < m_k );
            return( *m_matrix[index] );
          }
        protected:
          std::string m_name;
          int m_k, m_rows, m_cols;
          std::vector<T**> m_data;
          T* the_data;
          std::vector<M*> m_matrix;
          bool external_data;
          matricesDotter m_dotter;
#ifdef MNL_MEMORY_VERBOSE_
          int m_memId;
#endif
          void setupMatrices()
          {
            destroyMatrices();
            for( int i=0;i<matrices();++i ) {
              std::stringstream s;
              s << "Matrix " << i << " from " << name();
              m_matrix.push_back(new M(s.str(),rows(),cols(),false,m_data[i]));
            }
          }

          void destroyMatrices()
          {
            for( int i=0;i<m_matrix.size();++i )
              delete m_matrix[i];

            m_matrix.clear();
          }
        private:
          void operator =(const basicMatrices& matrices)
          {
          }
      };

    class Matrices : public basicMatrices<Real,Matrix> {
      public:
        Matrices(const std::string& n_name, int n_rows, int n_cols, int n_k, bool clear=true, Real* n_data=NULL);
        Matrices(const Matrices&);
        Matrices(const std::string& n_name, const Matrices&);
        virtual ~Matrices();

        Matrices& operator =(const Matrices&);
        Matrices& operator =(const Real);

        void operator +=(const Matrices&);
        void operator +=(const Real);
        void operator -=(const Matrices&);
        void operator -=(const Real);

        void operator *=(const Real);

        void axpy(const Real alpha, const Matrices& x);
        Real dot(const basicMatrices<Real,Matrix>& x) const;

        inline const Matrix& operator[](int index) const
        {
          assert( index < m_k );
          return( *m_matrix[index] );
        }	
        inline Matrix& operator[](const int index)
        {
          assert( index < m_k );
          return( *m_matrix[index] );
        }
    };

    class complexMatrices : public basicMatrices<std::complex<Real>,complexMatrix> {
      public:
        complexMatrices(const std::string& n_name, int n_rows, int n_cols, int n_k, bool clear=true);
        complexMatrices(const complexMatrices&);
        complexMatrices(const std::string& name, const complexMatrices&);
        virtual ~complexMatrices();

        Matrices real() const;
        Matrices imag() const;

        void axpy(const std::complex<Real> alpha, const complexMatrices& x);
        void axpy(const Real alpha, const complexMatrices& x);
        std::complex<Real> dot(const basicMatrices<std::complex<Real>,complexMatrix>& x) const;

        complexMatrices& operator =(const complexMatrices&);
        complexMatrices& operator =(const Matrices&);

        void operator -=(const complexMatrices&);
        void operator -=(const Matrices&);
        void operator -=(const Real mean);

        void operator +=(const complexMatrices&);
        void operator += (const Matrices&);

        void operator *=(const std::complex<Real>);
        void operator *=(const Real);

        complexMatrices operator -(const complexMatrices&) const;
    };

    void multPointwise(complexMatrices& result, const complexMatrices& phi, 
        const complexMatrices& u);
    void multPointwise(Matrices& result, const Matrices& phi, const Matrices& u);
    void multPointwise(Matrices& result, const Matrices& phi);
    void multPointwiseReal(complexMatrices& result, const complexMatrices& phi,
        const complexMatrices& u);
    void multPointwise(complexMatrices& phi, const complexMatrices& u);
    void multPointwiseReal(complexMatrices& phi, const complexMatrices& u);
    void applyLocalGlobal(basics::Matrices& result, const basics::Matrices& victim, 
        const basics::Matrix& op, char transA, char transB, int ofs=0, Real scaleA=Real(1), Real scaleB=0);


  } // namespace basics
} // namespace mnl

#endif

