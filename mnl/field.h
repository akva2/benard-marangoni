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
 *   but WITHOUT ANY WARRANTm_Y; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef MNL_FIELD_H_
#define MNL_FIELD_H_

#include "config.h"
#include "matrix.h"
#include "matrices.h"

namespace mnl {
  namespace basics {
    template<class T>
      class Field2 {
        public:
          Field2(std::string n_name, int n_rows, int n_cols, bool clear=true)
          {
            m_X = new T(n_name+"-X",n_rows,n_cols,clear);
            m_Y = new T(n_name+"-Y",n_rows,n_cols,clear);
            mine = true;
          }

          Field2(T* n_X, T* n_Y)
          {
            m_X = n_X;
            m_Y = n_Y;
            mine = false;
          }

          ~Field2()
          {
            if( mine ) {
              delete m_X;
              delete m_Y;
            }
          }

          void print(int precision=5) const
          {
            m_X->print(precision);
            m_Y->print(precision);
          }

          void save(const std::string& name, const Matrix::fileFormat format=Matrix::ASCII) const
          {
            m_X->save(name+"-X.asc",format);
            m_Y->save(name+"-Y.asc",format);
          }

          void setName(const std::string& name)
          {
            m_X->setName(name+"-X");
            m_Y->setName(name+"-Y");
          }

          inline Field2<T>& operator=(const Field2<T>& rhs)
          {
            *m_X = rhs.X();
            *m_Y = rhs.Y();

            return *this;
          }

          inline void operator *=(const Real scale)
          {
            (*m_X) *= scale;
            (*m_Y) *= scale;
          }

          inline void operator +=(const Field2<T>& rhs)
          {
            X() += rhs.X();
            Y() += rhs.Y();
          }

          inline void operator -=(const Field2<T>& rhs)
          {
            X() -= rhs.X();
            Y() -= rhs.Y();
          }

          inline void axpy(const Real alpha, const Field2<T>& rhs)
          {
            X().axpy(alpha,rhs.X());
            Y().axpy(alpha,rhs.Y());
          }

          inline void clear()
          {
            X().clear();
            Y().clear();
          }

          inline T& X()
          {
            return *m_X;
          }

          inline const T& X() const
          {
            return *m_X;
          }

          inline T& Y()
          {
            return *m_Y;
          }

          inline const T& Y() const
          {
            return *m_Y;
          }

          inline T& operator[](const int index)
          {
            if( index == 0 )
              return X();
            if( index == 1 )
              return Y();
            assert(0);
          }

          inline const T& operator[](const int index) const
          {
            if( index == 0 )
              return X();
            if( index == 1 )
              return Y();
            assert(0);
          }

          inline const int& rows() const
          {
            return( m_X->rows() );
          }

          inline const int& cols() const
          {
            return( m_X->cols() );
          }
        protected:
          T* m_X;
          T* m_Y;
          bool mine;
      };

    typedef Field2<Matrix> Field2D;
    typedef Field2<complexMatrix> complexField2D;

    template<class T>
      class Field3 {
        public:
          Field3(std::string n_name, int n_rows, int n_cols, int n_matrices, bool clear=true)
          {
            m_X = new T(n_name+"-X",n_rows,n_cols,n_matrices,clear);
            m_Y = new T(n_name+"-Y",n_rows,n_cols,n_matrices,clear);
            m_Z = new T(n_name+"-Z",n_rows,n_cols,n_matrices,clear);
            mine = true;
          }

          Field3(T* n_X, T* n_Y, T* n_Z)
          {
            m_X = n_X;
            m_Y = n_Y;
            m_Z = n_Z;
            mine = false;
          }

          ~Field3()
          {
            if( mine ) {
              delete m_X;
              delete m_Y;
              delete m_Z;
            }
          }

          void print(int precision=5) const
          {
            m_X->print(precision);
            m_Y->print(precision);
            m_Z->print(precision);
          }

          void save(const std::string name, const Matrix::fileFormat format=Matrix::ASCII) const
          {
            m_X->save(name+"-X");
            m_Y->save(name+"-Y");
            m_Z->save(name+"-Z");
          }

          void setName(const std::string& name)
          {
            m_X->setName(name+"-X");
            m_Y->setName(name+"-Y");
            m_Z->setName(name+"-Z");
          }

          inline Field3<T>& operator=(const Field3<T>& rhs)
          {
            *m_X = rhs.X();
            *m_Y = rhs.Y();
            *m_Z = rhs.Z();

            return *this;
          }

          inline Field3<T>& operator=(const Real alpha)
          {
            *m_X = alpha;
            *m_Y = alpha;
            *m_Z = alpha;

            return *this;
          }

          inline void operator *=(const Real scale)
          {
            (*m_X) *= scale;
            (*m_Y) *= scale;
            (*m_Z) *= scale;
          }

          inline void operator +=(const Field3<T>& rhs)
          {
            X() += rhs.X();
            Y() += rhs.Y();
            Z() += rhs.Z();
          }

          inline void operator -=(const Field3<T>& rhs)
          {
            X() -= rhs.X();
            Y() -= rhs.Y();
            Z() -= rhs.Z();
          }

          inline void axpy(const Real alpha, const Field3<T>& rhs)
          {
            X().axpy(alpha,rhs.X());
            Y().axpy(alpha,rhs.Y());
            Z().axpy(alpha,rhs.Z());
          }

          inline void clear()
          {
            X().clear();
            Y().clear();
            Z().clear();
          }

          inline T& X()
          {
            return *m_X;
          }

          inline const T& X() const
          {
            return *m_X;
          }

          inline T& Y()
          {
            return *m_Y;
          }

          inline const T& Y() const
          {
            return *m_Y;
          }

          inline T& Z()
          {
            return *m_Z;
          }

          inline const T& Z() const
          {
            return *m_Z;
          }

          inline T& operator[](const int index)
          {
            if( index == 0 )
              return X();
            if( index == 1 )
              return Y();
            if( index == 2 )
              return Z();
            assert(0);
            return X(); // to shut the compiler up
          }

          inline const T& operator[](const int index) const
          {
            if( index == 0 )
              return X();
            if( index == 1 )
              return Y();
            if( index == 2 )
              return Z();
            assert(0);
            return X(); // to shut the compiler up
          }

          inline const int& matrices() const
          {
            return( m_X->matrices() );
          }

          inline const int& rows() const
          {
            return( m_X->rows() );
          }

          inline const int& cols() const
          {
            return( m_X->cols() );
          }
        protected:
          T* m_X;
          T* m_Y;
          T* m_Z;
          bool mine;
      };

    typedef Field3<Matrices> Field3D;
    typedef Field3<complexMatrices> complexField3D;
  } // namespace basics
} // namespace mnl

#endif

