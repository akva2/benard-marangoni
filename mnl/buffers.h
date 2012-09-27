/***************************************************************************
 *   Copyright (C) 2005-2008 by Arne Morten Kvarving                       *
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
#ifndef MNL_BUFFERS_H_
#define MNL_BUFFERS_H_
#include <vector>
#include <stdio.h>

#include "function.h"

namespace mnl {
  namespace utilities {
    template<class T>
      class ringBuffer {
        public:
          ringBuffer()
          {
          }

          void add(T& input)
          {
            m_entry.push_back(&input);
          }

          T& get(int n, int pos=0)
          {
            if (pos < 0)
              pos += m_entry.size();
            int npos=(n+pos)%m_entry.size();

            return *m_entry[npos];
          }

          const T& get(int n, int pos=0) const
          {
            if (pos < 0)
              pos += m_entry.size();
            int npos=(n+pos)%m_entry.size();

            return *m_entry[npos];
          }

          int size() const
          {
            return m_entry.size();
          }
        protected:
          std::vector<T*> m_entry;
        private:
          ringBuffer(const ringBuffer& foo)
          {
          }
      };

    enum EXTRAPOLANT {
      ZEROTH_ORDER              =  0,
      FIRST_ORDER        =  1,
      SECOND_ORDER       =  2,
      HALF_ORDER         =  3,
      ZEROTH_LSQ         =  4,
      MEAN_FIRST_ZEROTH_LSQ  =  5,
      HALF_LSQ        =  6,
      FIRST_LSQ        =  7,
      SECOND_LSQ        =  8,
      MEAN_FIRST_SECOND    =  9,
      MEAN_FIRST_LSQ_SECOND  = 10,
      THIRD_ORDER        = 11
    };

    template <class T>
      void extrapolate(T& dest, const utilities::ringBuffer<T>& input, int n, EXTRAPOLANT iorder)
      {
        /* setup pressure extrapolant */
        switch( iorder ) {
          case ZEROTH_ORDER:
            /* zeroth order extrapolation */
            dest = input.get(n);
            break;
          case FIRST_ORDER:
            /* first order extrapolation */
            dest = input.get(n,-1);
            dest *= Real(-1);
            dest.axpy(2,input.get(n));
            break;
          default:
          case SECOND_ORDER:
            /* second order extrapolation */
            dest = input.get(n,-2);
            dest.axpy(-3,input.get(n,-1));
            dest.axpy(3,input.get(n));
            break;
          case THIRD_ORDER:
            dest = input.get(n);
            dest *= 4;
            dest.axpy(-6,input.get(n,-1));
            dest.axpy( 4,input.get(n,-2));
            dest.axpy(-1,input.get(n,-3));
            break;
          case HALF_ORDER:
            /* half order extrapolation */
            dest = input.get(n);
            dest *= Real(3)/Real(2);
            dest.axpy(Real(-1)/Real(2),input.get(n,-1));
            break;
          case ZEROTH_LSQ:
            /* zeroth order least squares extrapolation */
            dest = input.get(n);
            dest += input.get(n,-1);
            dest *= Real(1)/Real(2);
            break;
          case MEAN_FIRST_ZEROTH_LSQ:
            /* mean of zeroth order least squares extrapolation and first order extrapolation */
            dest = input.get(n);
            dest *= Real(5)/Real(4);
            dest.axpy(Real(-1)/Real(4),input.get(n,-1));
            break;
          case HALF_LSQ:
            /* half order least squares extrapolation */
            dest = input.get(n,-2);
            dest *= Real(-1)/Real(3);
            dest.axpy(Real(5)/Real(12),input.get(n,-1));
            dest.axpy(Real(11)/Real(12),input.get(n));
            break;
          case FIRST_LSQ:
            /* first order least-squares extrapolation */
            dest = input.get(n,-2);
            dest *= Real(-2)/Real(3);
            dest.axpy(Real(1)/Real(3),input.get(n,-1));
            dest.axpy(Real(4)/Real(3),input.get(n));
            break;
          case SECOND_LSQ:
            dest = input.get(n,-3);
            dest *= Real(3)/Real(4);
            dest.axpy(Real(-5)/Real(4),input.get(n,-2));
            dest.axpy(Real(-3)/Real(4),input.get(n,-1));
            dest.axpy(Real( 9)/Real(4),input.get(n));
            break;
          case MEAN_FIRST_SECOND:
            dest = input.get(n,-2);
            dest *= Real(1)/2;
            dest.axpy(-2,input.get(n,-1));
            dest.axpy(Real(5)/2,input.get(n));
            break;
          case MEAN_FIRST_LSQ_SECOND:
            dest = input.get(n,-2);
            dest *= Real(1)/6;
            dest.axpy(-Real(4)/3,input.get(n,-1));
            dest.axpy(Real(13)/6,input.get(n));
            break;
        }
      }
  }
  namespace basics {
    template<class T>
      class matrixStackT {
        public:
          class stackDotter {
            public:
              virtual Real operator()(const matrixStackT<T>& A, const matrixStackT<T>& B) const
              {
                return A.dot(B);
              }
          };

          matrixStackT(const std::string& name) :
            m_name(name)
        {
          m_collection.clear();
          m_mine = false;
        }

          matrixStackT(const std::string& name, int rows,
              int cols, int size) :
            m_name(name)
        {
          T templ(name,rows,cols);
          add(templ,size);
        }

          matrixStackT(const std::string& name, std::vector<T*> collection) :
            m_name(name)
        {
          m_collection = collection;
          m_mine = false;
        }

          virtual ~matrixStackT()
          {
            if ( m_mine )
              for( int i=0;i<m_collection.size();++i )
                delete m_collection[i];

            m_collection.clear();
          }

          matrixStackT(const matrixStackT& in)
          {
            if( in.m_mine )
              add(*in.m_collection[0],in.size());
            else
              m_collection = in.m_collection;

            m_name = in.m_name;
            m_mine = in.m_mine;
          }

          const std::string& name() const
          {
            return m_name;
          }

          void setName(const std::string& name)
          {
            m_name = name;
          }

          void add(const T& input, int n)
          {
            for( int i=0;i<n;++i ) {
              T* foo = new T(input);
              assert( foo );
              char temp[1024];
              sprintf(temp,"%s (part %zu/%i)",m_name.c_str(),m_collection.size()+1,n);
              foo->setName(temp);
              foo->clear();
              m_collection.push_back(foo);
            }
            m_mine = true;
          }

          void add(std::vector<T*>& input)
          {
            m_collection = input;
            m_mine = false;
          }

          int size() const
          {
            return( m_collection.size() );
          }

          int length() const
          {
            return( m_collection[0]->length()*m_collection.size() );
          }

          Real sum() const
          {
            Real result=0;
            int max=size();
#pragma omp parallel for reduction (+:result) schedule(static)
            for( int i=0;i<max;++i )
              result += m_collection[i]->sum();

            return( result );
          }

          Real max() const
          {
            Real result=0;
            for( int i=0;i<size();++i )
              if( m_collection[i]->max() > result )
                result = m_collection[i]->max();

            return( result );
          }

          void clear()
          {
            int max=size();
#pragma omp parallel for schedule(static)
            for( int i=0;i<max;++i )
              m_collection[i]->clear();
          }

          void fill(const Real& in)
          {
            int max=size();
#pragma omp parallel for schedule(static)
            for( int i=0;i<max;++i )
              m_collection[i]->fill(in);
          }

          const stackDotter& getDotter() const
          {
            return m_dotter;
          }

          Real dot(const matrixStackT<T>& B) const
          {
            Real result=0;
            int max=size();
#pragma omp parallel for reduction (+:result) schedule(static)
            for( int i=0;i<max;++i )
              result += (*this)[i].dot(B[i]);

            return( result );
          }

          inline T& operator[](unsigned int index)
          {
            assert( index < m_collection.size() );

            return *m_collection[index];
          }

          inline const T& operator[](unsigned int index) const
          {
            assert( index < m_collection.size() );

            return *m_collection[index];
          }

          inline void operator -=(Real alpha)
          {
            int max=size();
#pragma omp parallel for schedule(static)
            for( int i=0;i<max;++i )
              (*m_collection[i]) -= alpha;
          }

          inline matrixStackT& operator=(const matrixStackT& in)
          {
            assert( size() == in.size() );
            int max=size();
#pragma omp parallel for schedule(static)
            for( int i=0;i<max;++i )
              (*m_collection[i]) = in[i];

            return( *this );
          }

          inline matrixStackT& operator=(const Real in)
          {
            int max=size();
#pragma omp parallel for schedule(static)
            for( int i=0;i<max;++i )
              (*m_collection[i]) = in;

            return( *this );
          }

          inline matrixStackT& operator+=(const Real in)
          {
            int max=size();
#pragma omp parallel for schedule(static)
            for( int i=0;i<max;++i )
              (*m_collection[i]) += in;

            return( *this );
          }

          inline void operator +=(const matrixStackT& in)
          {
            int max=size();
#pragma omp parallel for schedule(static)
            for( int i=0;i<max;++i )
              (*m_collection[i]) += in[i];
          }

          inline void operator -=(const matrixStackT& in)
          {
            int max=size();
#pragma omp parallel for schedule(static)
            for( int i=0;i<max;++i )
              (*m_collection[i]) -= in[i];
          }

          inline void operator *=(const Real in)
          {
            int max=size();
#pragma omp parallel for schedule(static)
            for( int i=0;i<max;++i )
              (*m_collection[i]) *= in;
          }

          inline void axpy(Real alpha, const matrixStackT& in)
          {
            int max=size();
#pragma omp parallel for schedule(static)
            for( int i=0;i<max;++i )
              m_collection[i]->axpy(alpha,in[i]);
          }

          void print(int precision=5) const
          {
            for( int i=0;i<size();++i )
              m_collection[i]->print(precision);
          }

          void info() const
          {
            std::cout << name() << " " << size() << " elements, ";
            m_collection[0]->info();
          }

          void save(const std::string& name, const std::vector<int>& r) const
          {
            for( int i=0;i<size();++i ) {
              std::stringstream str;
              str << name;
              str << "-";
              str << r[i];
              m_collection[i]->save(str.str());
            }
          }

          void save(const std::string& name) const
          {
            std::vector<int> r;
            for( int i=0;i<size();++i )
              r.push_back(i);
            save(name,r);
          }
        protected:
          std::string m_name;
          std::vector<T*> m_collection;
          stackDotter m_dotter;
          bool m_mine;
      };

    typedef matrixStackT<basics::Matrix> matrixStack;

    class matricesStack : public matrixStackT<basics::Matrices> {
      public:
        using matrixStackT<basics::Matrices>::operator=;

        matricesStack(const std::string& name) :
          matrixStackT<basics::Matrices>(name)
      {
      }

        matricesStack(const std::string& name, int rows,
                      int cols, int matrices, int size) :
          matrixStackT<basics::Matrices>(name)
      {
        Matrices templ("template",rows,cols,matrices,size);
        matricesStack::add(templ,size);
      }

        matricesStack(const std::string& name,
                      std::vector<basics::Matrices*> collection) :
          matrixStackT<basics::Matrices>(name,collection)
      {
        /* setup the 'fake' matrix stack */
        for( int l=0;l<(*m_collection[0]).matrices();++l ) {
          std::vector<Matrix*> vec;
          for( int i=0;i<size();++i )
            vec.push_back(&((*m_collection[i])[l]));
          matrixStack mat(name+" (fake matrix stack)");
          mat.add(vec);
          m_stacks.push_back(mat);
        }
      }

        virtual ~matricesStack()
        {
        }

        void add(const basics::Matrices& input, int n)
        {
          matrixStackT<Matrices>::add(input,n);
          /* setup the 'fake' matrix stack */
          for( int l=0;l<input.matrices();++l ) {
            std::vector<Matrix*> vec;
            for( int i=0;i<n;++i )
              vec.push_back(&((*m_collection[i])[l]));
            matrixStack mat(name()+" (fake matrix stack)");
            mat.add(vec);
            m_stacks.push_back(mat);
          }
        }

        matrixStack& at(int l)
        {
          return m_stacks[l];
        }

        const matrixStack& at(int l) const
        {
          return m_stacks[l];
        }
        std::vector<matrixStack> m_stacks;
    };

    template<class T, class T2 >
      class componentView {
        public:
          componentView(const std::vector<T2*>& stack)
          {
            for( int j=0;j<stack[0]->size();++j ) {
              std::vector<T*> vec;
              for( int i=0;i<stack.size();++i )
                vec.push_back(&(*stack[i])[j]);
              m_stacks.push_back(new T2("view",vec));
            }
          }

          T2& operator[](const int index)
          {
            return *m_stacks[index];
          }

          ~componentView()
          {
            for( int i=0;i<m_stacks.size();++i )
              delete m_stacks[i];
          }
        protected:
          std::vector<T2*> m_stacks;
      };
  }
}

#endif

