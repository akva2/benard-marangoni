/***************************************************************************
 *   Copyright (C) 2006 by Arne Morten Kvarving                            *
 *   arnemort@fourier04
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
#ifndef MNL_MAPPING_H_
#define MNL_MAPPING_H_

/** -- muu numerics library --
* - this is the mapping class, for use with element methods 
*/
#include "config.h"
#include "range.h"
#include "vector.h"
#include "matrix.h"

#include <map>
#include <vector>
#include <assert.h>

namespace mnl {
  namespace utilities {
    class Element {
      public:
        Element(int N, Real x0, Real size, const basics::Vector& grid, const basics::Vector& weight) : m_grid(grid), m_weight(weight), m_op("elemental operator",N,N), m_rhs("elemental rhs",N)
      {   
        m_x0 = x0;
        m_size = size;
      }

        virtual ~Element() {};

        inline basics::Matrix& getOperator()
        {
          return m_op;
        }	

        inline basics::Vector& getRHS()
        {
          return m_rhs;
        }

        inline basics::Vector& getGrid()
        {
          return m_grid;
        }

        inline basics::Vector& getWeight()
        {
          return m_weight;
        }

        inline Real size()
        {
          return m_size;
        }

        inline Real x0()
        {
          return m_x0;
        }
      protected:
        Real m_x0;
        Real m_size;
        basics::Vector m_grid;
        basics::Vector m_weight;
        basics::Matrix m_op;
        basics::Vector m_rhs;
    };

    class Element2D {
      public:
        Element2D(int N, int M, Real x0, Real y0, Real width, Real height, const basics::Vector& gridx, const basics::Vector& gridy, const basics::Vector& weightx, const basics::Vector& weighty) : m_gridx(gridx), m_gridy(gridy), m_weightx(weightx), m_weighty(weighty), m_op("elemental operator",N*M,N*M), m_rhs("elemental rhs",N*M), m_1Dop("dummy",1,1)
      {   
        m_x0 = x0;
        m_y0 = y0;
        m_width = width;
        m_height = height;
        m_N = N;
        m_M = M;
      }

        Element2D(int N, int M, int K, Real x0, Real y0, Real width, Real height, const basics::Vector& gridx, const basics::Vector& gridy, const basics::Vector& weightx, const basics::Vector& weighty) : m_gridx(gridx), m_gridy(gridy), m_weightx(weightx), m_weighty(weighty), m_op("elemental operator",N*M,N*M), m_rhs("elemental rhs",N*M), m_1Dop("1d operator",K,K) // K = dim of 1D op
      {
        m_x0 = x0;
        m_y0 = y0;
        m_width = width;
        m_height = height;
        m_N = N;
        m_M = M;
      }

        virtual ~Element2D() {};

        inline basics::Matrix& getOperator()
        {
          return m_op;
        }	

        inline basics::Vector& getRHS()
        {
          return m_rhs;
        }

        inline basics::Vector& getGridX()
        {
          return m_gridx;
        }

        inline basics::Vector& getGridY()
        {
          return m_gridy;
        }

        inline basics::Vector& getWeightX()
        {
          return m_weightx;
        }

        inline basics::Vector& getWeightY()
        {
          return m_weighty;
        }

        inline Real width()
        {
          return m_width;
        }

        inline Real height()
        {
          return m_height;
        }

        inline Real x0()
        {
          return m_x0;
        }

        inline Real y0()
        {
          return m_y0;
        }

        inline int N()
        {
          return m_N;
        }

        inline int M()
        {
          return m_M;
        }

        virtual void assembleOperator(const Range& nodes, const basics::Matrix& op)
        {
        }

        inline void set1Dop(basics::Matrix& mat)
        {
          m_1Dop = mat;
        }

        inline basics::Matrix& get1Dop()
        {
          return m_1Dop;
        }
      protected:
        Real m_x0;
        Real m_y0;
        Real m_width;
        Real m_height;
        int m_N;
        int m_M;
        basics::Vector m_gridx;
        basics::Vector m_gridy;
        basics::Vector m_weightx;
        basics::Vector m_weighty;
        basics::Matrix m_op;
        basics::Matrix m_1Dop;
        basics::Vector m_rhs;
    };

    class Mapping {
      public:
        Mapping();
        ~Mapping();

        inline const Range& operator[](int index) 
        {
          std::map<int, std::pair<Range,Element> >::iterator iter=m_mapping.find(index);

          assert( iter != m_mapping.end() );

          return( iter->second.first );
        }

        void addElement(int index, const Range& nodes, const Element& element);

        inline Element& getElement(int index) 
        {
          std::map<int, std::pair<Range,Element> >::iterator iter=m_mapping.find(index);

          assert( iter != m_mapping.end() );

          return( iter->second.second );
        }

        void assemble(basics::Matrix& op);
        void assembleRHS(basics::Vector& vec);
        void assembleGrid(basics::Vector& grid);

        void print();

      protected:
        std::map<int,std::pair<Range,Element> > m_mapping;
    };

    class Mapping2D {
      public:
        Mapping2D();
        ~Mapping2D();

        inline const Range2D& operator[](int index) 
        {
          std::map<int, std::pair<Range2D,Element2D> >::iterator iter=m_mapping.find(index);

          assert( iter != m_mapping.end() );

          return( iter->second.first );
        }

        inline Element2D& getElement(int index) 
        {
          std::map<int, std::pair<Range2D,Element2D> >::iterator iter=m_mapping.find(index);

          assert( iter != m_mapping.end() );

          return( iter->second.second );
        }

        void addElement(const int, const Range2D& nodes, const Element2D& element);

        void assemble(basics::Matrix& op);
        void assembleRHS(basics::Vector& vec);
        void assembleGrid(basics::Vector& grid);

        void print();

      protected:
        std::map<int,std::pair<Range2D,Element2D> > m_mapping;
    };
  }
}

#endif

