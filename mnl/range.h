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
#ifndef MNL_RANGE_H_
#define MNL_RANGE_H_

#include <vector>
#include <assert.h>

namespace mnl {
  namespace utilities {
    /** main indices class **/
    class Range {
      public:
        Range()
        {
        }

        Range(int size)
        {
          m_index.reserve(size);
        }
        ~Range()
        {
        }

        /* assignment */
        void setData(const int size, const int *data); // only copies data !
        void setData(const std::vector<int>& index);

        int size() const;
        void size(int n_size);
        void print() const;

        /* generators  */
        static Range colon(int start, int end, const int step=1);

        /* operators */
        inline int operator[](int i) const
        {
          assert( i < m_index.size() );
          return( m_index[i] );
        } 

        inline int& operator[](int i)
        {
          if (i >= size() )
            m_index.resize(i+1);

          return( m_index[i] );
        } 
      protected:
        std::vector<int> m_index;
    };

    class Range2D {
      public:
        Range2D();
        Range2D(const Range2D& range);
        Range2D(int n_rows, int n_cols);
        ~Range2D();

        void size(int n_rows, int n_cols);

        inline const int* operator[](int col) const
        {
          assert( m_indices && col < m_cols );

          return m_indices[col];
        }

        inline int* operator[](int col)
        {
          assert( m_indices && col < m_cols );

          return m_indices[col];
        }

        inline int cols() const
        {
          return m_cols;
        }

        inline int rows() const
        {
          return m_rows;
        }

        void print();
      protected:
        static int** allocate(int n_rows, int n_cols, bool clear);
        static void deallocate(int**);
        int** m_indices; 
        int m_cols;
        int m_rows;
    };
  }
}

#endif

