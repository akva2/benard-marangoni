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

#include "range.h"

#include <assert.h>
#include <iostream>
#include <memory.h>

namespace mnl {
  namespace utilities {
    void Range::setData(const int size, const int *data)
    {
      for( int i=0;i<size;++i )
        m_index.push_back(data[i]);
    }	

    void Range::setData(const std::vector<int>& data)
    {
      m_index = data;
    }

    int Range::size() const
    {
      return( m_index.size() );
    }

    void Range::size(int n_size)
    {
      m_index.resize(n_size);
    }

    void Range::print() const
    {
      std::cout << "[";
      for( int i=0;i<m_index.size()-1;++i )
        std::cout << m_index[i] << " ";
      std::cout << m_index[m_index.size()-1] << "]" << std::endl;
    }

    Range Range::colon(int start, int end, const int step)
    {
      Range result;
      if( start > end ) { // for fast colon(N,1,-1) as colon(N,1)
        int i=start;
        start = end;
        end = i;
      }
      result.size((end-start)/step+1);
      for( int i=0;i<result.size();++i )
        result[i] = start+i*step;

      return( result );
    } 

    Range2D::Range2D()
    {
      m_rows = -1;
      m_cols = -1;
      m_indices = NULL;
    }

    Range2D::Range2D(int n_rows, int n_cols)
    {
      m_indices = NULL;
      size(n_rows,n_cols);
    }

    Range2D::Range2D(const Range2D& range)
    {
      m_indices = NULL;
      size(range.m_rows,range.m_cols);
      memcpy(m_indices[0],range.m_indices[0],m_rows*m_cols*sizeof(int));
    }

    Range2D::~Range2D()
    {
      deallocate(m_indices);
    }

    void Range2D::size(int n_rows, int n_cols)
    {
      m_rows = n_rows;
      m_cols = n_cols;
      if( m_indices )
        deallocate(m_indices);

      m_indices = allocate(n_rows,n_cols,true);
    }

    int** Range2D::allocate(int n_rows, int n_cols, bool clear)
    {
      int** result = new int*[n_cols];
      result[0] = new int[n_rows*n_cols];
      for( int i=1;i<n_cols;i++ )
        result[i] = (int*)(result[0]+i*n_rows); // manually ensure nice linear memory alignment in fortran order!

      if( clear )
        memset(result[0],0,n_rows*n_cols*sizeof(int));

      return( result );
    }

    void Range2D::deallocate(int** array)
    {
      delete[] array[0];
      delete[] array;
    }

    void Range2D::print()
    {
      for( int i=0;i<m_rows;i++ ) {
        std::cout << "[ ";
        for( int j=0;j<m_cols;j++ )
          std::cout << m_indices[j][i] << " ";
        std::cout << "]" << std::endl;
      }
    }
  }
}

