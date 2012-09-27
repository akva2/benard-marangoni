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

#include "grid.h"

namespace mnl {
  namespace utilities {
    basics::Vector equidistant(Real start, Real end, int N)
    {
      basics::Vector result("Equidistant grid",N,false);
      for( int i=0;i<N;++i )
        result[i] = start+(end-start)*static_cast<Real>(i)/static_cast<Real>(N);

      return( result );
    }

    void equidistant(basics::Vector& result, Real start, Real end)
    {
      for( int i=0;i<result.length();++i )
        result[i] = start+(end-start)*static_cast<Real>(i)/static_cast<Real>(result.length());
    }

    basics::Vector wavenumbers(int N)
    {
      basics::Vector result("wavenumbers",N);
      wavenumbers(result);

      return( result );
    }

    void wavenumbers(basics::Vector& result)
    {
      for( int i=0;i<result.length()/2;++i )
        result[i] = (Real)i;
      int j=result.length()/2;
      result[j++] = 0;
      for( int i=-result.length()/2+1;i<0;++i )
        result[j++] = (Real)i;
    }

    void evaluate(basics::Vector& result, MNL_FUNC_1D& f, const basics::Vector& gridX, const Real t)
    {
      assert( result.length() == gridX.length() );

      for( int k=0;k<gridX.length();++k )
        result[k] = f(gridX[k],t);
    }

    void evaluate(basics::Matrix& result, MNL_FUNC_2D& f, const basics::Vector& gridX, const basics::Vector& gridY, const Real t)
    {
      assert( (result.rows() == gridX.length()) && (result.cols() == gridY.length()) );

      for( int j=0;j<gridY.length();++j )
        for( int k=0;k<gridX.length();++k )
          result[j][k] = f(gridX[k],gridY[j],t);
    }

    void evaluate(basics::Matrices& result, MNL_FUNC_3D& f, const basics::Vector& gridX, const basics::Vector& gridY, const basics::Vector& gridZ, const Real t)
    {
      assert( (result.rows() == gridX.length()) && (result.cols() == gridY.length()) && (result.matrices() == gridZ.length()) );

      for( int l=0;l<gridZ.length();++l ) {
        for( int j=0;j<gridY.length();++j )
          for( int k=0;k<gridX.length();++k )
            result[l][j][k] = f(gridX[k],gridY[j],gridZ[l],t);
      }
    }

    void evaluate(basics::Matrix& result, const basics::function2D& f, const basics::Vector& gridX, const basics::Vector& gridY, const Real t)
    {
      assert( (result.rows() == gridX.length()) && (result.cols() == gridY.length()) );

      for( int j=0;j<gridY.length();++j )
        for( int k=0;k<gridX.length();++k )
          result[j][k] = f.val(gridX[k],gridY[j],t);
    }

    void evaluate(basics::Matrices& result, const basics::function3D& f, const basics::Vector& gridX, const basics::Vector& gridY, const basics::Vector& gridZ, const Real t)
    {
      assert( (result.rows() == gridX.length()) && (result.cols() == gridY.length()) && (result.matrices() == gridZ.length()) );

      for( int l=0;l<gridZ.length();++l ) {
        for( int j=0;j<gridY.length();++j )
          for( int k=0;k<gridX.length();++k )
            result[l][j][k] = f.val(gridX[k],gridY[j],gridZ[l],t);
      }
    }

    void evaluate(basics::complexVector& result, MNL_FUNC_1D& f, const basics::Vector& gridX, const Real t)
    {
      assert( result.length() == gridX.length() );

      for( int k=0;k<gridX.length();++k )
        result[k] = f(gridX[k],t);
    }

    void evaluate(basics::complexMatrix& result, MNL_FUNC_2D& f, const basics::Vector& gridX, const basics::Vector& gridY, const Real t)
    {
      assert( (result.rows() == gridX.length()) && (result.cols() == gridY.length()) );

      for( int j=0;j<gridY.length();++j )
        for( int k=0;k<gridX.length();++k )
          result[j][k] = f(gridX[k],gridY[j],t);
    }

    void evaluate(basics::complexMatrices& result, MNL_FUNC_3D& f, const basics::Vector& gridX, const basics::Vector& gridY, const basics::Vector& gridZ, const Real t)
    {
      assert( (result.rows() == gridX.length()) && (result.cols() == gridY.length()) && (result.matrices() == gridZ.length()) );

      for( int l=0;l<gridZ.length();++l ) {
        for( int j=0;j<gridY.length();++j )
          for( int k=0;k<gridX.length();++k )
            result[l][j][k] = f(gridX[k],gridY[j],gridZ[l],t);
      }
    }

    void evaluate(basics::complexMatrix& result, const basics::function2D& f, const basics::Vector& gridX, const basics::Vector& gridY, const Real t)
    {
      assert( (result.rows() == gridX.length()) && (result.cols() == gridY.length()) );

      for( int j=0;j<gridY.length();++j )
        for( int k=0;k<gridX.length();++k )
          result[j][k] = f.val(gridX[k],gridY[j],t);
    }

    void evaluate(basics::complexMatrices& result, const basics::function3D& f, const basics::Vector& gridX, const basics::Vector& gridY, const basics::Vector& gridZ, const Real t)
    {
      assert( (result.rows() == gridX.length()) && (result.cols() == gridY.length()) && (result.matrices() == gridZ.length()) );

      for( int l=0;l<gridZ.length();++l ) {
        for( int j=0;j<gridY.length();++j )
          for( int k=0;k<gridX.length();++k )
            result[l][j][k] = f.val(gridX[k],gridY[j],gridZ[l],t);
      }
    }
  }
}

