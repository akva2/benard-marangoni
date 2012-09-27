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

#include "gll.h"

#include <vector>
#include <algorithm>
#include <assert.h>

namespace mnl {
  namespace utilities {
    Real GLL::Legendre(const Real x, const int N)
    {
      std::vector<Real> Ln;
      Ln.resize(N+1);
      Ln[0] = 1.f;
      Ln[1] = x;
      if( N > 1 )
        for( int i=1;i<N;i++ )
          Ln[i+1] = (2*(Real)i+1)/((Real)i+1)*x*Ln[i]-i/((Real)i+1)*Ln[i-1];

      return( Ln[N] );
    }

    Real GLL::LegendreDerivative(const Real x, const int N)
    {
      std::vector<Real> Ln;
      Ln.resize(N+1);

      Ln[0] = 1.f; Ln[1] = x;

      if( (x == 1) || (x == -1) )
        return( pow(x,N-1)*(Real)N*((Real)N+1.f)/2.f );
      else {
        for( int i=1;i<N;i++ )
          Ln[i+1] = (2.f*(Real)i+1.f)/((Real)i+1.f)*x*Ln[i]-(Real)i/((Real)i+1)*Ln[i-1];
        return( (Real)N/(1.f-x*x)*Ln[N-1]-(Real)N*x/(1-x*x)*Ln[N] );
      }
    }

    Real GLL::Lagrange(const Real x, const int nr, const basics::Vector& grid)
    {
      assert( nr < grid.length() );

      Real result = 1.f;
      for( int j=0;j<grid.length();++j )
        if( j != nr ) 
          result *= (x-grid[j])/(grid[nr]-grid[j]);

      return( result );
    }

    void GLL::GaussLegendreGrid(basics::Vector& grid)
    {
      basics::Matrix A("GaussLegendre matrix",grid.length(),grid.length(),true);
      A[0][1] = 1.f;
      for( int i=1;i<grid.length()-1;i++ ) {
        A[i][i-1] = (Real)i/(2.f*(Real)(i+1.f)-1.f);
        A[i][i+1] = (Real)(i+1.f)/(2.f*(Real)(i+1.f)-1.f);
      }
      A[grid.length()-1][grid.length()-2] = ((Real)grid.length()-1.f)/(2*(Real)grid.length()-1);
      grid = A.eigenValues().real();
      std::sort(grid.data(),grid.data()+grid.length());
    }

    void GLL::GaussLegendreWeights(basics::Vector& weight, const basics::Vector& grid)
    {
      for( int j=0;j<grid.length();j++ ) {
        Real p = 1.f-grid[j]*grid[j];
        Real p2 = LegendreDerivative(grid[j],grid.length());
        weight[j] = 2.f/(p*p2*p2);
      }
    }

    void GLL::GaussLobattoLegendreGrid(basics::Vector& grid)
    {
#define TOLERANCE 1.e-15f

      grid[0] = -1.f;
      grid[grid.length()-1] = 1.f; // Lobatto

      if( grid.length() == 2 )
        return;

      basics::Vector GLgrid("Gauss Legendre Grid",grid.length()-1,false);
      GaussLegendreGrid(GLgrid);

      for( int i=1;i<grid.length()-1;i++ ) {
        grid[i] = (GLgrid[i-1]+GLgrid[i])/2.f;
        Real old = 0.f;
        while( fabs(old-grid[i]) > TOLERANCE ) {
          old = grid[i];
          Real L = Legendre(old,grid.length()-1);
          Real Ld = LegendreDerivative(old,grid.length()-1);
          grid[i] += (1.f-old*old)*Ld/(((Real)grid.length()-1.f)*(Real)grid.length()*L);
        }
      }
    }

    void GLL::GaussLobattoLegendreGridPeriodic(basics::Vector& grid)
    {
      basics::Vector tempGrid("temp",grid.length()+1);
      GaussLobattoLegendreGrid(tempGrid);
      grid = tempGrid[utilities::Range::colon(0,grid.length()-1)];
    }

    void GLL::GaussLobattoLegendreGridHomogenous(basics::Vector& grid)
    {
      basics::Vector tempGrid("temp",grid.length()+2);
      GaussLobattoLegendreGrid(tempGrid);
      grid = tempGrid[utilities::Range::colon(1,grid.length())];
    }

    void GLL::GaussLobattoLegendreWeights(basics::Vector& weight, const basics::Vector& grid)
    {
      assert( grid.length() == weight.length() );

      int N = grid.length();
      for( int i=0;i<N;++i ) {
        Real L=Legendre(grid[i],N-1);
        weight[i] = 2.f/(((Real)N-1)*N*L*L);
      }
    }

    void GLL::GaussLobattoLegendreWeightsPeriodic(basics::Vector& weight, const basics::Vector& grid)
    {
      basics::Vector gridTemp("temp",grid.length()+1);
      BLASRPFX(copy,BLASINT(grid.length()),BLASPREAL(grid.data()),BLASINT(mnlIntOne),BLASPREAL(gridTemp.data()),BLASINT(mnlIntOne));
      gridTemp[grid.length()] = 1;

      GaussLobattoLegendreWeights(weight,gridTemp);
    }

    void GLL::GaussLobattoLegendreWeightsHomogenous(basics::Vector& weight, const basics::Vector& grid)
    {
      basics::Vector temp("temp",weight.length()+2);
      basics::Vector gridTemp("temp",grid.length()+2);
      gridTemp[0] = -1;
      BLASRPFX(copy,BLASINT(grid.length()),BLASPREAL(grid.data()),BLASINT(mnlIntOne),BLASPREAL(gridTemp.data()+1),BLASINT(mnlIntOne));
      gridTemp[grid.length()+1] = 1;

      GaussLobattoLegendreWeights(temp,gridTemp);
      weight = temp[utilities::Range::colon(1,grid.length())];
    }

    basics::Matrix GLL::LagrangeDerivativeMatrix(int N)
    {
      basics::Matrix result("Lagrange Derivative Matrix",N,N);
      LagrangeDerivativeMatrix(result);

      return( result );
    }

    void GLL::LagrangeDerivativeMatrix(basics::Matrix& D)
    {
      basics::Vector grid("Gauss Lobatto Legendre grid",D.rows(),false);
      GaussLobattoLegendreGrid(grid);
      LagrangeDerivativeMatrix(D,grid);
    }

    void GLL::LagrangeDerivativeMatrix(basics::Matrix& A, const basics::Vector& grid)
    {
      assert( A.rows() == A.cols() && A.rows() == grid.length());

      int N=grid.length();
      for( int i=0;i<N;++i )
        for( int j=0;j<N;++j ) {
          if( i == j )
            A[j][i] = 0;
          else
            A[j][i] = Legendre(grid[i],N-1)/(Legendre(grid[j],N-1)*(grid[i]-grid[j]));
        }
      A[0][0] = -N*(N-1)/4.f;
      A[N-1][N-1] = N*(N-1)/4.f;
    }

    void GLL::LagrangeDerivativeMatrixPeriodic(basics::Matrix& D)
    {
      basics::Vector grid("Gauss Lobatto Legendre grid",D.rows(),false);
      GaussLobattoLegendreGridPeriodic(grid);
      LagrangeDerivativeMatrixPeriodic(D,grid);
    }

    void GLL::LagrangeDerivativeMatrixPeriodic(basics::Matrix& A, const basics::Vector& grid)
    {
      basics::Matrix M2 = LagrangeDerivativeMatrix(A.cols()+1);
      A = M2.submatrix(utilities::Range::colon(0,M2.rows()-2),utilities::Range::colon(0,M2.cols()-2));
      A[0] += M2[M2.cols()-1][utilities::Range::colon(0,M2.rows()-2)];
    }

    void GLL::LagrangeDerivativeMatrixHomogenous(basics::Matrix& A, const basics::Vector& grid)
    {
      basics::Matrix M2 = LagrangeDerivativeMatrix(A.rows()+2);
      A = M2.submatrix(utilities::Range::colon(1,M2.rows()-2),utilities::Range::colon(1,M2.cols()-2));
    }

    /* returns a interpolation matrix from gridFrom to M GLL points */
    basics::Matrix GLL::interpolationMatrix(const basics::Vector& gridFrom, int M)
    {
      basics::Vector gridTo("GLL grid M",M);
      GaussLobattoLegendreGrid(gridTo);

      return( GLL::interpolationMatrix(gridFrom,gridTo) );
    }

    /* returns a Lagrange interpolant matrix from grid gridFrom to grid gridTo */
    basics::Matrix GLL::interpolationMatrix(const basics::Vector& gridFrom, const basics::Vector& gridTo)
    {
      basics::Matrix result("Lagrange interpolant",gridTo.length(),gridFrom.length());
      for( int j=0;j<gridTo.length();++j )
        for( int i=0;i<gridFrom.length();++i )
          result.data()[i][j] = Lagrange(gridTo[j],i,gridFrom);

      return( result );
    }

    /* returns a interpolant matrix from periodic grid gridFrom to grid gridTo */
    basics::Matrix GLL::periodicInterpolationMatrix(const basics::Vector& gridFrom, const basics::Vector& gridTo)
    {
      basics::Matrix result("GLL interpolant",gridTo.length(),gridFrom.length());
      basics::Vector temp("tempgrid",gridFrom.length()+1);
      BLASRPFX(copy,BLASINT(gridFrom.length()),BLASPREAL(gridFrom.data()),BLASINT(mnlIntOne),BLASPREAL(temp.data()),BLASINT(mnlIntOne));
      temp[temp.length()-1] = gridFrom[gridFrom.length()-1]+temp[1]-temp[0];
      for( int j=0;j<gridTo.length();++j )
        for( int i=0;i<gridFrom.length();++i )
          if (i == 0)
            result.data()[i][j] = Lagrange(gridTo[j],0,temp)+Lagrange(gridTo[j],gridFrom.length(),temp);
          else
            result.data()[i][j] = Lagrange(gridTo[j],i,temp);

      return( result );
    }
  } // namespace utilities
} // namespace mnl

