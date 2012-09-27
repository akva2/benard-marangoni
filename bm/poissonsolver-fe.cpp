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

#include "poissonsolver-fe.h"
#include "mnl/range.h"
#include "mnl/gll.h"
#include "mnl/matrices.h"
#include "mnl/memtracker.h"

#include <numeric>

using namespace mnl;
using namespace std;

namespace linearElements {
  poissonSolver::poissonSolver(const BC& h, const basics::Vector& grid,
      const basics::Vector& scaleX, const basics::Vector& scaleY,
      bool L2, bool doZ, const BC& bz) :
    poissonSolverT<basics::Matrix,basics::Matrices>(h,1,1,0),
    m_buffer(NULL), m_L2(L2), m_bz(bz)
  {
    m_eigsX = m_eigsY = m_eigsZ = NULL;
    m_gridZ = m_weightZ = m_weightX = m_weightY = NULL;
    m_gridX = m_gridY = m_gridZ = NULL;
    m_Qx = m_Qy = m_Qz = m_Az = NULL;
    m_Dx = m_Dy = m_Dz = NULL;
    m_Bx = m_By = m_Bz = NULL;

    if( L2 )
      setupL2(h,grid,scaleX,scaleY,doZ);
    else
      setupH1(h,grid,scaleX,scaleY,doZ,bz);
  }

  poissonSolver::~poissonSolver()
  {
    delete m_Bx;
    delete m_By;
    delete m_Bz;

    utilities::g_manager.unlock(m_buffer);
  }

  basics::Vector* poissonSolver::setupL2Grid(const BC& h,
      const basics::Vector& grid,
      const basics::Vector& scale)
  {
    int len = scale.length()*grid.length();

    if( h == Homogenous )
      len += 2;
    if( h == HomogenousLeft )
      len += 1;

    Real L = scale.sum();

    basics::Vector* result = new basics::Vector("FE grid",len);
    Real start = 0;
    int index = (h==HomogenousLeft?0:1);
    for( int i=0;i<scale.length();++i ) {
      result->assign(utilities::Range::colon(index,index+grid.length()-1),
          (grid+1)*scale[i]/(2*L)+start);
      index += grid.length();

      start += scale[i]/L;
    }
    if( h == Homogenous )
      (*result)[0] = -(*result)[2];
    if( h == Homogenous || h == HomogenousLeft )
      (*result)[len-1] = 1+(1-(*result)[len-3]);

    return result;
  }

  void poissonSolver::setupL2(const BC& h, const basics::Vector& grid,
      const basics::Vector& scaleX,
      const basics::Vector& scaleY, bool doZ)
  {
    m_gridX = setupL2Grid(Homogenous,grid,scaleX);
    m_gridY = setupL2Grid(h==Homogenous?h:HomogenousLeft,grid,scaleY);

    assembleAndMask(m_Ax,m_Bx,*m_gridX,Homogenous);
    assembleAndMask(m_Ay,m_By,*m_gridY,h==Homogenous?h:HomogenousRight);
    diagonalizeOperator(m_Ax,m_Qx,m_eigsX,*m_Bx);
    diagonalizeOperator(m_Ay,m_Qy,m_eigsY,*m_By);
    int lenZ = 1;
    if( doZ ) {
      m_gridZ = new basics::Vector(grid);
      *m_gridZ += 1;
      *m_gridZ /= 2;
      assembleAndMask(m_Az,m_Bz,*m_gridZ,Nonhomogenous);
      lenZ = m_Az->rows();
      diagonalizeOperator(m_Az,m_Qz,m_eigsZ,*m_Bz);
      m_bz = Nonhomogenous;
    }

    m_buffer = utilities::g_manager.aquireMatrices("temp",m_Ax->rows(),m_Ay->rows(),lenZ);
  }

  void poissonSolver::setupH1(const BC& h, const basics::Vector& grid,
      const basics::Vector& scaleX,
      const basics::Vector& scaleY, bool doZ, const BC& bz)
  {
    m_gridX = setupH1Grid(Homogenous,grid,scaleX);
    m_gridY = setupH1Grid(h==Homogenous?h:HomogenousLeft,grid,scaleY);

    assembleAndMask(m_Ax,m_Bx,*m_gridX,Homogenous);
    assembleAndMask(m_Ay,m_By,*m_gridY,h==Nonhomogenous?HomogenousRight:Homogenous);
    diagonalizeOperator(m_Ax,m_Qx,m_eigsX,*m_Bx);
    diagonalizeOperator(m_Ay,m_Qy,m_eigsY,*m_By);
    int lenZ = 1;
    if( doZ ) {
      m_gridZ = new basics::Vector(grid);
      *m_gridZ += 1;
      *m_gridZ /= 2;
      assembleAndMask(m_Az,m_Bz,*m_gridZ,bz);
      lenZ = m_Az->rows();
      diagonalizeOperator(m_Az,m_Qz,m_eigsZ,*m_Bz);
    }

    m_buffer = utilities::g_manager.aquireMatrices("temp",m_Ax->rows(),m_Ay->rows(),lenZ);
  }


  basics::Vector* poissonSolver::setupH1Grid(const BC& h, const basics::Vector& grid,
      const basics::Vector& scale)
  {
    int len = scale.length()*grid.length();
    len -= scale.length()-1; // remove shared points between element boundaries
    if( h == Homogenous )
      len += 2;
    if( h == HomogenousLeft )
      len += 1;

    Real L = scale.sum();

    basics::Vector* result = new basics::Vector("FE grid",len);
    Real start = 0;
    int index = (h==HomogenousLeft?0:1);
    for( int i=0;i<scale.length();++i ) {
      result->assign(utilities::Range::colon(index,index+grid.length()-1),
          (grid+1)*scale[i]/(2*L)+start);
      index += grid.length()-1;

      start += scale[i]/L;
    }
    if( h == Homogenous )
      (*result)[0] = -(*result)[2];
    if( h == Homogenous || h == HomogenousLeft )
      (*result)[len-1] = 1+(1-(*result)[len-3]);

    return result;
  }

#define ADD_ELEMENT_MATRIX(A,Ak,i,h) \
  A[  i][  i] += h*Ak[0][0]; \
  A[  i][i+1] += h*Ak[0][1]; \
  A[i+1][  i] += h*Ak[1][0]; \
  A[i+1][i+1] += h*Ak[1][1];

  void poissonSolver::assemble(basics::Matrix& A, basics::Matrix& B,
      const basics::Vector& grid)
  {
    const Real Ak[2][2] = {{1.,-1.f},		   {-1.f,1.f}};
    const Real Bk[2][2] = {{1.f/3.f,1.f/6.f},{1.f/6.f,1.f/3.f}};

    int len = grid.length()-1;
    for( int i=0;i<len;++i ) {
      Real hk = grid[i+1]-grid[i];
      ADD_ELEMENT_MATRIX(A,Ak,i,Real(1)/hk)
        ADD_ELEMENT_MATRIX(B,Bk,i,hk)
    }
  }

  void poissonSolver::assembleAndMask(basics::Matrix*& A,
      basics::Matrix*& B,
      const basics::Vector& grid, BC h)
  {
    basics::Matrix At("temp",grid.length(),grid.length());
    basics::Matrix Bt("temp",grid.length(),grid.length());
    assemble(At,Bt,grid);
    if( h == Homogenous ) {
      A = new basics::Matrix(At.submatrix(utilities::Range::colon(1,At.rows()-2),
            utilities::Range::colon(1,At.cols()-2)));
      B =	new basics::Matrix(Bt.submatrix(utilities::Range::colon(1,Bt.rows()-2),
            utilities::Range::colon(1,Bt.cols()-2)));
    }
    if( h == HomogenousLeft ) {
      A = new basics::Matrix(At.submatrix(utilities::Range::colon(1,At.rows()-1),
            utilities::Range::colon(1,At.cols()-1)));
      B =	new basics::Matrix(Bt.submatrix(utilities::Range::colon(1,Bt.rows()-1),
            utilities::Range::colon(1,Bt.cols()-1)));
    }
    if( h == HomogenousRight ) {
      A = new basics::Matrix(At.submatrix(utilities::Range::colon(0,At.rows()-2),
            utilities::Range::colon(0,At.cols()-2)));
      B =	new basics::Matrix(Bt.submatrix(utilities::Range::colon(1,Bt.rows()-1),
            utilities::Range::colon(1,Bt.cols()-1)));
    }
    if( h == Nonhomogenous ) {
      A = new basics::Matrix(At.submatrix(utilities::Range::colon(0,At.rows()-1),
            utilities::Range::colon(0,At.cols()-1)));
      B =	new basics::Matrix(Bt.submatrix(utilities::Range::colon(0,Bt.rows()-1),
            utilities::Range::colon(0,Bt.cols()-1)));
    }
  }

  void poissonSolver::solve(basics::Matrix& F, const Real mu, 
      const Real nu, const Real scale) const
  {
    solve(F,(*m_buffer)[0],mu,nu,scale);
  }

  void poissonSolver::solve(basics::Matrix& F, basics::Matrix& buffer, const Real mu,
      const Real Lx, const Real Ly) const
  {
    LOG("running poissonSolver2D-FE");

    basics::multTranspose(buffer,F,*m_Qy,'N','N');
    basics::multTranspose(F,*m_Qx,buffer,'T','N');

    Real scaleX = Ly/Lx;
    Real scaleY = Lx/Ly;
    Real smu = Lx*Ly*mu;

    /* do the algebra */
    for( int j=0;j<F.cols();++j )
      for( int k=0;k<F.rows();++k )
        F[j][k] /= (scaleY*(*m_eigsY)[j]+scaleX*(*m_eigsX)[k]+smu);

    basics::multTranspose(buffer,F,*m_Qy,'N','T');
    basics::multTranspose(F,*m_Qx,buffer,'N','N');
  }

  void poissonSolver::solve(basics::Matrices& F, const Real mu, const Real nu) const
  {
    solve(F,*m_buffer,mu,nu,1);
  }

  void poissonSolver::solve(basics::Matrices& F, basics::Matrices& buffer, const Real Lx,
      const Real Ly, const Real Lz, const Real mu) const
  {
    LOG("running poissonSolver3D-FE");

    int index = (m_L2?0:1);
    Real scaleX = Ly*Lz/(Lx);
    Real scaleY = Lx*Lz/(Ly);
    Real scaleZ = Lx*Ly/(Lz);
    Real scaleM = Lx*Ly*Lz*mu;

    int skip=m_L2?0:2;
    int l=index;
    int max=F.matrices()-index;
    if( m_bz == HomogenousLeft ) {
      skip = l = index = 1;
      max = F.matrices();
    }
    if( m_bz == HomogenousRight ) {
      skip = 3;
      l = 0;
      index = 0;
      max = F.matrices()-1;
    }
    if( m_bz == Nonhomogenous ) {
      skip = 0;
      max = F.matrices();
      index = 0;
    }
    basics::applyLocalGlobal(buffer,F,*m_Qz,'N','N',skip);
    for( l;l<max;++l ) {
      basics::multTranspose(F[l],buffer[l],*m_Qy,'N','N');
      basics::multTranspose(buffer[l],*m_Qx,F[l],'T','N');

      /* do the algebra */
      for( int j=0;j<F.cols();++j )
        for( int k=0;k<F.rows();++k ) {
          buffer[l][j][k] /= (scaleX*(*m_eigsX)[k]+
              scaleY*(*m_eigsY)[j]+
              scaleZ*(*m_eigsZ)[l-index]+scaleM);
        }

      basics::multTranspose(F[l],buffer[l],*m_Qy,'N','T');
      basics::multTranspose(buffer[l],*m_Qx,F[l],'N','N');
    }
    basics::applyLocalGlobal(F,buffer,*m_Qz,'N','T',skip);
  }
}

