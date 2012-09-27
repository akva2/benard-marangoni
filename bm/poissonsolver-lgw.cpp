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

#include "poissonsolver-lgw.h"
#include "mnl/range.h"
#include "mnl/gll.h"
#include "mnl/matrices.h"
#include "mnl/memtracker.h"

#include "legendrelegendrew.h"

#undef M_PI
#define M_PI M_PIl

using namespace mnl;
using namespace std;

/**
 * poissonsolver.cpp - solves the three dimensional helmholtz/poisson problem 
 * with periodic, homogenous and periodic b.c's
 *
 *		- this class solves the problem
 *
 *		-nu^2(nabla)^2 u + mu*u= f in (0,2*pi)x(-1,1)x(0,2*pi)
 *		with 	u(x,-1,z) 	= u(x,1,z) 		= 0 
 *				u(2*pi,y,z) = u(2*pi,y,z) 	= 0
 *				u(x,y,2*pi) = u(x,y,2*pi)		= 0
 *		or without boundary conditions
 *	the algorithm used is Legendre Galerkin 
 **/

namespace legendreLegendreW {

  static void generateA(basics::Matrix& A, basics::Matrix& D, basics::Vector& weight,
      basics::Vector& grid, const basics::Vector& grid1, Real scale, int elems)
  {
    A.clear();
    int startcol=0;
    Real L;
    if( scale == -1.f && elems==1 )
      L = 2;
    else
      L = scale/elems;
    for (int e=0;e<elems;++e ) {
      for( int j=0;j<D.cols();++j ) {
        for( int k=0;k<D.rows();++k )
          for( int l=0;l<D.cols();++l )
            A[k+startcol][j+startcol] += weight[l]*D[j][l]*D[k][l]*2/L;
        if( scale == -1.f && elems==1 )
          grid[j+startcol] = grid1[j];
        else
          grid[j+startcol] = (grid1[j]+1)/2*L+L*e;
      }
      startcol += weight.length()-1;
    }
  }

  static void generateB(basics::Matrix& B, basics::Vector& weight, Real scale, int elems)
  {
    B.clear();
    int startcol=0;
    Real L;
    if( scale == -1.f && elems==1 )
      L = 2.f;
    else
      L = scale/elems;
    for (int e=0;e<elems;++e ) {
      for( int j=0;j<weight.length();++j )
        B[j+startcol][j+startcol] += weight[j]*L/2;
      startcol += weight.length()-1;
    }
  }

  poissonSolver::poissonSolver(const int Nx, const int Ny, const int Nz, const BC h,
      bool generateOperator, Real scaleX, Real scaleY, Real scaleZ,
      int elemX, int elemY, int elemZ) :
    poissonSolverT<basics::Matrix,basics::Matrices>(h,scaleX,scaleY,scaleZ)
  {
    if( generateOperator ) {
      LOG("Setting up grids");

      m_eigsX = m_eigsY = m_eigsZ = NULL;
      m_Qx = m_Qy = m_Qz = NULL;
      m_gridX = new basics::Vector("Grid in X",(Nx+1)*elemX-(elemX-1));
      m_gridY = new basics::Vector("Gauss Lobatto Legendre grid",(Ny+1)*elemY-(elemY-1));
      if( Nz )
        m_gridZ = new basics::Vector("Grid in Z",Nz+1);
      else
        m_gridZ = NULL;
      m_weightX = new basics::Vector("Gauss Lobatto Legendre weights",Nx+1);
      m_weightY = new basics::Vector("Gauss Lobatto Legendre weights",Ny+1);
      if( Nz )
        m_weightZ = new basics::Vector("Gauss Lobatto Legendre weights",Nz+1);
      else
        m_weightZ = NULL;

      m_Dx = new basics::Matrix("Lagrange Derivative Matrix",Nx+1,Nx+1);
      m_Dy = new basics::Matrix("Lagrange Derivative Matrix",Ny+1,Ny+1);
      if( Nz )
        m_Dz = new basics::Matrix("Lagrange Derivative Matrix",Nz+1,Nz+1);
      else {
        m_Az = m_Dz = NULL;
      }

      LOG("Generating GLL quadrature");
      basics::Vector gridX("fool",Nx+1);
      utilities::GLL::GaussLobattoLegendreGrid(gridX);
      utilities::GLL::GaussLobattoLegendreWeights(*m_weightX,gridX);
      utilities::GLL::LagrangeDerivativeMatrix(*m_Dx);
      utilities::GLL::GaussLobattoLegendreGrid(*m_gridY);
      utilities::GLL::GaussLobattoLegendreWeights(*m_weightY,*m_gridY);
      utilities::GLL::LagrangeDerivativeMatrix(*m_Dy,*m_gridY);

      if( Nz ) {
        utilities::GLL::GaussLobattoLegendreGrid(*m_gridZ);
        utilities::GLL::GaussLobattoLegendreWeights(*m_weightZ,*m_gridZ);
        utilities::GLL::LagrangeDerivativeMatrix(*m_Dz,*m_gridZ);
      }

      LOG("Generating Laplacian operators");
      basics::Matrix Ax("protoA X",(Nx+1)*elemX-(elemX-1),(Nx+1)*elemX-(elemX-1));
      basics::Matrix Bx("protoB X",(Nx+1)*elemX-(elemX-1),(Nx+1)*elemX-(elemX-1));
      generateA(Ax,*m_Dx,*m_weightX,*m_gridX,gridX,scaleX,elemX);
      generateB(Bx,*m_weightX,scaleX,elemX);

      basics::Matrix Ay("protoA Y",(Ny+1)*elemY-(elemY-1),(Ny+1)*elemY-(elemY-1));
      basics::Matrix By("protoB Y",(Ny+1)*elemY-(elemY-1),(Ny+1)*elemY-(elemY-1));
      basics::Vector gridY("fool",Ny+1);
      utilities::GLL::GaussLobattoLegendreGrid(gridY);
      generateA(Ay,*m_Dy,*m_weightY,*m_gridY,gridY,scaleY,elemY);
      generateB(By,*m_weightY,scaleY,elemY);
      basics::Matrix Az("protoA Z",(Nz+1)*elemZ-(elemZ-1),(Nz+1)*elemZ-(elemZ-1));
      basics::Matrix Bz("protoB Z",(Nz+1)*elemZ-(elemZ-1),(Nz+1)*elemZ-(elemZ-1));
      if(	Nz ) {
        //            generateA(Az,*m_Dz,*m_weightZ,elemZ);
        //            generateB(Bz,*m_weightZ,scaleZ,elemZ);
      }
      {
        basics::Matrix* B=NULL;
        if( h == Nonhomogenous || h == DoNothing ) {
          m_Ax = new basics::Matrix(Ax);
          B = new basics::Matrix(Bx);
        }
        if( h == Homogenous  ) {
          utilities::Range r = utilities::Range::colon(1,Ax.rows()-2);
          m_Ax = new basics::Matrix(Ax.submatrix(r,r));
          B = new basics::Matrix(Bx.submatrix(r,r));

        }
        if( h == HomogenousLeft ) {
          utilities::Range r = utilities::Range::colon(1,Ax.rows()-1);
          m_Ax = new basics::Matrix(Ax.submatrix(r,r));
          B = new basics::Matrix(Bx.submatrix(r,r));
        }
        diagonalizeOperator(m_Ax,m_Qx,m_eigsX,*B);
        delete B;
      }
      {
        basics::Matrix* B=NULL;
        if( h == Nonhomogenous || h == DoNothing ) {
          m_Ay = new basics::Matrix(Ay);
          B = new basics::Matrix(By);
        }
        if( h == Homogenous  ) {
          utilities::Range r = utilities::Range::colon(1,Ay.rows()-2);
          m_Ay = new basics::Matrix(Ay.submatrix(r,r));
          B = new basics::Matrix(By.submatrix(r,r));
        }
        if( h == HomogenousLeft ) {
          utilities::Range r = utilities::Range::colon(1,Ay.rows()-1);
          m_Ay = new basics::Matrix(Ay.submatrix(r,r));
          B = new basics::Matrix(By.submatrix(r,r));
        }
        diagonalizeOperator(m_Ay,m_Qy,m_eigsY,*B);
        delete B;
      }
      if( Nz ) {
        basics::Matrix* B=NULL;
        if( h == Nonhomogenous || h == DoNothing ) {
          m_Az = new basics::Matrix(Az);
          B = new basics::Matrix(Bz);
        }
        if( h == Homogenous  ) {
          utilities::Range r = utilities::Range::colon(1,Az.rows()-2);
          m_Az = new basics::Matrix(Az.submatrix(r,r));
          B = new basics::Matrix(Bz.submatrix(r,r));
        }
        if( h == HomogenousLeft ) {
          utilities::Range r = utilities::Range::colon(1,Az.rows()-1);
          m_Az = new basics::Matrix(Az.submatrix(r,r));
          B = new basics::Matrix(Bz.submatrix(r,r));
        }
        diagonalizeOperator(m_Az,m_Qz,m_eigsZ,*B);
        delete B;
      }
      m_buffer = utilities::g_manager.aquireMatrices("temp",Ax.cols(),Ay.cols(),Nz?Nz:1);
    } else {
      m_buffer = utilities::g_manager.aquireMatrices("temp",Nx+1,Ny+1,Nz?Nz:1);
      m_gridX = m_gridY = m_gridZ = m_weightX = m_weightY 
        = m_weightZ = m_eigsX = m_eigsY = m_eigsZ = NULL;
      m_Ax    = m_Ay = m_Az = m_Qx = m_Qy = m_Qz = m_Dx = m_Dy = m_Dz = NULL;
    }
  }

  poissonSolver::~poissonSolver()
  {
    utilities::g_manager.unlock(m_buffer);
  }

  void poissonSolver::solve(basics::Matrix& F, const Real mu, const Real nu, const Real scale) const
  {
    LOG("running poissonSolver2D-LGW");

    basics::Matrix* temp = &((*m_buffer)[0]);

    int startindex = (m_b==Homogenous?1:0); // we still have the zeros at the end points to ease optimizations

    BLASRPFX(gemm,BLASCHAR(mnlNoTrans),            BLASCHAR(mnlNoTrans),
        BLASINT(F.rows()),               BLASINT(m_Qy->rows()),
        BLASINT(m_Qy->cols()),           BLASREAL(mnlRealOne),
        BLASPREAL(F.data()[startindex]), BLASINT(F.rows()),
        BLASPREAL(m_Qy->data()[0]),      BLASINT(m_Qy->cols()),
        BLASREAL(mnlRealZero),           BLASPREAL(temp->data()[0]),
        BLASINT(F.rows()));
    if( m_b != Homogenous )
      BLASRPFX(gemm, BLASCHAR(mnlTrans),         BLASCHAR(mnlNoTrans),
          BLASINT(m_Qx->cols()),      BLASINT(temp->cols()),
          BLASINT(m_Qx->rows()),      BLASREAL(mnlRealOne),
          BLASPREAL(m_Qx->data()[0]), BLASINT(m_Qx->rows()),
          BLASPREAL(temp->data()[0]), BLASINT(m_Qx->rows()),
          BLASREAL(mnlRealZero),      BLASPREAL(F.data()[0]),
          BLASINT(m_Qx->cols()));

    /* do the algebra */
    for( int j=startindex;j<F.cols()-startindex;++j ) {
      if( m_b == Homogenous )
        BLASRPFX(gemv, BLASCHAR(mnlTrans),                      BLASINT(m_Qx->rows()),
            BLASINT(m_Qx->cols()),                   BLASREAL(mnlRealOne),
            BLASPREAL(m_Qx->data()[0]),              BLASINT(m_Qx->rows()),
            BLASPREAL(temp->data()[j-startindex]+1), BLASINT(mnlIntOne),
            BLASREAL(mnlRealZero),                   BLASPREAL(F.data()[j]+1),
            BLASINT(mnlIntOne));
      for( int k=startindex;k<F.rows()-startindex;++k ) {
        if (m_b == Nonhomogenous && k == 0 && j == 0 )
          F[j][k] = 0;
        else
          F[j][k] /= (nu*(*m_eigsY)[j-startindex]+nu*(*m_eigsX)[k-startindex]+mu);
      }
    }
    BLASRPFX(gemm, BLASCHAR(mnlNoTrans),            BLASCHAR(mnlTrans),
        BLASINT(F.rows()),               BLASINT(m_Qy->rows()),
        BLASINT(m_Qy->cols()),           BLASREAL(mnlRealOne),
        BLASPREAL(F.data()[startindex]), BLASINT(F.rows()),
        BLASPREAL(m_Qy->data()[0]),      BLASINT(m_Qy->cols()),
        BLASREAL(mnlRealZero),           BLASPREAL(temp->data()[0]),
        BLASINT(F.rows()));

    if( m_b == Homogenous ) {
      for( int j=startindex;j<F.cols()-startindex;++j ) {
        BLASRPFX(gemv, BLASCHAR(mnlNoTrans),                            BLASINT(m_Qx->rows()),
            BLASINT(m_Qx->cols()),                           BLASREAL(mnlRealOne),
            BLASPREAL(m_Qx->data()[0]),                      BLASINT(m_Qx->rows()),
            BLASPREAL(temp->data()[j-startindex]+1),         BLASINT(mnlIntOne),
            BLASREAL(mnlRealZero), BLASPREAL(F.data()[j]+1), BLASINT(mnlIntOne));
      }
    } else
      BLASRPFX(gemm, BLASCHAR(mnlNoTrans),       BLASCHAR(mnlNoTrans),
          BLASINT(m_Qx->rows()),      BLASINT(temp->cols()),
          BLASINT(temp->rows()),      BLASREAL(mnlRealOne),
          BLASPREAL(m_Qx->data()[0]), BLASINT(m_Qx->rows()), 
          BLASPREAL(temp->data()[0]), BLASINT(temp->rows()),
          BLASREAL(mnlRealZero),      BLASPREAL(F.data()[0]),
          BLASINT(m_Qx->rows()));

    if( m_b == Homogenous )
      applyHomogenousDirichletBC(F);
  }

  void poissonSolver::solveMult(basics::Matrix& F, const Real mu, const Real nu, const Real scale) const
  {
    LOG("running poissonSolver2D-LGW");

    basics::Matrix* temp = &((*m_buffer)[0]);

    basics::multTranspose(*temp,F,*m_Qy,'N','N');
    basics::multTranspose(F,*m_Qx,*temp,'T','N');

    /* do the algebra */
    for( int j=0;j<F.cols();++j ) {
      for( int k=0;k<F.rows();++k ) {
        if(k && j )
          //            if( (*m_eigsY)[j] > 1.e-14 && (*m_eigsX)[k] > 1.e-14 )
          F[j][k] /= ((*m_eigsY)[j]*(*m_eigsX)[k]);
        else
          F[j][k] = 0;
      }
    }

    basics::multTranspose(*temp,*m_Qx,F,'N','N');
    basics::multTranspose(F,*temp,*m_Qy,'N','T');
  }

  void poissonSolver::solve(basics::Matrices& F, const Real mu, const Real nu) const
  {
    LOG("running poissonSolver3D-LGW");

    int startindex = (m_b==Homogenous?1:0); // we still have the zeros at the end points to ease optimizations

    basics::applyLocalGlobal(*m_buffer,F,*m_Qz,'N','N',0);
    for( int l=0;l<F.matrices();++l ) {
      BLASRPFX(gemm,BLASCHAR(mnlNoTrans),            BLASCHAR(mnlNoTrans),
          BLASINT(F.rows()),               BLASINT(m_Qy->rows()),
          BLASINT(m_Qy->cols()),           BLASREAL(mnlRealOne),
          BLASPREAL((*m_buffer)[l].data()[startindex]), BLASINT(F.rows()),
          BLASPREAL(m_Qy->data()[0]),      BLASINT(m_Qy->cols()),
          BLASREAL(mnlRealZero),           BLASPREAL(F[l].data()[0]),
          BLASINT(F.rows()));
      BLASRPFX(gemm,BLASCHAR(mnlTrans),         BLASCHAR(mnlNoTrans),
          BLASINT(m_Qx->cols()),      BLASINT(F.cols()),
          BLASINT(m_Qx->rows()),      BLASREAL(mnlRealOne),
          BLASPREAL(m_Qx->data()[0]), BLASINT(m_Qx->rows()),
          BLASPREAL(F[l].data()[0]),  BLASINT(m_Qx->rows()),
          BLASREAL(mnlRealZero),      BLASPREAL((*m_buffer)[l].data()[0]),
          BLASINT(m_Qx->cols()));

      /* do the algebra */
      for( int j=startindex;j<F.cols()-startindex;++j ) {
        for( int k=startindex;k<F.rows()-startindex;++k ) {
          if (m_b == Nonhomogenous && ((*m_eigsX)[k] < 1.e-14 && 
                (*m_eigsY)[j] < 1.e-14)&&
              (*m_eigsZ)[l] < 1.e-14   ) 
            (*m_buffer)[l][j][k] = 0;
          else
            (*m_buffer)[l][j][k] /= (nu*(*m_eigsY)[j-startindex]+nu*(*m_eigsX)[k-startindex]+nu*(*m_eigsZ)[l]+mu);
        }
      }
      BLASRPFX(gemm, BLASCHAR(mnlNoTrans),            BLASCHAR(mnlTrans),
          BLASINT(F.rows()),               BLASINT(m_Qy->rows()),
          BLASINT(m_Qy->cols()),           BLASREAL(mnlRealOne),
          BLASPREAL((*m_buffer)[l].data()[startindex]), BLASINT(F.rows()),
          BLASPREAL(m_Qy->data()[0]),      BLASINT(m_Qy->cols()),
          BLASREAL(mnlRealZero),           BLASPREAL(F[l].data()[0]),
          BLASINT(F.rows()));

      BLASRPFX(gemm, BLASCHAR(mnlNoTrans),       BLASCHAR(mnlNoTrans),
          BLASINT(m_Qx->rows()),      BLASINT(F.cols()),
          BLASINT(F.rows()),      	   BLASREAL(mnlRealOne),
          BLASPREAL(m_Qx->data()[0]), BLASINT(m_Qx->rows()), 
          BLASPREAL(F[l].data()[0]),  BLASINT(F.rows()),
          BLASREAL(mnlRealZero),      BLASPREAL((*m_buffer)[l].data()[0]),
          BLASINT(m_Qx->rows()));
    }	
    basics::applyLocalGlobal(F,*m_buffer,*m_Qz,'N','T',0);
  }

}

