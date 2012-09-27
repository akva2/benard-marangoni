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

#include "poissonsolver-lg.h"
#include "mnl/range.h"
#include "mnl/gll.h"
#include "mnl/matrices.h"
#include "mnl/memtracker.h"
#include "legendrelegendre.h"

using namespace mnl;

/**
 * poissonsolver.cpp - solves the three dimensional helmholtz/poisson problem with periodic, homogenous and periodic b.c's
 *
 *		- this class solves the problem
 *
 *		-nu^2(nabla)^2 u + mu*u= f in (0,2*pi)x(-1,1)x(0,2*pi)
 *		with 	u(x,-1,z) = u(x,1,z) = 0 or without boundary conditions
 *			u(0,y,z) = u(2*pi,y,z)
 *			u(x,y,0) = u(x,y,2*pi)
 *
 *	the algorithm used is Legendre Galerkin 
 **/

namespace legendreLegendre {
  poissonSolver::poissonSolver(const int Nx, const int Ny, const int Nz, 
      const BC h, bool generateOperator,
      Real scaleX, Real scaleY) :
    poissonSolverT<basics::Matrix,basics::Matrices>(h)
  {
    if( generateOperator ) {
      LOG("Setting up grids");

      int n_Nz = Nz;
      if( n_Nz == 0 )
        n_Nz = 1;

      m_eigsX = m_eigsY = m_eigsZ = NULL;
      m_Qx = m_Qy = m_Qz = NULL;
      m_gridX = new basics::Vector("Grid in X",Nx,false);
      m_gridY = new basics::Vector("Gauss Lobatto Legendre grid",Ny+1);
      if( Nz )
        m_gridZ = new basics::Vector("Grid in Z",Nz);
      else
        m_gridZ = NULL;
      m_weightX = new basics::Vector("Gauss Lobatto Legendre weights",Nx+1);
      m_weightY = new basics::Vector("Gauss Lobatto Legendre weights",Ny+1);
      if( Nz )
        m_weightZ = new basics::Vector("Gauss Lobatto Legendre weights",Nz);
      else
        m_weightZ = NULL;

      m_Ax = new basics::Matrix("Legendre Galerkin operator",Nx,Nx);
      m_Dx = new basics::Matrix("Lagrange Derivative Matrix",Nx+1,Nx+1);
      m_Ay = new basics::Matrix("Legendre Galerkin operator",Ny-1+2*h,Ny-1+2*h);
      m_Dy = new basics::Matrix("Lagrange Derivative Matrix",Ny+1,Ny+1);
      if( Nz ) {
        m_Az = new basics::Matrix("Legendre Galerkin operator",Nz,Nz);
        m_Dz = new basics::Matrix("Lagrange Derivative Matrix",Nz,Nz);
      } else {
        m_Az = m_Dz = NULL;
      }

      LOG("Generating GLL quadrature");
      utilities::GLL::GaussLobattoLegendreGridPeriodic(*m_gridX);
      utilities::GLL::GaussLobattoLegendreWeightsPeriodic(*m_weightX,*m_gridX);
      utilities::GLL::LagrangeDerivativeMatrix(*m_Dx);
      *m_gridX += 1;
      *m_gridX *= M_PI;

      utilities::GLL::GaussLobattoLegendreGrid(*m_gridY);
      utilities::GLL::GaussLobattoLegendreWeights(*m_weightY,*m_gridY);
      utilities::GLL::LagrangeDerivativeMatrix(*m_Dy,*m_gridY);

      if( Nz ) {
        utilities::GLL::GaussLobattoLegendreGridPeriodic(*m_gridZ);
        utilities::GLL::GaussLobattoLegendreWeightsPeriodic(*m_weightZ,*m_gridZ);
        utilities::GLL::LagrangeDerivativeMatrixPeriodic(*m_Dz,*m_gridZ);
        *m_gridZ += 1;
        *m_gridZ *= M_PI;
      }

      LOG("Generating Laplacian operators");
      {
        basics::Matrix Ax("temp",Nx+1,Nx+1);
        for( int j=0;j<Nx+1;++j )
          for( int k=0;k<Nx+1;++k ) {
            for( int l=0;l<Nx+1;l++ )
              Ax[k][j] += (*m_weightX)[l]*(*m_Dx)[j][l]*(*m_Dx)[k][l];
          }
        *m_Ax = makePeriodic(Ax,true);
        basics::Matrix B = makePeriodic(basics::Matrix::diag((*m_weightX)),true);
        B *= M_PI;
        *m_Ax *= Real(1)/M_PI;
        diagonalizeOperator(m_Ax,m_Qx,m_eigsX,B);
      }
      {
        for( int j=0;j<Ny-1+2*h;++j )
          for( int k=0;k<Ny-1+2*h;++k ) {
            (*m_Ay)[k][j] = 0;
            for( int l=0;l<Ny+1;l++ )
              (*m_Ay)[k][j] += (*m_weightY)[l]*(*m_Dy)[j+1-h][l]*(*m_Dy)[k+1-h][l];
          }
        basics::Matrix B = basics::Matrix::diag((*m_weightY)[utilities::Range::colon(1-m_b,Ny-1+m_b)]);
        diagonalizeOperator(m_Ay,m_Qy,m_eigsY,B);
      }
      if( Nz ) {
        for( int j=0;j<Nz;++j )
          for( int k=0;k<Nz;++k ) {
            (*m_Az)[k][j] = 0;
            for( int l=0;l<Nz;l++ ) {
              Real jbid = (*m_Dz)[j][l];
              Real kbid = (*m_Dz)[k][l];

              (*m_Az)[k][j] += (*m_weightZ)[l]*jbid*kbid;
            }
          }
        basics::Matrix B = basics::Matrix::diag(*m_weightZ);
        B *= M_PI;
        *m_Az *= Real(1)/M_PI;
        diagonalizeOperator(m_Az,m_Qz,m_eigsZ,B);
      }
    } else {
      m_gridX = m_gridY = m_gridZ = m_weightX = m_weightY = m_weightZ = m_eigsX = m_eigsY = m_eigsZ = NULL;
      m_Ax = m_Ay = m_Az = m_Qx = m_Qy = m_Qz = m_Dx = m_Dy = m_Dz = NULL;
    }
  }

  poissonSolver::~poissonSolver()
  {
  }

  void poissonSolver::solve(basics::Matrix& F, const Real mu, const Real nu, const Real scale) const
  {
    LOG("running poissonSolver2D-LG");

    basics::Matrix* temp  = utilities::g_manager.aquireMatrix("temp",F.rows(),F.cols()-2+2*m_b);
    basics::Matrix* temp2 = utilities::g_manager.aquireMatrix("temp",F.rows(),F.cols()-2+2*m_b);

    int startindex = (m_b==Homogenous?1:0); // we still have the zeros at the end points to ease optimizations
    for( int i=0;i<m_Qy->cols();++i )
      (*temp2)[i] = F[startindex+i];

    basics::multTranspose(*temp,*temp2,*m_Qy,'N','N');
    basics::multTranspose(*temp2,*m_Qx,*temp,'T','N');

    /* do the algebra */
    for( int k=0;k<temp2->rows();++k ) {
      for( int j=0;j<temp2->cols();++j ) {
        if (m_b == Nonhomogenous && (k == 0 && j == 0))
          (*temp2)[j][k] = 0;
        else
          (*temp2)[j][k] /= (nu*(*m_eigsY)[j]+nu*(*m_eigsX)[k]+mu);
      }
    }

    basics::multTranspose(*temp,*temp2,*m_Qy,'N','T');
    basics::multTranspose(*temp2,*m_Qx,*temp,'N','N');

    for( int i=0;i<m_Qy->cols();++i )
      F[startindex+i] = (*temp2)[i];

    if( m_b == Homogenous )
      applyHomogenousDirichletBC(F);

    utilities::g_manager.unlock(temp);
    utilities::g_manager.unlock(temp2);
  }

  //    void poissonSolver::solve(basics::Matrix& F, const Real mu, const Real nu, const Real scale) const
  //    {
  //        LOG("running poissonSolver2D-LG");

  //        basics::Matrix* temp  = utilities::g_manager.aquireMatrix("temp",F.rows(),F.cols()-2+2*m_b);
  //        basics::Matrix* temp2 = utilities::g_manager.aquireMatrix("temp",F.rows(),F.cols()-2+2*m_b);
  //        
  //        int startindex = (m_b==Homogenous?1:0); // we still have the zeros at the end points to ease optimizations
  //        for( int i=0;i<m_Qy->cols();++i )
  //            (*temp2)[i] = F[startindex+i];
  //    
  //        BLASRPFX(gemm, BLASCHAR(mnlNoTrans),                  BLASCHAR(mnlNoTrans), 
  //                       BLASINT(F.rows()),                     BLASINT(m_Qy->rows()),
  //                       BLASINT(m_Qy->cols()),                 BLASREAL(mnlRealOne),
  //                       BLASPREAL(F.data()[startindex]),       BLASINT(F.rows()),
  //                       BLASPREAL(m_Qy->data()[0]),            BLASINT(m_Qy->cols()),
  //                       BLASREAL(mnlRealZero),temp->data()[0], BLASINT(F.rows()));

  //        BLASRPFX(gemm, BLASCHAR(mnlTrans),         BLASCHAR(mnlNoTrans),
  //                       BLASINT(m_Qx->cols()),      BLASINT(temp->cols()),
  //                       BLASINT(m_Qx->rows()),      BLASREAL(mnlRealOne),
  //                       BLASPREAL(m_Qx->data()[0]), BLASINT(m_Qx->rows()),
  //                       BLASPREAL(temp->data()[0]), BLASINT(m_Qx->rows()),
  //                       BLASREAL(mnlRealZero),      BLASPREAL(F.data()[startindex]),
  //                       BLASINT(m_Qx->cols()));

  //        /* do the algebra */
  //        for( int k=0;k<F.rows();++k ) {
  //            for( int j=startindex;j<F.cols()-startindex;++j ) {
  //                if (m_b == Nonhomogenous && (k == 0 && j == 0))
  //                    F[j][k] = 0;
  //                else
  //                    F[j][k] /= (nu*(*m_eigsY)[j-startindex]+nu*(*m_eigsX)[k]+mu);
  //            }
  //        }

  //        BLASRPFX(gemm, BLASCHAR(mnlNoTrans),            BLASCHAR(mnlTrans),
  //                       BLASINT(F.rows()),               BLASINT(m_Qy->rows()),
  //                       BLASINT(m_Qy->cols()),           BLASREAL(mnlRealOne),
  //                       BLASPREAL(F.data()[startindex]), BLASINT(F.rows()),
  //                       BLASPREAL(m_Qy->data()[0]),      BLASINT(m_Qy->cols()),
  //                       BLASREAL(mnlRealZero),           BLASPREAL(temp->data()[0]),
  //                       BLASINT(F.rows()));

  //        BLASRPFX(gemm, BLASCHAR(mnlNoTrans),       BLASCHAR(mnlNoTrans),
  //                       BLASINT(m_Qx->rows()),      BLASINT(temp->cols()),
  //                       BLASINT(temp->rows()),      BLASREAL(mnlRealOne),
  //                       BLASPREAL(m_Qx->data()[0]), BLASINT(m_Qx->rows()),
  //                       BLASPREAL(temp->data()[0]), BLASINT(temp->rows()),
  //                       BLASREAL(mnlRealZero),      BLASPREAL(F.data()[startindex]),
  //                       BLASINT(m_Qx->rows()));

  //        if( m_b == Homogenous )
  //            applyHomogenousDirichletBC(F);

  //        utilities::g_manager.unlock(temp);
  //    }

  void poissonSolver::solve(basics::Matrices& F, const Real mu, const Real nu) const
  {
    int startindex = (m_b+1)%2; // we still have the zeros at the end points to ease optimizations

#pragma omp parallel for schedule(static)
    for( int l=0;l<F.matrices();++l ) {
      basics::Matrix* temp = utilities::g_manager.aquireMatrix("temp",F.rows(),F.cols()-2+2*m_b);

      BLASRPFX(gemm, BLASCHAR(mnlNoTrans),               BLASCHAR(mnlNoTrans),
          BLASINT(F.rows()),                  BLASINT(m_Qy->rows()),
          BLASINT(m_Qy->cols()),              BLASREAL(mnlRealOne),
          BLASPREAL(F[l].data()[startindex]), BLASINT(F.rows()),
          BLASPREAL(m_Qy->data()[0]),         BLASINT(m_Qy->cols()),
          BLASREAL(mnlRealZero),              BLASPREAL(temp->data()[0]),
          BLASINT(F.rows()));
      BLASRPFX(gemm, BLASCHAR(mnlTrans),         BLASCHAR(mnlNoTrans),
          BLASINT(m_Qx->cols()),      BLASINT(temp->cols()),
          BLASINT(m_Qx->rows()),      BLASREAL(mnlRealOne),
          BLASPREAL(m_Qx->data()[0]), BLASINT(m_Qx->rows()),
          BLASPREAL(temp->data()[0]), BLASINT(m_Qx->rows()),
          BLASREAL(mnlRealZero),      BLASPREAL(F[l].data()[startindex]),
          BLASINT(m_Qx->cols()));

      /* do the algebra */
      for( int k=0;k<F.rows();++k ) {
        for( int p=startindex;p<F.cols()-startindex;++p )
          F[l][p][k] /= (nu*(*m_eigsY)[p-startindex]+nu*(*m_eigsX)[k]+nu*(*m_eigsZ)[l]+mu);
      }

      BLASRPFX(gemm, BLASCHAR(mnlNoTrans),               BLASCHAR(mnlTrans),
          BLASINT(F.rows()),                  BLASINT(m_Qy->rows()),
          BLASINT(m_Qy->cols()),              BLASREAL(mnlRealOne),
          BLASPREAL(F[l].data()[startindex]), BLASINT(F.rows()),
          BLASPREAL(m_Qy->data()[0]),         BLASINT(m_Qy->cols()),
          BLASREAL(mnlRealZero),              BLASPREAL(temp->data()[0]),
          BLASINT(F.rows()));

      BLASRPFX(gemm, BLASCHAR(mnlNoTrans),       BLASCHAR(mnlNoTrans),
          BLASINT(m_Qx->rows()),      BLASINT(temp->cols()),
          BLASINT(temp->rows()),      BLASREAL(mnlRealOne),
          BLASPREAL(m_Qx->data()[0]), BLASINT(m_Qx->rows()),
          BLASPREAL(temp->data()[0]), BLASINT(temp->rows()),
          BLASREAL(mnlRealZero),      BLASPREAL(F[l].data()[startindex]),
          BLASINT(m_Qx->rows()));

      utilities::g_manager.unlock(temp);
      //            if( m_b == Homogenous )
      //                applyHomogenousDirichletBC(F);
    }
  }
}

