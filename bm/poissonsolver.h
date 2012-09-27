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

#ifndef POISSON_SOLVER_H_
#define POISSON_SOLVER_H_

#include "../mnl/vector.h"
#include "../mnl/matrices.h"
#include "../mnl/matrices.h"
#include "../mnl/function.h"

template <class T, class T3>
class poissonSolverT {
  public:
    enum BC {
      Homogenous 		= 0,
      Nonhomogenous 	= 1,
      DoNothing 		= 2,
      HomogenousLeft  = 3,
      HomogenousRight = 4
    };

    poissonSolverT(BC h=Homogenous, Real scaleX=2*M_PI, Real scaleY=-1, Real scaleZ=2*M_PI) :
      m_b(h), m_scaleX(scaleX), m_scaleY(scaleY), m_scaleZ(scaleZ)
  {
  }

    virtual ~poissonSolverT()
    {
      delete m_gridX;
      delete m_gridY;
      delete m_gridZ;
      delete m_weightX;
      delete m_weightY;
      delete m_weightZ;
      delete m_eigsX;
      delete m_eigsY;
      delete m_eigsZ;
      delete m_Ax;
      delete m_Ay;
      delete m_Az;
      delete m_Qx;
      delete m_Qy;
      delete m_Qz;
      delete m_Dx;
      delete m_Dy;
      delete m_Dz;
    }

    void source(T& F, const mnl::basics::function2D& f, const Real t=0, const Real nu=1.f, const Real mu=0.f) const
    {
      assert( F.rows() == m_gridX->length() && F.cols() == m_gridY->length() );

      for( int j=0;j<m_gridY->length();++j ) {
        Real y = (*m_gridY)[j];
        for( int k=0;k<m_gridX->length();++k ) {
          Real x = (*m_gridX)[k];
          F[j][k] = mu*mu*f.val(x,y,t)-nu*(f.diff2x(x,y,t)+f.diff2y(x,y,t));
        }
      }
    }

    void source(T3& F, const mnl::basics::function3D& f, const Real t=0, const Real nu=1.f, const Real mu=0.f) const
    {
      assert( F.rows() == m_gridX->length() && 
          F.cols() == m_gridY->length() && 
          F.matrices() == m_gridZ->length());

      for( int l=0;l<m_gridZ->length();++l )
        for( int j=0;j<m_gridY->length();++j )
          for( int k=0;k<m_gridX->length();++k )
            F[l][j][k] =  mu*mu*f.val((*m_gridX)[k],(*m_gridY)[j],(*m_gridZ)[l],t)
              -nu*( f.diff2x((*m_gridX)[k],(*m_gridY)[j],(*m_gridZ)[l],t)
                  +f.diff2y((*m_gridX)[k],(*m_gridY)[j],(*m_gridZ)[l],t)
                  +f.diff2z((*m_gridX)[k],(*m_gridY)[j],(*m_gridZ)[l],t));
    }

    virtual void solve(T& F, const Real mu=Real(0), const Real nu=Real(1), const Real scale=Real(1)) const = 0;
    virtual void solve(T3& F, const Real mu=Real(0), const Real nu=Real(1)) const = 0;

    inline void solve(mnl::basics::Field2<T>& F, const Real mu=Real(0), const Real nu=Real(1)) const
    {
      solve(F.X(),mu,nu);
      solve(F.Y(),mu,nu);
    }

    inline mnl::basics::Vector& gridX() 
    { 
      assert( m_gridX );
      return( *m_gridX ); 
    }

    inline const mnl::basics::Vector& gridX() const
    { 
      assert( m_gridX );
      return( *m_gridX ); 
    }

    inline mnl::basics::Vector& gridY() 
    { 
      assert( m_gridY );
      return( *m_gridY ); 
    }

    inline const mnl::basics::Vector& gridY() const
    { 
      assert( m_gridY );
      return( *m_gridY ); 
    }

    inline mnl::basics::Vector& gridZ() 
    { 
      assert( m_gridZ );
      return( *m_gridZ ); 
    }

    inline const mnl::basics::Vector& gridZ() const
    { 
      assert( m_gridZ );
      return( *m_gridZ ); 
    }

    inline mnl::basics::Vector& weightX() 
    { 
      assert( m_weightX );
      return( *m_weightX ); 
    }

    inline const mnl::basics::Vector& weightX() const
    { 
      assert( m_weightX );
      return( *m_weightX ); 
    }

    inline const mnl::basics::Vector& eigX() const
    { 
      assert( m_eigsX );
      return( *m_eigsX ); 
    }

    inline mnl::basics::Vector& weightY() 
    { 
      assert( m_weightY );
      return( *m_weightY ); 
    }

    inline const mnl::basics::Vector& weightY() const
    { 
      assert( m_weightY );
      return( *m_weightY ); 
    }

    inline const mnl::basics::Vector& eigY() const
    { 
      assert( m_eigsY );
      return( *m_eigsY ); 
    }

    inline mnl::basics::Vector& weightZ() 
    { 
      assert( m_weightZ );
      return( *m_weightZ ); 
    }

    inline const mnl::basics::Vector& weightZ() const
    { 
      assert( m_weightZ );
      return( *m_weightZ ); 
    }

    inline const mnl::basics::Vector& eigZ() const
    { 
      assert( m_eigsZ );
      return( *m_eigsZ ); 
    }

    inline mnl::basics::Matrix& Dx() 
    { 
      assert( m_Dx );
      return( *m_Dx ); 
    }

    inline const mnl::basics::Matrix& Dx() const
    { 
      assert( m_Dx );
      return( *m_Dx ); 
    }

    inline mnl::basics::Matrix& Dy() 
    { 
      assert( m_Dy );
      return( *m_Dy ); 
    }

    inline const mnl::basics::Matrix& Dy() const
    { 
      assert( m_Dy );
      return( *m_Dy ); 
    }

    inline mnl::basics::Matrix& Dz() 
    { 
      assert( m_Dz );
      return( *m_Dz ); 
    }

    inline const mnl::basics::Matrix* Dz() const
    { 
      assert( m_Dz );
      return( *m_Dz ); 
    }

    inline const mnl::basics::Matrix* Qx() const
    { 
      assert( m_Qx );
      return( m_Qx ); 
    }

    inline const mnl::basics::Matrix* Qy() const
    { 
      assert( m_Qy );
      return( m_Qy ); 
    }

    inline const mnl::basics::Matrix* Qz() const
    { 
      assert( m_Qz );
      return( m_Qz ); 
    }

    inline mnl::basics::Matrix& Ax() 
    { 
      assert( m_Ax );
      return( *m_Ax ); 
    }

    inline const mnl::basics::Matrix& Ax() const
    { 
      assert( m_Ax );
      return( *m_Ax ); 
    }

    inline mnl::basics::Matrix& Ay() 
    { 
      assert( m_Ay );
      return( *m_Ay ); 
    }

    inline const mnl::basics::Matrix& Ay() const
    { 
      assert( m_Ay );
      return( *m_Ay ); 
    }

    inline mnl::basics::Matrix& Az() 
    { 
      assert( m_Az );
      return( *m_Az ); 
    }

    inline const mnl::basics::Matrix& Az() const
    { 
      assert( m_Az );
      return( m_Az ); 
    }

    inline BC& bc() 
    { 
      return( m_b ); 
    }

    inline const BC& bc() const
    { 
      return( m_b ); 
    }

    void setOperatorX(const mnl::basics::Matrix& A, const mnl::basics::Matrix& B)
    {
      if( m_Ax )
        delete m_Ax;
      m_Ax = new mnl::basics::Matrix(A);
      diagonalizeOperator(m_Ax,m_Qx,m_eigsX,B);
    }

    void setOperatorY(const mnl::basics::Matrix& A, const mnl::basics::Matrix& B)
    {
      if( m_Ay )
        delete m_Ay;
      m_Ay = new mnl::basics::Matrix(A);
      diagonalizeOperator(m_Ay,m_Qy,m_eigsY,B);
    }

    void setOperatorZ(const mnl::basics::Matrix& A, const mnl::basics::Matrix& B)
    {
      if( m_Az )
        delete m_Az;
      m_Az = new mnl::basics::Matrix(A);
      diagonalizeOperator(m_Az,m_Qz,m_eigsZ,B);
    }

    void setGridX(const mnl::basics::Vector& gridX)
    {
      if( m_gridX )
        delete m_gridX;
      m_gridX = new mnl::basics::Vector(gridX);
    }

    void setGridY(const mnl::basics::Vector& gridY)
    {
      if( m_gridY )
        delete m_gridY;
      m_gridY = new mnl::basics::Vector(gridY);
    }

    void setGridZ(const mnl::basics::Vector& gridZ)
    {
      if( m_gridZ )
        delete m_gridZ;
      m_gridZ = new mnl::basics::Vector(gridZ);
    }

    void setWeightX(const mnl::basics::Vector& weight)
    {
      if( m_weightX )
        delete m_weightX;
      m_weightX = new mnl::basics::Vector(weight);
    }

    void setWeightY(const mnl::basics::Vector& weight)
    {
      if( m_weightY )
        delete m_weightY;
      m_weightY = new mnl::basics::Vector(weight);
    }

    void setWeightZ(const mnl::basics::Vector& weight)
    {
      if( m_weightZ )
        delete m_weightZ;
      m_weightZ = new mnl::basics::Vector(weightZ);
    }
  protected:
    BC m_b;
    mnl::basics::Vector* m_gridX;
    mnl::basics::Vector* m_gridY; 
    mnl::basics::Vector* m_gridZ;
    mnl::basics::Vector* m_weightX;
    mnl::basics::Vector* m_weightY;
    mnl::basics::Vector* m_weightZ;
    mnl::basics::Vector* m_eigsX;
    mnl::basics::Vector* m_eigsY;
    mnl::basics::Vector* m_eigsZ;
    mnl::basics::Matrix* m_Ax;
    mnl::basics::Matrix* m_Qx;
    mnl::basics::Matrix* m_Dx;
    mnl::basics::Matrix* m_Ay;
    mnl::basics::Matrix* m_Qy;
    mnl::basics::Matrix* m_Dy;
    mnl::basics::Matrix* m_Az;
    mnl::basics::Matrix* m_Qz;
    mnl::basics::Matrix* m_Dz;
    Real m_scaleX, m_scaleY, m_scaleZ;

    void diagonalizeOperator(const mnl::basics::Matrix* A, mnl::basics::Matrix*& Q, mnl::basics::Vector*& eigs, const mnl::basics::Matrix& B)
    {
      if( eigs )
        delete eigs;
      eigs = new mnl::basics::Vector("eigenvalues",B.rows());
      if( Q )
        delete Q;
      Q = new mnl::basics::Matrix("Generalized eigenvectors",B.rows(),B.rows(),false);
      A->generalizedEigenvaluesSPD(B,*Q,*eigs);

      /* obtain correct scaling of eigenvectors! */
      mnl::basics::Vector vec("Temporary",Q->rows(),false);
      for( int j=0;j<Q->cols();j++ ) {
        vec = mnl::basics::multTranspose(B,(*Q)[j],'N');
        Real scale3 = vec.dot((*Q)[j]);
        (*Q)[j] /= sqrt(scale3);
      }
    }
};
#endif

