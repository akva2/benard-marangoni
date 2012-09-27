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

#include "geometry.h"
#include "matrix.h"

using namespace std;

namespace mnl {
  namespace basics {
    void spectralElement2D::mass(Matrix& res, const Vector& weight) const
    {
      for( int j=0;j<res.cols();++j )
        for( int k=0;k<res.rows();++k )
          res[j][k] *= weight[j]*weight[k]*getGH().getJacobian()[j][k];
    }

    void spectralElement2D::invMass(Matrix& res, const Vector& weight) const
    {
      for( int j=0;j<res.cols();++j )
        for( int k=0;k<res.rows();++k )
          res[j][k] /= weight[j]*weight[k]*getGH().getJacobian()[j][k];
    }

    void spectralElement2D::evaluate(Matrix& res, const Vector& grid, 
        const Vector& ggrid, const function2D& source,
        Real t, Real Lz) const

    {
      for( int alpha=0;alpha<grid.length();++alpha ) {
        for( int beta=0;beta<grid.length();++beta ) {
          pair<Real,Real> point = getGH().evaluate(grid[beta],grid[alpha],ggrid);
          res[alpha][beta] = source.val(point.first,point.second,t);
        }
      }
    }

    void spectralElement2D::evaluate(Matrix& res, const Vector& grid, const Vector& ggrid, 
        const function3D& source, Real z, Real t) const
    {
      for( int alpha=0;alpha<grid.length();++alpha ) {
        for( int beta=0;beta<grid.length();++beta ) {
          pair<Real,Real> point = getGH().evaluate(grid[beta],grid[alpha],ggrid);
          res[alpha][beta] = source.val(point.first,point.second,z,t);
        }
      }
    }

    void spectralElement2D::evaluate(Matrices& res, const Vector& grid, const Vector& ggrid, 
        const function3D& source, Real t, Real Lz) const
    {
      for( int alpha=0;alpha<grid.length();++alpha ) {
        for( int beta=0;beta<grid.length();++beta ) {
          pair<Real,Real> point;
          if( ggrid.length() == grid.length() )
            point = make_pair(getGH().getMapping().X()[alpha][beta],getGH().getMapping().Y()[alpha][beta]);
          else
            point = getGH().evaluate(grid[beta],grid[alpha],ggrid);
          for( int l=0;l<grid.length();++l ) 
            res[l][alpha][beta] = source.val(point.first,point.second,(grid[l]+1)/2*Lz,t);
        }
      }
    }

    Real spectralElement2D::getVolume(const Vector& weight) const
    {
      Matrix* temp  = utilities::g_manager.clone(getGH().getMapping().X());
      Matrix* temp2 = utilities::g_manager.clone(getGH().getMapping().X());

      *temp = 1;
      *temp2 = 1;
      mass(*temp,weight);
      Real area = temp2->dot(*temp);

      utilities::g_manager.unlock(temp);
      utilities::g_manager.unlock(temp2);

      return( area );
    }

    pair<Real,Real> spectralElement2D::getSize(const Vector& weight, const Vector& grid) const
    {

      pair<Real,Real> result(0.f,0.f);
      /* use straight lines */
      //            pair<Real,Real> bleft  = getGH().evaluate(-1,-1,grid);
      //            pair<Real,Real> bright = getGH().evaluate( 1, -1,grid);
      //            pair<Real,Real> tright = getGH().evaluate( 1, 1,grid);
      //            pair<Real,Real> tleft  = getGH().evaluate(-1, 1,grid);
      //            result.first  = fabs((bright.first-bleft.first+tright.first-tleft.first)/2);
      //            result.second = fabs((tright.second-bright.second+tleft.second-bleft.second)/2);

      /* use curve length */
      int N = getGH().getMapping().X().rows();
      const matrixStack& D = getGH().getGeometryDerivatives();
      for( int k=0;k<N;++k) {
        result.first   +=  weight[k]*( sqrt(pow(D[0][  0][  k],2) + pow(D[2][  0][  k],2)) 
            +sqrt(pow(D[0][N-1][  k],2) + pow(D[2][N-1][  k],2)))/2;
        result.second  +=  weight[k]*( sqrt(pow(D[1][  k][  0],2) + pow(D[3][  k][  0],2)) 
            +sqrt(pow(D[1][  k][N-1],2) + pow(D[3][  k][N-1],2)))/2;
      }

      return( result );
    }

    pair<Real,Real> spectralElement2D::getAdjustedSize(const Vector& weight, const Vector& grid) const
    {
      Real area = getVolume(weight);	
      pair<Real,Real> result = getSize(weight,grid);

      //                    /* keep the width, adjust the height */
      //                    result.second 		= area/result.first;

      //                    /* adjust the width, keep the height */
      //                    result.first 	= area/result.second;

      /* adjust both */
      Real scale = sqrt(area/(result.first*result.second));
      result.first  *= scale;
      result.second *= scale;

      return( result );
    }

    void spectralElement3D::mass(Matrices& res, const Vector& weight) const
    {
      for( int l=0;l<res.matrices();++l )
        for( int j=0;j<res.cols();++j )
          for( int k=0;k<res.rows();++k )
            res[l][j][k] *= weight[l]*weight[j]*weight[k]*getGH().getJacobian()[l][j][k];
    }

    void spectralElement3D::evaluate(Matrices& res, const function3D& source, Real t) const
    {
      int N = getGH().getMapping().X().rows();
      for( int alpha=0;alpha<N;++alpha ) {
        for( int beta=0;beta<N;++beta ) {
          for( int gamma=0;gamma<N;++gamma )
            res[alpha][beta][gamma] = source.val(getGH().getMapping().X()[alpha][beta][gamma],
                getGH().getMapping().Y()[alpha][beta][gamma],
                getGH().getMapping().Z()[alpha][beta][gamma],t);
        }
      }
    }

    Real spectralElement3D::getVolume(const Vector& weight) const
    {
      Matrices* temp  = utilities::g_manager.clone(getGH().getMapping().X());
      Matrices* temp2 = utilities::g_manager.clone(getGH().getMapping().X());

      *temp  = 1;
      *temp2 = 1;
      mass(*temp,weight);
      Real volume = temp2->dot(*temp);

      utilities::g_manager.unlock(temp);
      utilities::g_manager.unlock(temp2);

      return( volume );
    }

    vector<Real> spectralElement3D::getSize(const Vector& weight,
        const Vector& grid) const
    {
      /* use curve length */
      int N = getGH().getMapping().X().rows();
      const matricesStack& D = getGH().getGeometryDerivatives();
      Real xl,yl,zl;
      vector<Real> result;
      for( int i=0;i<3;++i )
        result.push_back(0);
      for( int k=0;k<N;++k) {
        xl = (sqrt(  pow(D[0][  0][  0][k],2) /* front bottom */
              +pow(D[3][  0][  0][k],2)
              +pow(D[6][  0][  0][k],2))
            +sqrt( pow(D[0][N-1][  0][k],2) /* front top */
              +pow(D[3][N-1][  0][k],2)
              +pow(D[6][N-1][  0][k],2))
            +sqrt( pow(D[0][  0][N-1][k],2) /* back bottom */
              +pow(D[3][  0][N-1][k],2)
              +pow(D[6][  0][N-1][k],2))
            +sqrt( pow(D[0][N-1][N-1][k],2) /* back top */
              +pow(D[3][N-1][N-1][k],2)
              +pow(D[6][N-1][N-1][k],2)))/4;

        yl = (sqrt(  pow(D[1][  0][k][  0],2) /* left bottom */
              +pow(D[4][  0][k][  0],2)
              +pow(D[7][  0][k][  0],2))
            +sqrt( pow(D[1][  0][k][N-1],2) /* right bottom */
              +pow(D[4][  0][k][N-1],2)
              +pow(D[7][  0][k][N-1],2))
            +sqrt( pow(D[1][N-1][k][  0],2) /* left top */
              +pow(D[4][N-1][k][  0],2)
              +pow(D[7][N-1][k][  0],2))
            +sqrt( pow(D[1][N-1][k][N-1],2) /* right top */
              +pow(D[4][N-1][k][N-1],2)
              +pow(D[7][N-1][k][N-1],2)))/4;

        zl = (sqrt(  pow(D[2][k][  0][  0],2) /* left front */
              +pow(D[5][k][  0][  0],2)
              +pow(D[8][k][  0][  0],2))
            +sqrt( pow(D[2][k][  0][N-1],2) /* right front */
              +pow(D[5][k][  0][N-1],2)
              +pow(D[8][k][  0][N-1],2))
            +sqrt( pow(D[2][k][N-1][  0],2) /* left back */
              +pow(D[5][k][N-1][  0],2)
              +pow(D[8][k][N-1][  0],2))
            +sqrt( pow(D[2][k][N-1][N-1],2) /* right back */
              +pow(D[5][k][N-1][N-1],2)
              +pow(D[8][k][N-1][N-1],2)))/4;

        result[0]   +=  weight[k]*xl;
        result[1]   +=  weight[k]*yl;
        result[2]   +=  weight[k]*zl;
      }

      return( result );
    }

    vector<Real> spectralElement3D::getAdjustedSize(const Vector& weight, const Vector& grid) const
    {
      Real volume = getVolume(weight);
      vector<Real> result = getSize(weight,grid);

      /* adjust all */
      Real scale = pow(volume/(result[0]*result[1]*result[2]),Real(1)/3);
      result[0] *= scale;
      result[1] *= scale;
      result[2] *= scale;

      return( result );
    }

    geometryStack::geometryStack(const Vector& grid, const Vector& weight) :
      m_dotter(*this), m_dotterP(*this,true), m_ggrid(grid), m_weight(weight), m_mine(true)
    {
      m_mult = NULL;
      m_mass = NULL;
      m_imass = NULL;
      m_Lz = 2;
    }

    geometryStack::~geometryStack()
    {
      if( m_mine )
        for( unsigned int i=0;i<m_grid.size();++i )
          delete m_grid[i];

      utilities::g_manager.unlock(m_mult);
      utilities::g_manager.unlock(m_mass);
      utilities::g_manager.unlock(m_imass);
    }

    void geometryStack::mass(matrixStack& res, bool dosum) const
    {
      int max=res.size();
#pragma omp parallel for schedule(static)
      for( int i=0;i<max;++i )
        basics::multPointwise(res[i],(*m_mass)[i]);

      if( dosum )
        dssum(res);
    }

    void geometryStack::mass(matrixStack& res, const vector<int>& desc) const
    {
      int max = desc.size();
#pragma omp parallel for schedule(static)
      for( int i=0;i<max;++i )
        basics::multPointwise(res[i],(*m_mass)[desc[i]]);
    }

    void geometryStack::mass(matricesStack& res, const vector<int>& desc,
        int start) const
    {
      for( int l=0;l<res[0].matrices();++l ) {
        mass(res.at(l),desc);
        res.at(l) *= m_weight[l+start]*m_Lz/2;
      }
    }

    void geometryStack::mass(matricesStack& u, bool dosum) const
    {
      for( int l=0;l<u[0].matrices();++l ) {
        mass(u.at(l),dosum);
        u.at(l) *= m_weight[l]*m_Lz/2;
      }
    }

    void geometryStack::invMass(matrixStack& res) const // enforces dssum through operator
    {
      int max=res.size();
#pragma omp parallel for schedule(static)
      for( int i=0;i<max;++i )
        basics::multPointwise(res[i],(*m_imass)[i]);
    }

    void geometryStack::invMassP(matrixStack& res) const // enforces dssum through operator
    {
      int max=res.size();
#pragma omp parallel for schedule(static)
      for( int i=0;i<max;++i )
        basics::multPointwise(res[i],(*m_imassP)[i]);
    }

    void geometryStack::invMass(matrixStack& res, const vector<int>& r) const // enforces dssum through operator
    {
      int max=res.size();
#pragma omp parallel for schedule(static)
      for( int i=0;i<max;++i )
        basics::multPointwise(res[i],(*m_imass)[r[i]]);
    }

    void geometryStack::invMassP(matrixStack& res, const vector<int>& r) const // enforces dssum through operator
    {
      int max=res.size();
#pragma omp parallel for schedule(static)
      for( int i=0;i<max;++i )
        basics::multPointwise(res[i],(*m_imassP)[r[i]]);
    }

    void geometryStack::invMass(matricesStack& res) const
    {
      for( int i=0;i<res[0].matrices();++i ) {
        invMass(res.at(i));
        res.at(i) *= Real(1)/(m_weight[i]*m_Lz/2);
      }
    }

    void geometryStack::invMassP(matricesStack& res) const
    {
      for( int i=0;i<res[0].matrices();++i ) {
        invMassP(res.at(i));
        res.at(i) *= Real(1)/(m_weight[i]*m_Lz/2);
      }
    }

    void geometryStack::invMass(matricesStack& res, const std::vector<int>& r) const
    {
      for( int l=0;l<res[0].matrices();++l ) {
        invMass(res.at(l),r);
        res.at(l) *= Real(1)/(m_weight[l]*m_Lz/2);
      }
    }

    void geometryStack::compute(Real t)
    {
      for( int i=0;i<m_grid.size();++i )
        m_grid[i]->getGH().compute(t);
    }

    void geometryStack::computeMultiplicities()
    {
      if( !m_mult ) {
        m_mult = utilities::g_manager.aquireMatrixStack("multiplicity",*this);
        m_multP = utilities::g_manager.aquireMatrixStack("periodic multiplicity",*this);
      }
      matrixStack& mult = *m_mult;
      for( unsigned int i=0;i<mult.size();++i ) {
        for( int j=0;j<mult[i].cols();++j ) {
          for( int k=0;k<mult[i].rows();++k )
            mult[i][j][k] = Real(1);
        }
      }
      *m_multP = mult;
      periodicDssum(*m_multP);
      dssum(mult);
      for( unsigned int i=0;i<mult.size();++i ) {
        for( int j=0;j<mult[i].cols();++j ) {
          for( int k=0;k<mult[i].rows();++k ) {
            mult[i][j][k] = Real(1)/mult[i][j][k];
            (*m_multP)[i][j][k] = Real(1)/(*m_multP)[i][j][k];
          }
        }
      }
    }

    void geometryStack::computeMass()
    {
      if( !m_mass) {
        m_mass   = utilities::g_manager.aquireMatrixStack("mass",*this);
        m_imassP = utilities::g_manager.aquireMatrixStack("mass",*this);
        m_imass  = utilities::g_manager.aquireMatrixStack("inverse mass",*this);
      }
      basics::matrixStack& mass  = *m_mass;
      basics::matrixStack& imass = *m_imass;
      int max=mass.size();
#pragma omp parallel for schedule(static)
      for( int i=0;i<max;++i ) {
        mass[i] = 1;
        m_grid[i]->mass(mass[i],m_weight);
      }

      imass = mass;
      *m_imassP = mass;
      dssum(imass);
      periodicDssum(*m_imassP);
#pragma omp parallel for schedule(static)
      for( int i=0;i<max;++i )
        for( int j=0;j<mass[i].cols();++j )
          for( int k=0;k<mass[i].rows();++k ) {
            imass[i][j][k] = Real(1)/imass[i][j][k];
            (*m_imassP)[i][j][k] = Real(1)/(*m_imassP)[i][j][k];
          }
    }

    Real geometryStack::dot(const matrixStack& in1, const matrixStack& in2, const matrixStack* mult) const
    {
      Real result=0;
      if( mult == NULL )
        mult = m_mult;

      int max=in1.size();
#pragma omp parallel for reduction(+:result) schedule(static)
      for( int i=0;i<max;++i ) {
        for( int j=0;j<(*mult)[i].cols();++j )
          for( int k=0;k<(*mult)[i].rows();++k )
            result += in1[i][j][k]*in2[i][j][k]*(*mult)[i][j][k];
      }

      return( result );
    }

    Real geometryStack::dot(const matricesStack& in1, const matricesStack& in2) const
    {
      Real result=0;
      const matrixStack& mult = *m_mult;
      int max=in1.size();
#pragma omp parallel for reduction(+:result) schedule(static)
      for( int i=0;i<max;++i ) {
        for( int l=0;l<in1[i].matrices();++l )
          for( int j=0;j<mult[i].cols();++j )
            for( int k=0;k<mult[i].rows();++k )
              result += in1[i][l][j][k]*in2[i][l][j][k]*mult[i][j][k];
      }

      return( result );
    }

    Real geometryStack::dot(const Matrices& in1, const Matrices& in2, int i) const
    {
      Real result=0;
      const matrixStack& mult = *m_mult;
      int max=in1.matrices();
#pragma omp parallel for reduction(+:result) schedule(static)
      for( int l=0;l<max;++l )
        for( int j=0;j<mult[i].cols();++j )
          for( int k=0;k<mult[i].rows();++k )
            result += in1[l][j][k]*in2[l][j][k]*mult[i][j][k];

      return( result );
    }

    Real geometryStack::dot(const Matrix& in1, const Matrix& in2, int i, const Matrix* mult) const
    {
      Real result=0;
      if( !mult )
        mult = &((*m_mult)[i]);
      for( int j=0;j<mult->cols();++j )
        for( int k=0;k<mult->rows();++k )
          result += in1[j][k]*in2[j][k]*(*mult)[j][k];

      return( result );
    }

    vector<matrixStack* > geometryStack::getInterpolatedGeometryDerivatives(const Matrix& GLL2G) const
    {
      /* setup interpolated geometry derivatives */
      vector<matrixStack*> result;
      for( int n=0;n<4;++n ) {
        stringstream str;
        str << "interpolated geometry derivatives ";
        str << n;
        result.push_back(new matrixStack(str.str(),GLL2G.rows(),
              GLL2G.rows(),size()));
        for( int i=0;i<size();++i )
          interpolate((*result[n])[i],
              m_grid[i]->getGH().getGeometryDerivatives()[n],
              GLL2G);
      }

      return( result );
    }

    void geometryStack::localToGlobal(basics::Vector& result,
        const basics::matrixStack& input,
        const basics::matrixStack& LG) const
    {
      for( int n=0;n<input.size();++n )
        for( int j=0;j<input[n].cols();++j )
          for( int k=0;k<input[n].rows();++k ) {
            if( LG[n][j][k] > -1 )
              result[LG[n][j][k]] = input[n][j][k];
          }
    }

    void geometryStack::globalToLocal(basics::matrixStack& result,
        const basics::Vector& input,
        const basics::matrixStack& LG) const
    {
      for( int n=0;n<result.size();++n )
        for( int j=0;j<result[n].cols();++j )
          for( int k=0;k<result[n].rows();++k )
            if( LG[n][j][k] > -1 )
              result[n][j][k] = input[LG[n][j][k]];
            else
              result[n][j][k] = 0;
    }

    void geometryStack::localToGlobal(basics::Vector& result,
        const basics::matricesStack& input,
        const basics::matricesStack& LG) const
    {
      for( int n=0;n<input.size();++n )
        for( int l=0;l<input[n].matrices();++l )
          for( int j=0;j<input[n].cols();++j )
            for( int k=0;k<input[n].rows();++k ) {
              if( LG[n][l][j][k] > -1 )
                result[LG[n][l][j][k]] = input[n][l][j][k];
            }
    }

    void geometryStack::globalToLocal(basics::matricesStack& result,
        const basics::Vector& input,
        const basics::matricesStack& LG) const
    {
      for( int n=0;n<result.size();++n )
        for( int l=0;l<result[n].matrices();++l )
          for( int j=0;j<result[n].cols();++j )
            for( int k=0;k<result[n].rows();++k )
              if( LG[n][l][j][k] > -1 )
                result[n][l][j][k] = input[LG[n][l][j][k]];
              else
                result[n][l][j][k] = 0;
    }

    geometryStack3D::geometryStack3D(const basics::Vector& weight) :
      m_dotter(*this), m_weight(weight), m_mult(NULL), m_mass(NULL), m_imass(NULL)
    {
    }

    geometryStack3D::~geometryStack3D()
    {
      for( unsigned int i=0;i<m_grid.size();++i )
        delete m_grid[i];
      for( unsigned int i=0;i<m_referenceDerivatives.size();++i )
        delete m_referenceDerivatives[i];

      utilities::g_manager.unlock(m_mult);
      utilities::g_manager.unlock(m_mass);
      utilities::g_manager.unlock(m_imass);
    }

    Real geometryStack3D::dot(const matricesStack& in1, const matricesStack& in2) const
    {
      Real result=0.f;
      int max=in1.size();
#pragma omp parallel for reduction(+:result) schedule(static)
      for( int i=0;i<max;++i )
        for( int l=0;l<in1[i].matrices();++l ) 
          for( int j=0;j<in1[i].cols();++j )
            for( int k=0;k<in1[i].rows();++k )
              result += in1[i][l][j][k]*in2[i][l][j][k]*(*m_mult)[i][l][j][k];

      return( result );
    }

    Real geometryStack3D::dot(const Matrices& in1, const Matrices& in2, int i) const
    {
      Real result=0;
      int max=in1.matrices();
#pragma omp parallel for reduction(+:result) schedule(static)
      for( int l=0;l<max;++l )
        for( int j=0;j<in1.cols();++j )
          for( int k=0;k<in1.rows();++k )
            result += in1[l][j][k]*in2[l][j][k]*(*m_mult)[i][l][j][k];

      return( result );
    }

    void geometryStack3D::computeMultiplicities()
    {
      int N = m_grid[0]->getGH().getMapping().X().rows();
      if( !m_mult )
        m_mult = utilities::g_manager.aquireMatricesStack("multiplicity",N,N,N,m_grid.size());
      matricesStack& mult = *m_mult;
      mult = 1;
      dssum(mult);
      int max = mult.size();
#pragma omp parallel for schedule(static)
      for( int i=0;i<max;++i )
        for( int l=0;l<mult[i].matrices();++l )
          for( int j=0;j<mult[i].cols();++j )
            for( int k=0;k<mult[i].rows();++k )
              mult[i][l][j][k] = Real(1)/mult[i][l][j][k];
    }

    void geometryStack3D::computeMass()
    {
      int N = m_grid[0]->getGH().getMapping().X().rows();
      if( !m_mass) {
        m_mass  = utilities::g_manager.aquireMatricesStack("mass",N,N,N,m_grid.size());
        m_imass = utilities::g_manager.aquireMatricesStack("inverse mass",N,N,N,m_grid.size());
      }
      basics::matricesStack& mass  = *m_mass;
      basics::matricesStack& imass = *m_imass;
      int max=mass.size();
#pragma omp parallel for schedule(static)
      for( int i=0;i<max;++i ) {
        mass[i] = 1;
        m_grid[i]->mass(mass[i],m_weight);
      }
      imass = mass;
      dssum(imass);
#pragma omp parallel for schedule(static)
      for( int i=0;i<max;++i )
        for( int l=0;l<mass[i].matrices();++l )
          for( int j=0;j<mass[i].cols();++j )
            for( int k=0;k<mass[i].rows();++k )
              imass[i][l][j][k] = Real(1)/imass[i][l][j][k];
    }

    const vector<matricesStack*>& geometryStack3D::getReferenceGeometryDerivatives() const
    {
      return m_referenceDerivatives;
    }

    vector<matricesStack*>& geometryStack3D::getReferenceGeometryDerivatives()
    {
      return m_referenceDerivatives;
    }

    void geometryStack3D::computeReferenceGeometryDerivatives()
    {
      vector<matricesStack*>& result = m_referenceDerivatives;
      for( int i=0;i<9;++i )
        result.push_back(new matricesStack("refderivative",
              (*m_mass)[0].rows(),
              (*m_mass)[0].cols(),
              (*m_mass)[0].matrices(),
              m_mass->size()));
      for( int i=0;i<size();++i) {
        const basics::matricesStack& GD = m_grid[i]->getGH().getGeometryDerivatives();
        const basics::Matrices& J = m_grid[i]->getGH().getJacobian();
        for( int l=0;l<J.matrices();++l ) {
          for( int j=0;j<J.cols();++j ) {
            for( int k=0;k<J.rows();++k ) {
              (*result[0])[i][l][j][k] = GD[4][l][j][k]*GD[8][l][j][k]-GD[5][l][j][k]*GD[7][l][j][k];
              (*result[1])[i][l][j][k] = GD[2][l][j][k]*GD[7][l][j][k]-GD[1][l][j][k]*GD[8][l][j][k];
              (*result[2])[i][l][j][k] = GD[1][l][j][k]*GD[5][l][j][k]-GD[2][l][j][k]*GD[4][l][j][k];
              (*result[3])[i][l][j][k] = GD[5][l][j][k]*GD[6][l][j][k]-GD[3][l][j][k]*GD[8][l][j][k];
              (*result[4])[i][l][j][k] = GD[0][l][j][k]*GD[8][l][j][k]-GD[2][l][j][k]*GD[6][l][j][k];
              (*result[5])[i][l][j][k] = GD[2][l][j][k]*GD[3][l][j][k]-GD[0][l][j][k]*GD[5][l][j][k];
              (*result[6])[i][l][j][k] = GD[3][l][j][k]*GD[7][l][j][k]-GD[4][l][j][k]*GD[6][l][j][k];
              (*result[7])[i][l][j][k] = GD[1][l][j][k]*GD[6][l][j][k]-GD[0][l][j][k]*GD[7][l][j][k];
              (*result[8])[i][l][j][k] = GD[0][l][j][k]*GD[4][l][j][k]-GD[1][l][j][k]*GD[3][l][j][k];
            }
          }
        }
      }
    }

    vector<matricesStack*> geometryStack3D::getInterpolatedGeometryDerivatives(const Matrix& GLL2G) const
    {
      vector<matricesStack*> result;
      for( int n=0;n<9;++n ) {
        stringstream str;
        str << "interpolated geometry derivatives ";
        str << n;
        result.push_back(new matricesStack(str.str(),GLL2G.rows(),
              GLL2G.rows(),GLL2G.rows(),size()));
        for( int i=0;i<size();++i )
          interpolate((*result[n])[i],
              m_grid[i]->getGH().getGeometryDerivatives()[n],
              GLL2G);
      }

      return( result );
    }

    vector<matricesStack*> geometryStack3D::getInterpolatedReferenceGeometryDerivatives(const Matrix& GLL2G) const
    {
      vector<matricesStack*> reference = getReferenceGeometryDerivatives();
      vector<matricesStack*> result;
      for( int n=0;n<9;++n ) {
        stringstream str;
        str << "interpolated reference geometry derivatives ";
        str << n;
        result.push_back(new matricesStack(str.str(),GLL2G.rows(),
              GLL2G.rows(),GLL2G.rows(),size()));
        for( int i=0;i<size();++i )
          interpolate((*result[n])[i],
              (*reference[n])[i],
              GLL2G);
      }

      return( result );
    }

    void geometryStack::interpolate(Matrix& result, 
        const Matrix& G,
        const Matrix& GLL2G)
    {
      Matrix temp ("temp",GLL2G.rows(),GLL2G.cols());
      multTranspose(temp,GLL2G,G,'N','N');
      multTranspose(result,temp,GLL2G,'N','T');
    }

    void geometryStack3D::interpolate(Matrices& result, 
        const Matrices& G,
        const Matrix& GLL2G)
    {
      Matrices temp2("temp",GLL2G.rows(),GLL2G.rows(),G.matrices());

      for( int l=0;l<G.matrices();++l )
        geometryStack::interpolate(temp2[l],G[l],GLL2G);
      applyLocalGlobal(result,temp2,GLL2G,'N','T');
    }

    void geometryStack3D::invMass(matricesStack& res) const
    {
      int max=size();
#pragma omp parallel for schedule(static)
      for( int i=0;i<max;++i )
        multPointwise(res[i],(*m_imass)[i]);
    }

    void geometryStack3D::invMass(matricesStack& res, const std::vector<int>& r) const
    {
      int max=r.size();
#pragma omp parallel for schedule(static)
      for( int i=0;i<max;++i )
        multPointwise(res[i],(*m_imass)[r[i]]);
    }

    void geometryStack3D::mass(matricesStack& res) const
    {
      int max=size();
#pragma omp parallel for schedule(static)
      for( int i=0;i<max;++i )
        multPointwise(res[i],(*m_mass)[i]);
    }

    void geometryStack3D::mass(matricesStack& res, const std::vector<int>& r) const
    {
      int max=r.size();
#pragma omp parallel for schedule(static)
      for( int i=0;i<max;++i )
        multPointwise(res[i],(*m_mass)[r[i]]);
    }


    vector<basics::matricesStack*>
      geometryStack3D::setupFakeStacks(basics::matricesStack& op, int level, int perlevel) const
      {
        if( level < 0 )
          level = m_level;
        if( perlevel < 0 )
          perlevel = m_perlevel;
        vector<basics::matricesStack*> result;
        if( level == 1 )
          result.push_back(&op);
        else {
          for( int i=0;i<level;++i ) {
            std::vector<basics::Matrices*> mats;
            for( int j=0;j<perlevel;++j )
              mats.push_back(&op[i*perlevel+j]);
            result.push_back(new basics::matricesStack("fake matrices stack",mats));
          }
        }

        return( result );
      }

    void geometryStack3D::killFakeStacks(vector<basics::matricesStack*> op) const
    {
      if( op.size() == 1 )
        return;

      for( int i=0;i<op.size();++i )
        delete op[i];
    }

    basics::Field3<basics::matricesStack>*
      geometryStack3D::localToGlobalZ(const basics::Field3<basics::matricesStack>& input,
          int rank, int size, bool L2, bool dosum) const
      {
        int N = input[0][0].matrices();
        int planes=N*m_level;
        if( !L2 )
          planes -= m_level-1;

        basics::Field3<basics::matricesStack>* result 
          = utilities::g_manager.aquireMatricesStackField("temp",input[0][0].rows(),
              input[0][0].cols(),planes,
              m_perlevel/size);
        for( int i=0;i<3;++i )
          localToGlobalZ(input[i],rank,size,L2,dosum,&(*result)[i]);

        return result;
      }

    basics::matricesStack* 
      geometryStack3D::localToGlobalZ(const basics::matricesStack& input,
          int rank, int size, bool L2, bool dosum,
          basics::matricesStack* result) const
      {
        std::vector<basics::matricesStack*> ops = setupFakeStacks(const_cast<basics::matricesStack&>(input),m_level,m_perlevel/size);
        int N = input[0].matrices();
        int planes=N*m_level;
        if( !L2 )
          planes -= m_level-1;

        if( !result )
          result = 
            utilities::g_manager.aquireMatricesStack("local in Z",input[0].rows(),
                input[0].cols(),planes,m_perlevel/size);
        int startplane = 0;
        if( dosum ) {
          result->clear();
          for( int i=0;i<m_level;++i ) {
            for( int l=0;l<N;++l )
              result->at(startplane+l) += (*ops[i]).at(l);
            startplane += N-(L2?0:1);
          }
        } else {
          for( int i=0;i<m_level;++i ) {
            for( int l=0;l<N;++l )
              result->at(startplane+l) = (*ops[i]).at(l);
            startplane += N-(L2?0:1);
          }
        }

        killFakeStacks(ops);

        return( result );
      }

    void geometryStack3D::globalToLocalZ(basics::Field3<basics::matricesStack>& result,
        basics::Field3<basics::matricesStack>* input,
        int rank, int size, bool L2) const
    {
      for( int i=0;i<3;++i )
        globalToLocalZ(result[i],&(*input)[i],rank,size,L2);

      delete input;
    }

    void geometryStack3D::globalToLocalZ(basics::matricesStack& result,
        basics::matricesStack* input,
        int rank, int size, bool L2) const
    {
      std::vector<basics::matricesStack*> ops = setupFakeStacks(result,m_level,m_perlevel/size);
      int N = result[0].matrices();
      int startplane = 0;
      for( int i=0;i<m_level;++i ) {
        for( int l=0;l<N;++l )
          ops[i]->at(l) = input->at(startplane+l);
        startplane += N-(L2?0:1);
      }

      killFakeStacks(ops);
      utilities::g_manager.unlock(input);
    }

  }
}

