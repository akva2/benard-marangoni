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

#include "legendrelegendrew.h"

#include "mnl/memtracker.h"
#include "mnl/grid.h"

#include "poissonsolver-lgw.h"

#include <sstream>

using namespace mnl;
namespace legendreLegendreW {

  /* 2D */
  void laplacianEvaluator2D::evaluate(basics::Matrix& mat, const basics::Matrix& mat2) const
  {
    /* A_x */
    basics::multTranspose(mat,m_Ax,mat2,'N','N',m_nu);
    massY(mat,m_weightY);

    basics::Matrix* T = utilities::g_manager.clone(mat);

    /* A_y */
    basics::multTranspose(*T,mat2,m_Ay,'N','T',m_nu); 
    massX(*T,m_weightX);
    mat += *T;

    if( m_mu ) {
      *T = mat2;
      massX(*T,m_weightX);
      massY(*T,m_weightY);
      mat.axpy(m_mu,*T);
    }

    utilities::g_manager.unlock(T);
  }

  basics::Matrix reconstruct(const basics::Matrix& u, const int n_Nx, const int n_Ny, const basics::Vector& gridX, const basics::Vector& gridY)
  {
    basics::Vector gridX2("grid twice the size",n_Nx);
    basics::Vector weightX2("weight twice the size",n_Nx);
    utilities::GLL::GaussLobattoLegendreGrid(gridX2);
    utilities::GLL::GaussLobattoLegendreWeights(weightX2,gridX2);
    gridX2 += Real(1);
    gridX2 *= M_PI;
    basics::Matrix interpolX = utilities::GLL::interpolationMatrix(gridX,gridX2);

    basics::Vector gridY2("grid twice the size",n_Ny);
    basics::Vector weightY2("weight twice the size",n_Ny);
    utilities::GLL::GaussLobattoLegendreGrid(gridY2);
    utilities::GLL::GaussLobattoLegendreWeights(weightY2,gridY2);
    basics::Matrix interpolY = utilities::GLL::interpolationMatrix(gridY,gridY2);

    basics::Matrix V2("reconstructed fully",n_Nx,n_Ny);
    basics::Matrix temp("temp",n_Nx,u.cols());
    basics::multTranspose(temp,interpolX,u,'N','N');
    basics::multTranspose(V2,temp,interpolY,'N','T');

    return V2;
  }

  Real L2norm(const basics::Matrix& u, const basics::Vector& weightX, const basics::Vector& weightY)
  {
    basics::Matrix temp(u);
    massX(temp,weightX);
    massY(temp,weightY);

    Real L2U2 = u.dot(temp);

    return( sqrt(L2U2) );
  }

  void removeHydrostaticMode(basics::Matrix& u, const basics::Vector& weightX, const basics::Vector& weightY)
  {
    basics::Matrix* e2 = utilities::g_manager.clone(u);

    *e2 = u;
    massX(*e2,weightX);
    massY(*e2,weightY);
    Real foo = e2->sum()/(Real(4)*M_PI);
    u -= foo; 

    utilities::g_manager.unlock(e2);
  }

  Real convergence(const basics::Matrix& u, const basics::Vector& gridX, const basics::Vector& gridY, const basics::function2D& exactFunc, const Real t, bool relative)
  {
    basics::Matrix u2 = reconstruct(u,2*u.rows()+1,2*u.cols()+1,gridX,gridY);

    basics::Vector gridX2("grid twice the size",u2.rows());
    basics::Vector weightX2("weight twice the size",u2.rows());
    utilities::GLL::GaussLobattoLegendreGrid(gridX2);
    utilities::GLL::GaussLobattoLegendreWeights(weightX2,gridX2);
    gridX2 += Real(1);
    gridX2 *= M_PI;

    basics::Vector gridY2("grid twice the size",u2.cols());
    basics::Vector weightY2("weight twice the size",u2.cols());
    utilities::GLL::GaussLobattoLegendreGrid(gridY2);
    utilities::GLL::GaussLobattoLegendreWeights(weightY2,gridY2);

    removeHydrostaticMode(u2,weightX2,weightY2);
    if( relative && exactFunc.id() > -1 ) {
      basics::Matrix exact("exact",u2.rows(),u2.cols());
      utilities::evaluate(exact,exactFunc,gridX2,gridY2,t);
      removeHydrostaticMode(exact,weightX2,weightY2);
      u2 -= exact;

      return( L2norm(u2,weightX2,weightY2)/L2norm(exact,weightX2,weightY2) );
    } else {
      return( L2norm(u2,weightX2,weightY2) );
    }
  }

  Real convergenceH1(const basics::Matrix& u, const basics::Vector& gridX, const basics::Vector& gridY, const basics::function2D& exactFunc, const Real t)
  {
    basics::Matrix u2 = reconstruct(u,2*u.rows()+1,2*u.cols()+1,gridX,gridY);

    basics::Vector gridX2("grid twice the size",u2.rows());
    basics::Vector weightX2("weight twice the size",u2.rows());
    utilities::GLL::GaussLobattoLegendreGrid(gridX2);
    utilities::GLL::GaussLobattoLegendreWeights(weightX2,gridX2);
    gridX2 += Real(1);
    gridX2 *= M_PI;

    basics::Vector gridY2("grid twice the size",u2.cols());
    basics::Vector weightY2("weight twice the size",u2.cols());
    utilities::GLL::GaussLobattoLegendreGrid(gridY2);
    utilities::GLL::GaussLobattoLegendreWeights(weightY2,gridY2);
    //    gridY2 += 1;

    basics::Matrix exact("exact",u2.rows(),u2.cols());
    utilities::evaluate(exact,exactFunc,gridX2,gridY2,t);

    poissonSolver solver(exact.rows()-1,exact.cols()-1,0,poissonSolver::Nonhomogenous);
    laplacianEvaluator2D laplace(solver.Ax(),solver.Ay(),weightX2,weightY2);
    u2 -= exact;

    basics::Matrix foo(u2);

    laplace.evaluate(foo,u2);
    Real norm1 = u2.dot(foo);
    laplace.evaluate(foo,exact);
    Real norm2 = exact.dot(foo);
    if( exactFunc.id() < 0 )
      norm2 = 1;

    return sqrt(norm1/norm2);
  }

  void axpyPreserveBorders(basics::Matrix& result, const basics::Matrix& phi, const Real factor)
  {
    int N=result.rows()-2;
    for( int j=1;j<result.cols()-1;++j )
      BLASRPFX(axpy,BLASINT(N),BLASREAL(factor),BLASPREAL(phi.data()[j]+1),BLASINT(mnlIntOne),BLASPREAL(result.data()[j]+1),BLASINT(mnlIntOne));
  }

  void scalePreserveBorders(basics::Matrix& result, const Real factor)
  {
    int N=result.rows()-2;
    for( int j=1;j<result.cols()-1;++j )
      BLASRPFX(scal,BLASINT(N),BLASREAL(factor),BLASPREAL(result.data()[j]+1),BLASINT(mnlIntOne));
  }

  void copyPreserveBorders(basics::Matrix& result, const basics::Matrix& phi)
  {
    int N=result.rows()-2;
    for( int j=1;j<result.cols()-1;++j )
      BLASRPFX(copy,BLASINT(N),BLASPREAL(phi.data()[j]+1),BLASINT(mnlIntOne),BLASPREAL(result.data()[j]+1),BLASINT(mnlIntOne));
  }

  void diffx(basics::Matrix& result, const basics::Matrix& D, Real scale)
  {
    LOG("find x derivative");

    basics::Matrix* buffer = utilities::g_manager.aquireMatrix("temp",result.rows(),result.cols());
    *buffer = result;

    basics::multTranspose(result,D,*buffer,'N','N',scale/M_PI);

    utilities::g_manager.unlock(buffer);
  }

  void diffx(basics::Matrix& result, const basics::Matrix& D, const basics::Matrix& input, const Real scaleA, const Real scaleB)
  {
    LOG("find x derivative");

    basics::multTranspose(result,D,input,'N','N',scaleA/M_PI,scaleB);
  }

  void diffy(basics::Matrix& result, const basics::Matrix& D, const Real scaleInput)
  {
    assert( result.cols() == D.cols() );
    LOG("find y derivative");

    basics::Matrix* T = utilities::g_manager.aquireMatrix("buffer", result.rows(),result.cols());

    *T = result; // cannot do this inplace due to the transpose
    basics::multTranspose(result,*T,D,'N','T',scaleInput);

    utilities::g_manager.unlock(T);
  }

  void diffy(basics::Matrix& result, const basics::Matrix& D, const basics::Matrix& input, const Real scaleOutput, const Real scaleInput)
  {
    assert( result.cols() == D.cols() && result.cols() == input.cols() && result.rows() == input.rows() );
    LOG("find y derivative");

    basics::multTranspose(result,input,D,'N','T',scaleInput,scaleOutput);
  }

  void interpolate(basics::Matrix& result, const basics::Matrix& Ix, const basics::Matrix& Iy, const basics::Matrix& U, Real alpha, Real beta)
  {
    basics::Matrix* temp = utilities::g_manager.aquireMatrix("temp",Ix.rows(),U.cols());

    interpolateX(*temp,Ix,U);
    interpolateY(result,Iy,*temp,alpha,beta);

    utilities::g_manager.unlock(temp);
  }

  void interpolateX(basics::Matrix& result, const basics::Matrix& I, const basics::Matrix& U, Real alpha, Real beta)
  {
    basics::multTranspose(result,I,U,'N','N',alpha,beta);
  }

  void interpolateY(basics::Matrix& result, const basics::Matrix& I, const basics::Matrix& U, Real alpha, Real beta)
  {
    assert( result.cols() == I.rows() && U.cols() == I.cols() );

    basics::multTranspose(result,U,I,'N','T',alpha,beta);
  }

  void massX(basics::Matrix& result, const basics::Vector& weight, Real scale)
  {
    for (int k=0;k<result.rows();++k )
      result.scaleRow(k,scale*weight[k]);
  }

  void invMassX(basics::Matrix& result, const basics::Vector& weight, Real scale)
  {
    for (int k=0;k<result.rows();++k )
      result.scaleRow(k,Real(1)/(scale*weight[k]));
  }

  void massY(basics::Matrix& result, const basics::Vector& weight, Real scale)
  {
    for( int j=0;j<result.cols();++j )
      result[j] *= scale*weight[j];
  }

  void invMassY(basics::Matrix& result, const basics::Vector& weight, Real scale)
  {
    for( int j=0;j<result.cols();++j )
      result[j] *= Real(1)/(scale*weight[j]);
  }

  std::string errorReport(const basics::Matrix& p, const basics::Field2D& u, const basics::Vector& gridX, const basics::Vector& gridY, const basics::function2D& pExact, const basics::function2D& uxExact, const basics::function2D& uyExact, const Real t)
  {
    std::stringstream stream;

    /* pressure, L^2 */
    Real pError = convergence(p,gridX,gridY,pExact,t);
    stream << "Error " << p.name() << ": " << pError << " (L^2)" << std::endl;

    /* velocity, H^1 */
    Real xError = convergenceH1(u.X(),gridX,gridY,uxExact,t);
    Real yError = convergenceH1(u.Y(),gridX,gridY,uyExact,t);
    stream << "Error " << u.X().name() << ": " << xError << std::endl;
    stream << "Error " << u.Y().name() << ": " << yError << std::endl;
    stream << "Mean error " << sqrt(xError*xError+yError*yError) << " (H^1)" << std::endl;

    /* velocity, L^2 */
    xError = convergence(u.X(),gridX,gridY,uxExact,t);
    yError = convergence(u.Y(),gridX,gridY,uyExact,t);
    stream << "Error " << u.X().name() << ": " << xError << std::endl;
    stream << "Error " << u.Y().name() << ": " << yError << std::endl;
    stream << "Mean error " << sqrt(xError*xError+yError*yError) << " (L^2)" << std::endl;

    /* divergence, L^2 */
    basics::Matrix* buffer = utilities::g_manager.clone(u.X());
    basics::Matrix* Dx = utilities::g_manager.aquireMatrix("temp",u.rows(),u.rows()); 
    basics::Matrix* Dy = utilities::g_manager.aquireMatrix("temp",u.cols(),u.cols()); 

    utilities::GLL::LagrangeDerivativeMatrix(*Dx);
    utilities::GLL::LagrangeDerivativeMatrix(*Dy);

    diffx(*buffer,*Dx,u.X());
    diffy(*buffer,*Dy,u.Y(),Real(1),Real(1));

    Real udivError = convergence(*buffer,gridX,gridY,pExact,t,false);
    stream << "Error (divergence): " << udivError << " (L^2)" << std::endl;

    /* vorticity, L^2 */
    diffx(*buffer,*Dx,u.Y());
    basics::multTranspose(*buffer,u.X(),*Dy,'N','T',Real(-1),Real(1));

    basics::Vorticity2D vorticity(uxExact,uyExact);
    Real vortError = convergence(*buffer,gridX,gridY,vorticity,t);
    stream << "Error (vorticity): " << vortError << " (L^2)" << std::endl;

    utilities::g_manager.unlock(buffer);
    utilities::g_manager.unlock(Dx);
    utilities::g_manager.unlock(Dy);

    return stream.str();
  }

  std::string errorReport(const basics::Field2D& u, const basics::Vector& gridX, const basics::Vector& gridY, const basics::function2D& uxExact, const basics::function2D& uyExact, const Real t)
  {
    std::stringstream stream;

    Real xError = convergenceH1(u.X(),gridX,gridY,uxExact,t);
    Real yError = convergenceH1(u.Y(),gridX,gridY,uyExact,t);
    stream << "Error " << u.X().name() << ": " << xError << std::endl;
    stream << "Error " << u.Y().name() << ": " << yError << std::endl;
    stream << "Mean error " << sqrt(xError*xError+yError*yError) << " (H^1)" << std::endl;

    return stream.str();
  }

  std::string errorReport(const basics::Matrix& u, const basics::Vector& gridX, const basics::Vector& gridY, const basics::function2D& uExact, const Real t)
  {
    std::stringstream stream;

    Real xError = convergenceH1(u,gridX,gridY,uExact,t);
    stream << "Error " << u.name() << ": " << xError << " (H^1)" << std::endl;

    return stream.str();
  }

  void applyHomogenousDirichletBC(basics::Matrix& u)
  {
    u.clearRow(0);
    for (int k=0;k<u.rows();++k )
      u[0][k] = u[u.cols()-1][k] = 0.f;
    u.clearRow(u.rows()-1);
  }

  /* 3D */
  void massX(basics::Matrices& result, const basics::Vector& weight, Real scale)
  {
    for (int l=0;l<result.matrices();++l )
      massX(result[l],weight,scale);
  }

  void invMassX(basics::Matrices& result, const basics::Vector& weight, Real scale)
  {
    for (int l=0;l<result.matrices();++l )
      invMassX(result[l],weight,scale);
  }

  void massY(basics::Matrices& result, const basics::Vector& weight, Real scale)
  {
    for (int l=0;l<result.matrices();++l )
      massY(result[l],weight,scale);
  }

  void invMassY(basics::Matrices& result, const basics::Vector& weight, Real scale)
  {
    for (int l=0;l<result.matrices();++l )
      invMassY(result[l],weight,scale);
  }

  void massZ(basics::Matrices& result, const basics::Vector& weight, Real scale)
  {
    for (int l=0;l<result.matrices();++l )
      result[l] *= scale*weight[l];
  }

  void invMassZ(basics::Matrices& result, const basics::Vector& weight, Real scale)
  {
    for (int l=0;l<result.matrices();++l )
      result[l] *= Real(1)/(scale*weight[l]);
  }

  void applyHomogenousDirichletBC(basics::Matrices& u)
  {
    u[0].clear();
    for( int l=1;l<u.matrices()-1;++l )
      applyHomogenousDirichletBC(u[l]);
    u[u.matrices()-1].clear();
  }

  //basics::Matrices reconstruct(const basics::Matrices& u, const int n_Nx, const int n_Ny, const int n_Nz, const basics::Vector& gridX, const basics::Vector& gridY, const basics::Vector& gridZ)
  //{
  //    basics::Vector gridX2("grid twice the size",n_Nx,false);
  //    basics::Vector weightX2("weight twice the size",n_Nx+1,false);
  //    utilities::GLL::GaussLobattoLegendreGridPeriodic(gridX2);
  //    utilities::GLL::GaussLobattoLegendreWeightsPeriodic(weightX2,gridX2);
  //    gridX2 += 1;
  //    gridX2 *= M_PI;
  //    basics::Matrix interpolX = utilities::GLL::periodicInterpolationMatrix(gridX,gridX2);

  //    basics::Vector gridY2("grid twice the size",n_Ny,false);
  //    basics::Vector weightY2("weight twice the size",n_Ny,false);
  //    utilities::GLL::GaussLobattoLegendreGrid(gridY2);
  //    utilities::GLL::GaussLobattoLegendreWeights(weightY2,gridY2);
  //    basics::Matrix interpolY = utilities::GLL::interpolationMatrix(gridY,gridY2);
  //    
  //    basics::Vector gridZ2("grid twice the size",n_Nz,false);
  //    basics::Vector weightZ2("weight twice the size",n_Nz,false);
  //    utilities::GLL::GaussLobattoLegendreGridPeriodic(gridZ2);
  //    utilities::GLL::GaussLobattoLegendreWeightsPeriodic(weightZ2,gridZ2);
  //    gridZ2 += 1;
  //    gridZ2 *= M_PI;
  //    basics::Matrix interpolZ = utilities::GLL::periodicInterpolationMatrix(gridZ,gridZ2);

  //    basics::Matrices V("reconstructed in X",n_Nx,u.cols(),u.matrices());
  //    interpolateX(V,interpolX,u);
  //    basics::Matrices V2("reconstructed in X and Y",n_Nx,n_Ny,u.matrices());
  //    interpolateY(V2,interpolY,V);
  //    basics::Matrices V3("reconstructed fully",n_Nx,n_Ny,n_Nz);
  //    interpolateZ(V3,interpolZ,V2);
  //    
  //    return V3;
  //}

  //Real L2norm(const basics::Matrices& u, const basics::Vector& weightX, const basics::Vector& weightY, const basics::Vector& weightZ)
  //{
  //    basics::Matrices temp(u);
  //  massX(temp,weightX);
  //  massY(temp,weightY);
  //    massZ(temp,weightZ);

  //    Real L2U2 = u.dot(temp);
  //    return( sqrt(L2U2) );
  //}

  //Real convergence(const basics::Matrices& u, const basics::Matrices& exact, const basics::Vector& gridX, const basics::Vector& gridY, const basics::Vector& gridZ)
  //{
  //    basics::Matrices u2 = reconstruct(u,exact.rows(),exact.cols(),exact.matrices(),gridX,gridY,gridZ);
  //    u2 -= exact;
  //    
  //    basics::Vector gridX2("grid twice the size",exact.rows(),false);
  //    basics::Vector weightX2("weight twice the size",exact.rows()+1,false);
  //    utilities::GLL::GaussLobattoLegendreGridPeriodic(gridX2);
  //    utilities::GLL::GaussLobattoLegendreWeightsPeriodic(weightX2,gridX2);
  //    
  //    basics::Vector gridY2("grid twice the size",exact.cols(),false);
  //    basics::Vector weightY2("weight twice the size",exact.cols(),false);
  //    utilities::GLL::GaussLobattoLegendreGrid(gridY2);
  //    utilities::GLL::GaussLobattoLegendreWeights(weightY2,gridY2);
  //    
  //    basics::Vector gridZ2("grid twice the size",exact.matrices(),false);
  //    basics::Vector weightZ2("weight twice the size",exact.matrices(),false);
  //    utilities::GLL::GaussLobattoLegendreGridPeriodic(gridZ2);
  //    utilities::GLL::GaussLobattoLegendreWeightsPeriodic(weightZ2,gridZ2);
  //    
  //    return( L2norm(u2,weightX2,weightY2,weightZ2)/L2norm(exact,weightX2,weightY2,weightZ2) );
  //}

  //void axpyPreserveBorders(basics::Matrices& result, const basics::Matrices& phi, const Real factor)
  //{
  //    int N=result.rows()*(result.cols()-2);
  //    for( int l=0;l<result.matrices();++l )  // need to preserve borders
  //        BLASRPFX(axpy,BLASINT(N),BLASREAL(factor),BLASPREAL(phi[l].data()[1]),BLASINT(mnlIntOne),BLASPREAL(result[l].data()[1]),BLASINT(mnlIntOne));
  //}

  //void scalePreserveBorders(basics::Matrices& result, const Real factor)
  //{
  //    int N=result.rows()*(result.cols()-2);
  //    for( int l=0;l<result.matrices();++l )  // need to preserve borders
  //        BLASRPFX(scal,BLASINT(N),BLASREAL(factor),BLASPREAL(result[l].data()[1]),BLASINT(mnlIntOne));
  //}

  //void copyPreserveBorders(basics::Matrices& result, const basics::Matrices& phi)
  //{
  //    int N=result.rows()*(result.cols()-2);
  //    for( int l=0;l<result.matrices();++l )
  //        BLASRPFX(copy,BLASINT(N),BLASPREAL(phi[l].data()[1]),BLASINT(mnlIntOne),BLASPREAL(result[l].data()[1]),BLASINT(mnlIntOne));
  //}

  //void clearBorders(basics::Matrices& result)
  //{
  //    for( int l=0;l<result.matrices();++l )
  //        for (int k=0;k<result.rows();++k )
  //            result[l][0][k] = result[l][result.cols()-1][k] = 0.f;
  //}

  //void diffx(basics::complexMatrices& result, const basics::Vector& waveNumber, Real scale)
  //{
  //    LOG("find x derivative");
  //#pragma omp parallel for schedule(static)
  //    for( int l=0;l<result.matrices();++l )
  //        for( int k=0;k<result.rows();++k ) {
  //            std::complex<Real> kx(0,scale*waveNumber[k]);
  //            BLASCPFX(scal,BLASINT(result.cols()),BLASCOMPLEX(kx),BLASPCOMPLEX(result[l].data()[0]+k),BLASINT(result.rows()));
  //        }
  //}

  //void curlcurl(basics::complexField2D& result, const basics::complexField2D& u, const basics::complexMatrix& D, const basics::complexMatrix& D2, const basics::Vector& Kx)
  //{
  //    basics::complexMatrix* buffer = utilities::g_manager.aquireComplexMatrix("buffer",u.cols(),u.rows());
  //    basics::complexMatrix* buffer2 = utilities::g_manager.aquireComplexMatrix("buffer",u.cols(),u.rows());

  //    int c=D.cols()-2;
  //    BLASCPFX(gemm,BLASCHAR(mnlNoTrans),BLASCHAR(mnlTrans),BLASINT(D.rows()),BLASINT(u.rows()),BLASINT(c),BLASCOMPLEX(mnlComplexOne),BLASPCOMPLEX(D.data()[1]),BLASINT(D.rows()),BLASPCOMPLEX(u.Y().data()[1]),BLASINT(u.rows()),BLASCOMPLEX(mnlComplexZero),BLASPCOMPLEX(buffer->data()[0]),BLASINT(buffer->rows())); // duy/dy
  //    BLASCPFX(gemm,BLASCHAR(mnlNoTrans),BLASCHAR(mnlTrans),BLASINT(D2.rows()),BLASINT(u.rows()),BLASINT(c),BLASCOMPLEX(mnlComplexOne),BLASPCOMPLEX(D2.data()[1]),BLASINT(D2.rows()),BLASPCOMPLEX(u.X().data()[1]),BLASINT(u.rows()),BLASCOMPLEX(mnlComplexZero),BLASPCOMPLEX(buffer2->data()[0]),BLASINT(buffer2->rows()));
  //    for( int k=0;k<u.rows();++k ) {
  //        result.X().copyRow(k,(*buffer2)[k]);
  //        result.X().axpyRow(k,std::complex<Real>(0,Kx[k]),(*buffer)[k]); // d2ux
  //    }

  //    BLASCPFX(gemm,BLASCHAR(mnlNoTrans),BLASCHAR(mnlTrans),BLASINT(D.rows()),BLASINT(u.rows()),BLASINT(c),BLASCOMPLEX(mnlComplexOne),BLASPCOMPLEX(D.data()[1]),BLASINT(D.rows()),BLASPCOMPLEX(u.X().data()[1]),BLASINT(u.rows()),BLASCOMPLEX(mnlComplexZero),BLASPCOMPLEX(buffer->data()[0]),BLASINT(buffer->rows())); // dux/dy
  //    for( int k=0;k<u.rows();++k ) {
  //        result.Y().copyRow(k,(*buffer)[k]); // d2ux/dxdy
  //        result.Y().scaleRow(k,std::complex<Real>(0,Kx[k]));
  //        result.Y().axpyRow(k,-Kx[k]*Kx[k],u.X()); // d2uy/dx2
  //    }

  //    utilities::g_manager.unlock(buffer);
  //}

  //void diffy(basics::complexMatrices& result, const basics::complexMatrix& D)
  //{
  //    assert( result.cols() == D.cols() );
  //    LOG("find y derivative");

  //#pragma omp parallel for schedule(static)
  //    for( int l=0;l<result.matrices();++l ) { // cannot do this inplace due to the transpose
  //        basics::complexMatrix* T = utilities::g_manager.aquireComplexMatrix("buffer", result.rows(),result.cols());
  //        BLASCPFX(gemm,BLASCHAR(mnlNoTrans),BLASCHAR(mnlTrans),BLASINT(result.rows()),BLASINT(D.rows()),BLASINT(result.cols()),BLASCOMPLEX(mnlComplexOne),BLASPCOMPLEX(result[l].data()[0]),BLASINT(result.rows()),BLASPCOMPLEX(D.data()[0]),BLASINT(D.cols()),BLASCOMPLEX(mnlComplexZero),BLASPCOMPLEX(T->data()[0]),BLASINT(result.rows())); 
  //        result[l] = *T;
  //        utilities::g_manager.unlock(T);
  //    }
  //}

  //void diffz(basics::complexMatrices& result, const basics::Vector& waveNumber, Real scale)
  //{
  //    LOG("find z derivative");
  //#pragma omp parallel for schedule(static)
  //    for( int l=0;l<result.matrices();++l )  {
  //        std::complex<Real> kz(0,waveNumber[l]*scale);
  //        result[l] *= kz;
  //    }
  //}

  //void interpolateX(basics::Matrices& result, const basics::Matrix& I, const basics::Matrices& U, Real alpha, Real beta)
  //{
  //#pragma omp parallel for schedule(static)
  //    for( int l=0;l<result.matrices();++l )
  //        basics::multTranspose(result[l],I,U[l],'N','N',alpha,beta);
  //}

  //void interpolateY(basics::Matrices& result, const basics::Matrix& I, const basics::Matrices& U, Real alpha, Real beta)
  //{
  //#pragma omp parallel for schedule(static)
  //    for( int l=0;l<result.matrices();++l )
  //        basics::multTranspose(result[l],U[l],I,'N','T',alpha,beta);
  //}

  //void interpolateZ(basics::Matrices& result, const basics::Matrix& I, const basics::Matrices& U, Real alpha, Real beta)
  //{
  //    int N=U.rows()*U.cols();
  //#pragma omp parallel for schedule(static)
  //    for( int j=0;j<result.cols();++j )
  //        for( int k=0;k<result.rows();++k )
  //            BLASRPFX(gemv,BLASCHAR(mnlNoTrans),BLASINT(I.rows()),BLASINT(I.cols()),BLASREAL(mnlRealOne),BLASPREAL(I.data()[0]),BLASINT(I.rows()),BLASPREAL(U[0].data()[j]+k),BLASINT(N),BLASREAL(mnlRealZero),BLASPREAL(result[0].data()[j]+k),BLASINT(N));
  //}
}

