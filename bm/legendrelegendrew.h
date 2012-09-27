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

#ifndef LEGENDRE_LEGENDREW_H_
#define LEGENDRE_LEGENDREW_H_

#include "mnl/gll.h"
#include "mnl/matrix.h"
#include "mnl/matrices.h"
#include "mnl/field.h"
#include "mnl/vector.h"
#include "mnl/cgsolver.h"
#include "mnl/function.h"

#include <string>

namespace legendreLegendreW {
  /* 2D */
  mnl::basics::Matrix reconstruct(const mnl::basics::Matrix& u, const int n_Nx, const int n_Ny, const mnl::basics::Vector& gridX, const mnl::basics::Vector& gridY);
  void removeHydrostaticMode(mnl::basics::Matrix& u, const mnl::basics::Vector& weightX, const mnl::basics::Vector& weightY);
  Real L2norm(const mnl::basics::Matrix& u, const mnl::basics::Vector& weightX, const mnl::basics::Vector& weightY);
  Real convergence(const mnl::basics::Matrix& u, const mnl::basics::Vector& gridX, const mnl::basics::Vector& gridY, const mnl::basics::function2D& exactFunc, const Real t, bool relative=true);
  Real convergenceH1(const mnl::basics::Matrix& u, const mnl::basics::Vector& gridX, const mnl::basics::Vector& gridY, const mnl::basics::function2D& exactFunc, const Real t);
  void copyPreserveBorders(mnl::basics::Matrix& result, const mnl::basics::Matrix& phi);
  inline void copyPreserveBorders(mnl::basics::Field2D& result, const mnl::basics::Field2D& phi)
  {
    copyPreserveBorders(result.X(),phi.X()); 
    copyPreserveBorders(result.Y(),phi.Y()); 
  }
  void axpyPreserveBorders(mnl::basics::Matrix& result, const mnl::basics::Matrix& phi, const Real factor=Real(1));
  inline void axpyPreserveBorders(mnl::basics::Field2D& result, const mnl::basics::Field2D& phi, const Real factor=Real(1))
  {
    axpyPreserveBorders(result.X(),phi.X(),factor);
    axpyPreserveBorders(result.Y(),phi.Y(),factor);
  }
  void scalePreserveBorders(mnl::basics::Matrix& result, const Real factor);
  inline void scalePreserveBorders(mnl::basics::Field2D& result, const Real factor)
  {
    scalePreserveBorders(result.X(),factor);
    scalePreserveBorders(result.Y(),factor);
  }
  void interpolate(mnl::basics::Matrix& result, const mnl::basics::Matrix& Ix, const mnl::basics::Matrix& Iy, const mnl::basics::Matrix& U, Real alpha=mnlRealOne, Real beta=mnlRealZero);
  void interpolateX(mnl::basics::Matrix& result, const mnl::basics::Matrix& I, const mnl::basics::Matrix& U, Real alpha=mnlRealOne, Real beta=mnlRealZero);
  void interpolateY(mnl::basics::Matrix& result, const mnl::basics::Matrix& I, const mnl::basics::Matrix& U, Real alpha=mnlRealOne, Real beta=mnlRealZero);
  void diffx(mnl::basics::Matrix& result, const mnl::basics::Matrix& D, Real scale=mnlRealOne);
  void diffx(mnl::basics::Matrix& result, const mnl::basics::Matrix& D, const mnl::basics::Matrix& input, const Real scaleA=mnlRealOne, const Real scaleB=mnlRealZero);
  void diffy(mnl::basics::Matrix& result, const mnl::basics::Matrix& D, const Real scaleInput=mnlRealOne);
  void diffy(mnl::basics::Matrix& result, const mnl::basics::Matrix& D, const mnl::basics::Matrix& input, const Real scaleOutput=mnlRealOne, const Real scaleInput=mnlRealOne);
  void massX(mnl::basics::Matrix& result, const mnl::basics::Vector& weight, Real scale=M_PI);
  inline void massX(mnl::basics::Field2D& result, const mnl::basics::Vector& weight, Real scale=M_PI)
  {
    massX(result.X(),weight,scale);
    massX(result.Y(),weight,scale);
  }
  void invMassX(mnl::basics::Matrix& result, const mnl::basics::Vector& weight, Real scale=M_PI);
  inline void invMassX(mnl::basics::Field2D& result, const mnl::basics::Vector& weight, Real scale=M_PI)
  {
    invMassX(result.X(),weight,scale);
    invMassX(result.Y(),weight,scale);
  }
  void massY(mnl::basics::Matrix& result, const mnl::basics::Vector& weight, Real scale=1);
  inline void massY(mnl::basics::Field2D& result, const mnl::basics::Vector& weight, Real scale=1)
  {
    massY(result.X(),weight,scale);
    massY(result.Y(),weight,scale);
  }
  inline void mass(mnl::basics::Matrix& result, const mnl::basics::Vector& weightX,
      const mnl::basics::Vector& weightY, Real scaleX=M_PI)
  {
    massX(result,weightX,scaleX);
    massY(result,weightY);
  }
  inline void mass(mnl::basics::Field2D& result, const mnl::basics::Vector& weightX,
      const mnl::basics::Vector& weightY, Real scaleX=M_PI)
  {
    mass(result.X(),weightX,weightY,scaleX);
    mass(result.Y(),weightX,weightY,scaleX);
  }
  void invMassY(mnl::basics::Matrix& result, const mnl::basics::Vector& weight, Real scale=1);
  inline void invMassY(mnl::basics::Field2D& result, const mnl::basics::Vector& weight, Real scale=1)
  {
    invMassY(result.X(),weight,scale);
    invMassY(result.Y(),weight,scale);
  }
  std::string errorReport(const mnl::basics::Matrix& p, const mnl::basics::Field2D& u, const mnl::basics::Vector& gridX, const mnl::basics::Vector& gridY, const mnl::basics::function2D& pExact, const mnl::basics::function2D& uxExact, const mnl::basics::function2D& uyExact, const Real t);
  std::string errorReport(const mnl::basics::Field2D& u, const mnl::basics::Vector& gridX, const mnl::basics::Vector& gridY, const mnl::basics::function2D& uxExact, const mnl::basics::function2D& uyExact, const Real t);
  std::string errorReport(const mnl::basics::Matrix& u, const mnl::basics::Vector& gridX, const mnl::basics::Vector& gridY, const mnl::basics::function2D& uExact, const Real t);
  void applyHomogenousDirichletBC(mnl::basics::Matrix& u);
  inline void applyHomogenousDirichletBC(mnl::basics::Field2D& u)
  {
    applyHomogenousDirichletBC(u.X());
    applyHomogenousDirichletBC(u.Y());
  }
  //     void curlcurl(mnl::basics::Field2D& result, const mnl::basics::Field2D& u, const mnl::basics::Matrix& Dx, const mnl::basics::Matrix& Dy, const mnl::basics::complexMatrix& D2, const mnl::basics::Vector& Kx);

  /* 3D */
  //    mnl::basics::Matrices reconstruct(const mnl::basics::Matrix& u, const int n_Nx, const int n_Ny, const int n_Nz, const mnl::basics::Vector& gridX, const mnl::basics::Vector& gridY, const mnl::basics::Vector& gridZ);
  //    Real L2norm(const mnl::basics::Matrices& u, const mnl::basics::Vector& weightX, const mnl::basics::Vector& weightY, const mnl::basics::Vector& weightZ);
  //    Real convergence(const mnl::basics::Matrices& u, const mnl::basics::Matrices& exact, const mnl::basics::Vector& gridX, const mnl::basics::Vector& gridY, const mnl::basics::Vector& gridZ);
  //    Real convergenceH1(const mnl::basics::Matrix& u, const mnl::basics::Matrix& exact, const mnl::basics::Vector& weightX, const mnl::basics::Vector& weightY);
  //    void copyPreserveBorders(mnl::basics::Matrices& result, const mnl::basics::Matrices& phi);
  //    void axpyPreserveBorders(mnl::basics::Matrices& result, const mnl::basics::Matrices& phi, const Real factor);
  //    void scalePreserveBorders(mnl::basics::Matrices& result, const Real factor);
  //    void interpolateX(mnl::basics::Matrices& result, const mnl::basics::Matrix& I, const mnl::basics::Matrices& U, Real alpha=mnlRealOne, Real beta=mnlRealZero);
  //    void interpolateY(mnl::basics::Matrices& result, const mnl::basics::Matrix& I, const mnl::basics::Matrices& U, Real alpha=mnlRealOne, Real beta=mnlRealZero);
  //    void interpolateZ(mnl::basics::Matrices& result, const mnl::basics::Matrix& I, const mnl::basics::Matrices& U, Real alpha=mnlRealOne, Real beta=mnlRealZero);
  //    void diffx(mnl::basics::Matrices& result, const mnl::basics::Matrix& D, Real scale=1.f);
  //    void diffx(mnl::basics::Matrices& result, const mnl::basics::Matrix& D, const mnl::basics::Matrices& input, Real scale=1.f);
  //    void diffy(mnl::basics::Matrices& result, const mnl::basics::Matrix& D);
  //    void diffy(mnl::basics::Matrices& result, const mnl::basics::Matrix& D, const mnl::basics::Matrices& input, const Real scaleOutput=mnlRealOne, const Real scaleInput=mnlRealOne);
  //    void diffz(mnl::basics::Matrices& result, const mnl::basics::Matrix& D);
  //    void diffz(mnl::basics::Matrices& result, const mnl::basics::Matrix& D, const mnl::basics::Matrices& input, const Real scaleOutput=mnlRealOne, const Real scaleInput=mnlRealOne);

  void massX(mnl::basics::Matrices& result, const mnl::basics::Vector& weight, Real scale=M_PI);
  inline void massX(mnl::basics::Field3D& result, const mnl::basics::Vector& weight, Real scale=M_PI)
  {
    massX(result.X(),weight,scale);
    massX(result.Y(),weight,scale);
    massX(result.Z(),weight,scale);
  }
  void invMassX(mnl::basics::Matrices& result, const mnl::basics::Vector& weight, Real scale=M_PI);
  inline void invMassX(mnl::basics::Field3D& result, const mnl::basics::Vector& weight, Real scale=M_PI)
  {
    invMassX(result.X(),weight,scale);
    invMassX(result.Y(),weight,scale);
    invMassX(result.Y(),weight,scale);
  }
  void massY(mnl::basics::Matrices& result, const mnl::basics::Vector& weight, Real scale=1);
  inline void massY(mnl::basics::Field3D& result, const mnl::basics::Vector& weight, Real scale=1)
  {
    massY(result.X(),weight,scale);
    massY(result.Y(),weight,scale);
    massY(result.Z(),weight,scale);
  }
  void invMassY(mnl::basics::Matrices& result, const mnl::basics::Vector& weight, Real scale=1);
  inline void invMassY(mnl::basics::Field3D& result, const mnl::basics::Vector& weight, Real scale=1)
  {
    invMassY(result.X(),weight,scale);
    invMassY(result.Y(),weight,scale);
    invMassY(result.Z(),weight,scale);
  }
  void massZ(mnl::basics::Matrices& result, const mnl::basics::Vector& weight, Real scale=1);
  inline void massZ(mnl::basics::Field3D& result, const mnl::basics::Vector& weight, Real scale=1)
  {
    massZ(result.X(),weight,scale);
    massZ(result.Y(),weight,scale);
    massZ(result.Z(),weight,scale);
  }
  void invMassZ(mnl::basics::Matrices& result, const mnl::basics::Vector& weight, Real scale=1);
  inline void invMassZ(mnl::basics::Field3D& result, const mnl::basics::Vector& weight, Real scale=1)
  {
    invMassZ(result.X(),weight,scale);
    invMassZ(result.Y(),weight,scale);
    invMassZ(result.Z(),weight,scale);
  }

  void applyHomogenousDirichletBC(mnl::basics::Matrices& u);
  inline void applyHomogenousDirichletBC(mnl::basics::Field3D& u)
  {
    applyHomogenousDirichletBC(u.X());
    applyHomogenousDirichletBC(u.Y());
    applyHomogenousDirichletBC(u.Z());
  }

  class laplacianEvaluator2D : public mnl::utilities::Evaluator<mnl::basics::Matrix> {
    public:
      laplacianEvaluator2D(const mnl::basics::Matrix& Ax, const mnl::basics::Matrix& Ay, const mnl::basics::Vector& weightX, const mnl::basics::Vector& weightY)
        : m_Ax(Ax), m_Ay(Ay), m_weightX(weightX), m_weightY(weightY), m_nu(1), m_mu(0)
      {
      }

      virtual ~laplacianEvaluator2D() {};

      void setNu(Real nu)
      {
        m_nu = nu;
      }

      void setMu(Real mu)
      {
        m_mu = mu;
      }

      virtual void evaluate(mnl::basics::Matrix& mat, const mnl::basics::Matrix& mat2) const;
      inline void evaluate(mnl::basics::Field2D& output, const mnl::basics::Field2D& input) const
      {
        evaluate(output.X(),input.X());
        evaluate(output.Y(),input.Y());
      }
    protected:
      const mnl::basics::Matrix& m_Ax;
      const mnl::basics::Matrix& m_Ay;
      const mnl::basics::Vector& m_weightX;
      const mnl::basics::Vector& m_weightY;
      Real m_nu;
      Real m_mu;
  };
}

#endif

