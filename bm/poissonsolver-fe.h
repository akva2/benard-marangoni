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

#ifndef POISSON_SOLVER_FE_H
#define POISSON_SOLVER_FE_H

#include "mnl/vector.h"
#include "mnl/matrices.h"
#include "mnl/matrices.h"
#include "mnl/function.h"

#include "poissonsolver.h"

namespace linearElements {
  class poissonSolver : public poissonSolverT<mnl::basics::Matrix,mnl::basics::Matrices> {
    public:
      /* Give Nz value 0 if only 2D is desired.. */
      poissonSolver(const BC& h, const mnl::basics::Vector& grid,
          const mnl::basics::Vector& scaleX, const mnl::basics::Vector& scaleY,
          bool L2=false, bool doZ=false, const BC& bz=Nonhomogenous);
      virtual ~poissonSolver();

      void solve(mnl::basics::Matrix& F, const Real mu=Real(0), 
          const Real nu=Real(1), const Real scale=Real(1)) const;
      void solve(mnl::basics::Matrix& F, mnl::basics::Matrix& buffer, const Real mu=Real(0), 
          const Real nu=Real(1), const Real scale=Real(1)) const;
      void solve(mnl::basics::Matrices& F, const Real mu=Real(0), const Real nu=Real(1)) const;
      void solve(mnl::basics::Matrices& F, mnl::basics::Matrices& buffer,
          const Real Lx=2, const Real Ly=2, const Real Lz=2, const Real mu=0) const;

      inline const mnl::basics::Matrix& Bx() const
      {
        return( *m_Bx );
      }
      inline const mnl::basics::Matrix& By() const
      {
        return( *m_By );
      }
    protected:
      mnl::basics::Matrix* m_Bx;
      mnl::basics::Matrix* m_By;
      mnl::basics::Matrix* m_Bz;
      mnl::basics::Matrices* m_buffer;
      bool m_L2;
      BC m_bz;
    private:
      void assemble(mnl::basics::Matrix& A, mnl::basics::Matrix& B,
          const mnl::basics::Vector& grid);
      void setupH1(const BC& h, const mnl::basics::Vector& grid,
          const mnl::basics::Vector& scaleX,
          const mnl::basics::Vector& scaleY, bool doZ, const BC& bz);
      void setupL2(const BC& h, const mnl::basics::Vector& grid,
          const mnl::basics::Vector& scaleX,
          const mnl::basics::Vector& scaleY, bool doZ);
      mnl::basics::Vector* setupH1Grid(const BC& h, const mnl::basics::Vector& grid,
          const mnl::basics::Vector& scale);
      mnl::basics::Vector* setupL2Grid(const BC& h, const mnl::basics::Vector& grid,
          const mnl::basics::Vector& scale);
      void assembleAndMask(mnl::basics::Matrix*& A,
          mnl::basics::Matrix*& B,
          const mnl::basics::Vector& grid, BC h);
  };
}
#endif

