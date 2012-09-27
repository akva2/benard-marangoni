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

#ifndef POISSON_SOLVER_LGW_H
#define POISSON_SOLVER_LGW_H

#include "mnl/vector.h"
#include "mnl/matrices.h"
#include "mnl/matrices.h"
#include "mnl/function.h"

#include "poissonsolver.h"

namespace legendreLegendreW {
  class poissonSolver : public poissonSolverT<mnl::basics::Matrix,mnl::basics::Matrices> {
    public:
      /* Give Nz value 0 if only 2D is desired.. */
      poissonSolver(const int Nx, const int Ny, const int Nz,
          const BC h=Homogenous, bool generateOperator=true,
          Real scaleX=2*M_PI, Real scaleY=-1, Real scaleZ=2*M_PI,
          int elemX=1, int elemY=1, int elemZ=1);
      ~poissonSolver();

      void solve(mnl::basics::Matrix& F, const Real mu=Real(0), const Real nu=Real(1), const Real scale=Real(1)) const;
      void solveMult(mnl::basics::Matrix& F, const Real mu=Real(0), const Real nu=Real(1), const Real scale=Real(1)) const;
      void solve(mnl::basics::Matrices& F, const Real mu=Real(0), const Real nu=Real(1)) const;
      //    protected:
      mnl::basics::Matrices* m_buffer;
  };
}
#endif

