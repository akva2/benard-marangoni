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
#ifndef MNL_GLL_H_
#define MNL_GLL_H_

#include "config.h"
#include "vector.h"
#include "matrix.h"

/** Muu numerics library covers the basics of numerical computations, including data structures such
	* as matrices and vectors and useful utilities such as GLL quadrature and finding Fourier transforms. */
namespace mnl {
  /// Namespace for useful functions integrated with mnl data types.
  namespace utilities {
    /**	GLL quadrature class 
      - ported straight from the matlab code by
      Kay Hansen-Zahl and Tormod Bjøntegaard (hope they don't mind..)
     **/
    class GLL {
      public:
        static Real Legendre(const Real x, const int N);
        static Real LegendreDerivative(const Real x, const int N);
        /* evalutes lagrangian cardinal polynomial nr at x using grid as gridpoints */
        static Real Lagrange(const Real x, const int nr, const basics::Vector& grid);
        static void GaussLegendreGrid(basics::Vector& grid);
        static void GaussLegendreWeights(basics::Vector& weight, const basics::Vector& grid);
        static void GaussLobattoLegendreGrid(basics::Vector& grid);
        static void GaussLobattoLegendreGridPeriodic(basics::Vector& grid);
        static void GaussLobattoLegendreGridHomogenous(basics::Vector& grid);
        static void GaussLobattoLegendreHomogenous(basics::Vector& grid);
        static void GaussLobattoLegendreWeights(basics::Vector& weight, const basics::Vector& grid);
        static void GaussLobattoLegendreWeightsPeriodic(basics::Vector& weight, const basics::Vector& grid);
        static void GaussLobattoLegendreWeightsHomogenous(basics::Vector& weight, const basics::Vector& grid);

        static basics::Matrix LagrangeDerivativeMatrix(int N);
        static void LagrangeDerivativeMatrix(basics::Matrix& A);
        static void LagrangeDerivativeMatrix(basics::Matrix& A, const basics::Vector& grid);
        static void LagrangeDerivativeMatrixPeriodic(basics::Matrix& A);
        static void LagrangeDerivativeMatrixPeriodic(basics::Matrix& A, const basics::Vector& grid);
        static void LagrangeDerivativeMatrixHomogenous(basics::Matrix& A, const basics::Vector& grid);

        static basics::Matrix interpolationMatrix(const basics::Vector& gridFrom, int M);
        static basics::Matrix interpolationMatrix(const basics::Vector& gridFrom, const basics::Vector& gridTo);
        static basics::Matrix periodicInterpolationMatrix(const basics::Vector& gridFrom, const basics::Vector& gridTo);
    };
  }
}

#endif

