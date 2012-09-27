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
#ifndef MNL_GRID_H_
#define MNL_GRID_H_

#include "vector.h"
#include "matrix.h"
#include "matrices.h"
#include "field.h"
#include "function.h"

typedef Real(MNL_FUNC_1D)(const Real x, const Real t);
typedef Real(MNL_FUNC_2D)(const Real x, const Real y, const Real t);
typedef Real(MNL_FUNC_3D)(const Real x, const Real y, const Real z, const Real t);

namespace mnl {
  namespace utilities {
    basics::Vector equidistant(Real start, Real end, int N);
    void equidistant(basics::Vector& result, Real start, Real end);

    basics::Vector wavenumbers(int N);
    void wavenumbers(basics::Vector& result);

    void evaluate(basics::Vector& result, MNL_FUNC_1D& f, const basics::Vector& gridX, const Real t);
    void evaluate(basics::Matrix& result, MNL_FUNC_2D& f, const basics::Vector& gridX, const basics::Vector& gridY, const Real t);
    void evaluate(basics::Matrices& result, MNL_FUNC_3D& f, const basics::Vector& gridX, const basics::Vector& gridY, const basics::Vector& gridZ, const Real t);
    void evaluate(basics::Matrix& result, const basics::function2D& f, const basics::Vector& gridX, const basics::Vector& gridY, const Real t);
    inline void evaluate(basics::Field2D& result, const basics::function2D& fX, const basics::function2D& fY, const basics::Vector& gridX, const basics::Vector& gridY, const Real t)
    {
      evaluate(result.X(),fX,gridX,gridY,t);
      evaluate(result.Y(),fY,gridX,gridY,t);
    }
    void evaluate(basics::Matrices& result, const basics::function3D& f, const basics::Vector& gridX, const basics::Vector& gridY, const basics::Vector& gridZ, const Real t);
    inline void evaluate(basics::Field3D& result, const basics::function3D& fX, const basics::function3D& fY, const basics::function3D& fZ, const basics::Vector& gridX, const basics::Vector& gridY, const basics::Vector& gridZ, const Real t)
    {
      evaluate(result.X(),fX,gridX,gridY,gridZ,t);
      evaluate(result.Y(),fY,gridX,gridY,gridZ,t);
      evaluate(result.Z(),fZ,gridX,gridY,gridZ,t);
    }

    void evaluate(basics::complexVector& result, MNL_FUNC_1D& f, const basics::Vector& gridX, const Real t);
    void evaluate(basics::complexMatrix& result, MNL_FUNC_2D& f, const basics::Vector& gridX, const basics::Vector& gridY, const Real t);
    void evaluate(basics::complexMatrices& result, MNL_FUNC_3D& f, const basics::Vector& gridX, const basics::Vector& gridY, const basics::Vector& gridZ, const Real t);
    void evaluate(basics::complexMatrix& result, const basics::function2D& f, const basics::Vector& gridX, const basics::Vector& gridY, const Real t);
    inline void evaluate(basics::complexField2D& result, const basics::function2D& fX, const basics::function2D& fY, const basics::Vector& gridX, const basics::Vector& gridY, const Real t)
    {
      evaluate(result.X(),fX,gridX,gridY,t);
      evaluate(result.Y(),fY,gridX,gridY,t);
    }
    void evaluate(basics::complexMatrices& result, const basics::function3D& f, const basics::Vector& gridX, const basics::Vector& gridY, const basics::Vector& gridZ, const Real t);
    inline void evaluate(basics::complexField3D& result, const basics::function3D& fX, const basics::function3D& fY, const basics::function3D& fZ, const basics::Vector& gridX, const basics::Vector& gridY, const basics::Vector& gridZ, const Real t)
    {
      evaluate(result.X(),fX,gridX,gridY,gridZ,t);
      evaluate(result.Y(),fY,gridX,gridY,gridZ,t);
      evaluate(result.Z(),fZ,gridX,gridY,gridZ,t);
    }
  }
}

#endif

