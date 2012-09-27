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
#ifndef MNL_GORDON_HALL_H_
#define MNL_GORDON_HALL_H_

#include "config.h"
#include "vector.h"
#include "matrix.h"
#include "matrices.h"
#include "function.h"
#include "buffers.h"

namespace mnl {
  namespace utilities {
    class GordonHall {
      public:
        GordonHall(int N, int M);

        GordonHall(int N, int M, basics::function1Dto2D& bottom, basics::function1Dto2D& right, 
            basics::function1Dto2D& top, basics::function1Dto2D& left);

        void compute(Real t);

        inline basics::Matrix& getJacobian()
        {
          return( m_J );
        }

        const basics::Matrix& getJacobian() const
        {
          return( m_J );
        }

        inline basics::matrixStack& getGeometryDerivatives()
        {
          return( m_D );
        }

        inline const basics::matrixStack& getGeometryDerivatives() const
        {
          return( m_D );
        }

        inline basics::Field2D& getMapping()
        {
          return( m_map );
        }

        inline const basics::Field2D& getMapping() const
        {
          return( m_map );
        }

        static void compute(basics::Matrix& x, basics::Matrix& y, const basics::Vector& grid); 

        /* evaluate the resulting polynomial  */
        static std::pair<Real,Real> evaluate(Real x, Real y, const basics::Matrix& X, const basics::Matrix& Y, const basics::Vector& grid);
        inline std::pair<Real,Real> evaluate(Real x, Real y, const basics::Field2D& map, const basics::Vector& grid) const
        {
          return( evaluate(x,y,map.X(),map.Y(),grid) );
        }
        inline std::pair<Real,Real> evaluate(Real x, Real y, const basics::Vector& grid) const
        {
          return( evaluate(x,y,getMapping(),grid) );
        }
        static void compute(basics::Field2D& map, const basics::function1Dto2D& g1, const basics::function1Dto2D& g2, const basics::function1Dto2D& g3, const basics::function1Dto2D& g4, Real t=0.f); // t=0 in function evaluations by default
        void computeJacobian(basics::Matrix& J, const basics::Field2D& map);
        inline void computeJacobian()
        {
          computeJacobian(m_J,m_map);
        }
      protected:
        static void compute(basics::Matrix& x, const basics::Vector& grid);

        basics::Matrix 			 m_J;
        basics::Matrix			 m_Dx;
        basics::Matrix			 m_Dy;
        basics::Field2D 		 m_map;
        basics::function1Dto2D*	 m_bottom;
        basics::function1Dto2D*	 m_right; 
        basics::function1Dto2D*	 m_top;
        basics::function1Dto2D*	 m_left;
        basics::matrixStack      m_D;
    };

    class GordonHall3D {
      public:
        GordonHall3D(int N, int M, int P);

        inline basics::Matrices& getJacobian()
        {
          return( m_J );
        }

        const basics::Matrices& getJacobian() const
        {
          return( m_J );
        }

        inline basics::matricesStack& getGeometryDerivatives()
        {
          return( m_D );
        }

        inline const basics::matricesStack& getGeometryDerivatives() const
        {
          return( m_D );
        }

        inline basics::Field3D& getMapping()
        {
          return( m_map );
        }

        inline const basics::Field3D& getMapping() const
        {
          return( m_map );
        }

        void computeJacobian(basics::Matrices& J, const basics::Field3D& map);
        inline void computeJacobian()
        {
          computeJacobian(m_J,m_map);
        }
      protected:
        basics::Matrices		 m_J;
        basics::Matrix		 	 m_Dx;
        basics::Field3D 		 m_map;
        basics::matricesStack	 m_D;
    };
  }
}

#endif

