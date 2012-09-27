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
#ifndef MNL_CGSOLVER_H_
#define MNL_CGSOLVER_H_

#include "vector.h"
#include "matrix.h"
#include "matrices.h"
#include "memtracker.h"
#include "buffers.h"
#include "geometry.h"
#include "function.h"
#include "timer.h"
#include "timer.h"

namespace mnl {
  namespace utilities {

    class CGSolver {
      public:
        CGSolver() {};
        ~CGSolver() {};

        template<class T, class DOTTER, class S>
          static int solve(T&, const Evaluator<T>&, const DOTTER& dot, S tol=1.e-10); // CG
        template<class T, class DOTTER, class S>
          static int solve(T&, const Evaluator<T>&, Evaluator<T>&, 
              const DOTTER& dot, S tol=1.e-10); // PCG
    };

    /* template implementation */
    template<class T, class DOTTER, class S>
      int CGSolver::solve(T& x, const Evaluator<T>& A, const DOTTER& dot, S tol) // x should be initialized to b
      {
        S rl=sqrt(dot(x,x));
        if( std::abs<Real>(rl) == 0 )
          return 0;

        T* r      = g_manager.clone(x);
        T* evaled = g_manager.clone(x);
        T* p      = g_manager.clone(x);

        *r = x;
        x.clear();
        int i=0;
        S dotp=1000;
        S rdr = dotp;
        while( i < x.length() && (std::abs<Real>(rdr) > std::abs<Real>(tol) /*|| std::abs<Real>(rdr) > std::abs<Real>(tol*rl)*/) ) {
          ++i;
          if( i == 1 ) {
            *p = *r;
            dotp = std::abs(dot(*r,*r));
          } else {
            S dotp2 = std::abs(dot(*r,*r));
            S beta = dotp2/dotp;
            dotp = dotp2;
            *p *= beta;
            *p += *r;
          }
          A.evaluate(*evaled,*p);
          S alpha = dotp/std::abs(dot(*p,*evaled));
          x.axpy(alpha,*p);
          r->axpy(-alpha,*evaled);
          rdr = sqrt(dot(*r,*r));
          //                    std::cout << rdr << std::endl;
          if( A.isL2() )
            A.filter(*r);
        }
        if( i == 0 ) {
          std::cout << "rdr " << rdr << std::endl;
          std::cout << "rl " << rl << std::endl;
          //                    Real max=0;
          //                    for( int i=0;i<r->size();++i ) {
          //                        Real lmax = *std::max_element((*r)[i].data()[0],(*r)[i].data()[0]+(*r)[i].length());
          //                        if( lmax > max )
          //                            max = lmax;
          //                    }
          //                    std::cout << "max pointwise residual " << max << std::endl;
        }

        g_manager.unlock(r);
        g_manager.unlock(p);
        g_manager.unlock(evaled);

        return( i );
      }

    template<class T, class DOTTER, class S>
      int CGSolver::solve(T& x, const Evaluator<T>& A, Evaluator<T>& P, 
          const DOTTER& dot, S tol) // x should be initialized to b
      {
        S rl=sqrt(dot(x,x));
        if( std::abs<Real>(rl) == 0 )
          return 0;

        T* r      = g_manager.clone(x);
        T* r2=NULL;
        if( P.buffer() )
          r2	  = g_manager.clone(x);
        T* evaled = g_manager.clone(x);
        T* z      = g_manager.clone(x);
        T* p      = g_manager.clone(x);

        *r = x;
        if( r2 )
          *r2 = *P.buffer();
        x.clear();
        int i=0;
        S dotp=1000;
        S rdr = dotp;
        while( i < x.length() && (std::abs<Real>(rdr) > std::abs<Real>(tol)/* || std::abs<Real>(rdr) > std::abs<Real>(tol*rl)*/) ) {
          P.evaluate(*z,*r);
          ++i;
          if( i == 1 ) {
            *p = *z;
            dotp = dot(*r,*z);
          } else {
            S dotp2 = dot(*r,*z);
            S beta = dotp2/dotp;
            dotp = dotp2;
            *p *= beta;
            *p += *z;
          }
          A.evaluate(*evaled,*p);
          S gamma = dot(*p,*evaled);
          S alpha = dotp/gamma;

          x.axpy(alpha,*p);
          r->axpy(-alpha,*evaled);
          if( r2 ) {
            r2->axpy(-alpha,*P.buffer());
            *P.buffer() = *r2;
          }
          //                    utilities::g_profiler.pause(3);
          rdr = sqrt(dot(*r,*r));
          //                    S rdr2 = rdr;
          //                    if( rdr2 < rdr && rdr < 100*tol ) {
          //                        break;
          //                    }
          //                    std::cout << rdr << std::endl;

          if( A.isL2() ) {
            //                        evaled->print();
            //                        std::cout << rdr << std::endl;
            //                        if( i > 300 ) {
            //                            r->save("foo");
            //                            exit(1);
            //                        }
            A.filter(*r);
          }
        }
        if( i == 0 ) {
          x = *r;
          i = 1;
          //                    std::cout << "rdr " << rdr << std::endl;
          //                    std::cout << dot(*r,*r) << std::endl;
          //                    std::cout << "rl " << rl << std::endl;
          //                    Real max=0;
          //                    for( int i=0;i<r->size();++i ) {
          //                        Real lmax = *std::max_element((*r)[i].data(),(*r)[i].data()+(*r)[i].length());
          //                        if( lmax > max )
          //                            max = lmax;
          //                    }
          //                    std::cout << "max pointwise residual " << max << std::endl;
        }


        g_manager.unlock(r);
        g_manager.unlock(evaled);
        g_manager.unlock(z);
        g_manager.unlock(p);
        if( r2 )
          g_manager.unlock(r2);

        return( i );
      }
  };
};

#endif

