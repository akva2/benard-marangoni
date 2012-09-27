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
#ifndef MNL_CONFIG_H_
#define MNL_CONFIG_H_

/* Define to work in double precision */
#define MNL_USE_DOUBLE_
/* Define this to get memory allocation/deallocation logging */
//#define MNL_MEMORY_VERBOSE_
/* Define make the routines be verbose */
//#define MNL_VERBOSE_

#include <complex>
#include <iostream>
#include <vector>

extern "C" {
#ifdef __sun__
#include <sunperf.h>
#endif

#ifdef __freebsd__
#include "/usr/local/include/blas.h"
#include "/usr/local/include/lapack.h"
#endif

#ifdef __INTEL_COMPILER
#include <mkl.h>
#elif defined(__linux__)
#include "f2c.h"
#include "blas.h"
#define LAPACKMANGLE(name) name##_
#include "clapack.h"
#endif

#ifdef _AIX
#define LAPACKMANGLE(name) name
namespace LAPACK {
#include "clapack.h"
}
#include <essl.h>
#endif
}

#ifdef MNL_USE_DOUBLE_
typedef double Real;
#ifdef _AIX
typedef std::complex<double> Complex;
#elif defined(__INTEL_COMPILER)
typedef MKL_Complex16 Complex;
#else
typedef doublecomplex Complex;
#endif
#ifdef __sun__
#define BLASRPFX(name,...) d##name(__VA_ARGS__)
#define BLASCPFX(name,...) z##name(__VA_ARGS__)
#define LAPACKRPFX(name,...) d##name##_(__VA_ARGS__)
#define LAPACKCPFX(name,...) z##name##_(__VA_ARGS__)
#define ZDOTC(res,...) doublecomplex zdotcresult; zdotcresult = BLASCPFX(dotc,__VA_ARGS__); memcpy(&res,&zdotcresult,2*sizeof(Real));
  #endif

  #ifdef _AIX
    #define BLASRPFX(name,...) d##name(__VA_ARGS__)
    #define BLASCPFX(name,...) z##name(__VA_ARGS__)
    #define LAPACKRPFX(name,...) LAPACK::d##name##(__VA_ARGS__)  
    #define LAPACKCPFX(name,...) LAPACK::z##name##(__VA_ARGS__)  
    #define ZDOTC(res,...) res = BLASCPFX(dotc,__VA_ARGS__)
  #endif

  #if defined(__INTEL_COMPILER)
    #define BLASRPFX(name,...) d##name(__VA_ARGS__)
    #define BLASCPFX(name,...) z##name(__VA_ARGS__)
    #define LAPACKRPFX(name,...) d##name(__VA_ARGS__)
    #define LAPACKCPFX(name,...) z##name(__VA_ARGS__)
    #define ZDOTC(res,...)   Complex zdotcresult; BLASCPFX(dotc,&zdotcresult,__VA_ARGS__); memcpy(&res,&zdotcresult,2*sizeof(Real));
  #elif defined(__freebsd__) || defined(__linux__)
    #define BLASRPFX(name,...) d##name##_(__VA_ARGS__)
    #define BLASCPFX(name,...) z##name##_(__VA_ARGS__)
    #define LAPACKRPFX(name,...) d##name##_(__VA_ARGS__)
    #define LAPACKCPFX(name,...) z##name##_(__VA_ARGS__)
    #define ZDOTC(res,...)   doublecomplex zdotcresult; BLASCPFX(dotc,&zdotcresult,__VA_ARGS__); memcpy(&res,&zdotcresult,2*sizeof(Real));
    /* this define isnt particularly beautiful but it works 
     * never use zdotcresult as a var name! :) 
    */
  #endif
#else
  typedef float Real;
  typedef complex Complex;
  #if defined(__sun__) || defined(__INTEL_COMPILER)
    #define BLASRPFX(name,...) s##name(__VA_ARGS__)
    #define BLASCPFX(name,...) c##name(__VA_ARGS__)
    #define LAPACKRPFX(name,...) s##name##_(__VA_ARGS__)
    #define LAPACKCPFX(name,...) c##name##_(__VA_ARGS__)
    #define ZDOTC(res,...) res = BLASCPFX(dotc,__VA_ARGS__)
  #else
    #define BLASRPFX(name,...) s##name##_(__VA_ARGS__)
    #define BLASCPFX(name,...) c##name##_(__VA_ARGS__)
    #define LAPACKRPFX(name,...) s##name##_(__VA_ARGS__)
    #define LAPACKCPFX(name,...) c##name##_(__VA_ARGS__)
    #define ZDOTC(res,...) BLASCPFX(dotc,&res,__VA_ARGS__)
  #endif // ifdef __sun__
#endif // ifdef MNL_USE_DOUBLE_

#ifdef __sun__
#define BLASINT(name) name
#define BLASREAL(name) name
#define BLASCHAR(name) name
#define BLASPCOMPLEX(name) const_cast<Complex*>(reinterpret_cast<const Complex*>(name))
#define BLASCOMPLEX(name) BLASPCOMPLEX(&name)
#endif
#ifdef _AIX

#define BLASINT(name) name
#define BLASREAL(name) name
#define BLASCHAR(name) &name
#define BLASCOMPLEX(name) const_cast<Complex&>(reinterpret_cast<const Complex&>(name))
#define BLASPCOMPLEX(name) const_cast<std::complex<Real>* >(name)
#endif

#if defined(__freebsd__) || defined(__linux__)
#define BLASINT(name) const_cast<int*>(&name)
#define BLASREAL(name) const_cast<Real*>(&name)
#define BLASCHAR(name) const_cast<char*>(&name)
#define BLASCOMPLEX(name) BLASPCOMPLEX(&name) 
#define BLASPCOMPLEX(name) const_cast<Complex*>(reinterpret_cast<const Complex*>(name))
#endif

#define BLASPREAL(name) const_cast<Real*>(reinterpret_cast<const Real*>(name))

static const int mnlIntOne=1;
static const int mnlIntTwo=2;
static const int mnlIntMinusOne=1;
static const int mnlIntMinusTwo=2;
static const Real mnlRealOne=1.f;
static const Real mnlRealMinusOne=-1.f;
static const Real mnlRealZero=0.f;
static const std::complex<Real> mnlComplexOne=std::complex<Real>(1.f,0.f);
static const std::complex<Real> mnlComplexMinusOne=std::complex<Real>(-1.f,0.f);
static const std::complex<Real> mnlComplexZero=std::complex<Real>(0.f,0.f);
static const char mnlNoTrans='N';
static const char mnlTrans='T';

#ifdef MNL_MEMORY_VERBOSE_
inline void ALLOCATELOG(const std::string name) { std::cerr << "Allocating " << name << std::endl; }
inline void ALLOCATELOG(const std::string name, int size) { std::cerr << "Allocating " << name << "(" << size << ")" << std::endl; }
inline void DEALLOCATELOG(const std::string name) { std::cerr << "Killing " << name << std::endl; }
#else
inline void ALLOCATELOG(std::string) {}
inline void ALLOCATELOG(std::string,int size) {}
inline void DEALLOCATELOG(std::string) {}
#endif

#ifdef MNL_VERBOSE_
inline void LOG(const std::string name) { std::cout << name << std::endl; }
#else
inline void LOG(const std::string name) {}
#endif

/* needed since it makes <fstream> cry - stupid #defines anyway */
#undef max
#undef min
#undef abs
#undef dabs
#undef dmax
#undef dmin

#endif

