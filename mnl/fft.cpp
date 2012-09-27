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

#include "fft.h"

#include <iostream>
#include <assert.h>

namespace mnl {
  namespace fftw {

    /**
     *       Calculates the DFT of an vector.
     * @param vector Input data
     * @return DFT of vector.
     */
    basics::complexVector DFT(const basics::complexVector& vector)
    {
      LOG("Finding DFT of "+vector.name());

      basics::complexVector result("DFT of "+vector.name(),vector.length(),false);
      FFT_PFX(complex*) out =(FFT_PFX(complex*))FFT_PFX(malloc)(sizeof(FFT_PFX(complex))*vector.length());

      FFT_PFX(plan) plan = FFT_PFX(plan_dft_1d)(vector.length(),reinterpret_cast<FFT_PFX(complex*)>(const_cast<std::complex<Real>*>(vector.data())),out,FFTW_FORWARD,FFTW_ESTIMATE);
      assert( plan != NULL );

      FFT_PFX(execute)(plan);

      result = reinterpret_cast<std::complex<Real>*>(out);

      FFT_PFX(destroy_plan)(plan);
      FFT_PFX(free)(out);

      return( result );
    }

    /**
     *       Calculates the DFT of a REAL vector.
     * @param vector The REAL indata. Destroyed during the transform.
     * @return The DFT of the indata packed as described in http://www.fftw.org/fftw3_doc/One_002dDimensional-DFTs-of-Real-Data.html#One_002dDimensional-DFTs-of-Real-Data
     */
    basics::complexVector DFT(const basics::Vector& vector) // REMEMBER: PACKED FORMAT!!
    {
      LOG("Finding REAL DFT of "+vector.name());

      basics::complexVector result("DFT of "+vector.name(),vector.length()/2+1,false);
      FFT_PFX(plan) plan = FFT_PFX(plan_dft_r2c_1d)(vector.length(),const_cast<Real*>(vector.data()),reinterpret_cast<FFT_PFX(complex*)>(result.data()),FFTW_ESTIMATE);
      assert( plan != NULL );

      FFT_PFX(execute)(plan);

      FFT_PFX(destroy_plan)(plan);

      return( result );
    }

    basics::complexMatrix DFT(const basics::complexMatrix& matrix) // columnwise
    {
      LOG("Finding DFT of the columns of "+matrix.name());

      basics::complexMatrix result("DFT of "+matrix.name(),matrix.rows(),matrix.cols(),false);
      int* n = new int[result.cols()];

      /* setup length-of-transforms array */
      for( int i=0;i<result.cols();i++ )
        n[i] = result.rows();

      FFT_PFX(plan) plan = FFT_PFX(plan_many_dft)(1,n,result.cols(),const_cast<FFT_PFX(complex*)>(reinterpret_cast<const FFT_PFX(complex*)>(matrix.data()[0])),NULL,1,result.rows(),reinterpret_cast<FFT_PFX(complex*)>(result.data()[0]),NULL,1,result.rows(),FFTW_FORWARD,FFTW_ESTIMATE);
      assert( plan != NULL );

      FFT_PFX(execute)(plan);

      FFT_PFX(destroy_plan)(plan);
      delete[] n;

      return( result );
    }

    basics::complexMatrices DFT2(const basics::complexMatrices& matrices)
    {
      LOG("Finding 2D DFT of "+matrices.name());

      basics::complexMatrices result("DFT of "+matrices.name(),matrices.rows(),matrices.cols(),matrices.matrices());

      FFT_PFX(iodim*) idim = new FFT_PFX(iodim)[2];
      idim[0].n = matrices.matrices();idim[0].is = idim[0].os = matrices.rows()*matrices.cols();
      idim[1].n = matrices.rows(); idim[1].is = idim[1].os = 1;
      FFT_PFX(iodim) hdim;
      hdim.n = matrices.cols(); hdim.is = hdim.os = matrices.rows();
      FFT_PFX(plan) plan = FFT_PFX(plan_guru_dft)(2,idim,1,&hdim,const_cast<FFT_PFX(complex*)>(reinterpret_cast<const FFT_PFX(complex*)>(matrices.data())),reinterpret_cast<FFT_PFX(complex*)>(result.data()),FFTW_FORWARD,FFTW_ESTIMATE);

      assert( plan != NULL );

      FFT_PFX(execute)(plan);

      FFT_PFX(destroy_plan)(plan);
      delete[] idim;

      return( result );
    }

    void DFTI(basics::complexMatrix& matrix) // columnwise, in-place
    {
      LOG("Finding DFT columns of "+matrix.name()+" (inplace)");

      int* n = new int[matrix.cols()];

      /* setup length-of-transforms array */
      for( int i=0;i<matrix.cols();i++ )
        n[i] = matrix.rows();

      FFT_PFX(plan) plan = FFT_PFX(plan_many_dft)(1,n,matrix.cols(),reinterpret_cast<FFT_PFX(complex*)>(matrix.data()[0]),NULL,1,matrix.rows(),reinterpret_cast<FFT_PFX(complex*)>(matrix.data()[0]),NULL,1,matrix.rows(),FFTW_FORWARD,FFTW_ESTIMATE);
      assert( plan != NULL );

      FFT_PFX(execute)(plan);

      FFT_PFX(destroy_plan)(plan);
      delete[] n;
    }

    void DFTTrunc(basics::complexMatrix& result, basics::complexMatrix& matrix, FFT_PFX(plan) plan)
    {
      if( plan )
        executeDFTI(plan,matrix);
      else
        DFTI(matrix);
      const int skipx = matrix.rows()-result.rows();
      for( int j=0;j<result.cols();++j ) { 
        int N = result.rows()/2;
        BLASCPFX(copy,BLASINT(N),BLASPCOMPLEX(matrix.data()[j]),BLASINT(mnlIntOne),BLASPCOMPLEX(result.data()[j]),BLASINT(mnlIntOne));
        BLASCPFX(copy,BLASINT(N),BLASPCOMPLEX(matrix.data()[j]+skipx+result.rows()/2),BLASINT(mnlIntOne),BLASPCOMPLEX(result.data()[j]+result.rows()/2),BLASINT(mnlIntOne));
      }
      result *= static_cast<Real>(result.rows())/static_cast<Real>(matrix.rows());
    }

    void DFT2I(basics::complexMatrices& matrices)
    {
      LOG("Finding 2D DFT columns of "+matrices.name()+" (inplace)");

      FFT_PFX(iodim*) idim = new FFT_PFX(iodim)[2];
      idim[0].n = matrices.matrices();idim[0].is = idim[0].os = matrices.rows()*matrices.cols();
      idim[1].n = matrices.rows(); idim[1].is = idim[1].os = 1;
      FFT_PFX(iodim) hdim;
      hdim.n = matrices.cols(); hdim.is = hdim.os = matrices.rows();
      FFT_PFX(plan) plan = FFT_PFX(plan_guru_dft)(2,idim,1,&hdim,reinterpret_cast<FFT_PFX(complex*)>(matrices.data()),reinterpret_cast<FFT_PFX(complex*)>(matrices.data()),FFTW_FORWARD,FFTW_ESTIMATE);

      assert( plan != NULL );

      FFT_PFX(execute)(plan);

      FFT_PFX(destroy_plan)(plan);
      delete[] idim;
    }

    void DFT2Trunc(basics::complexMatrices& result, basics::complexMatrices& matrices, FFT_PFX(plan) plan)
    {
      if( plan )
        executeDFT2I(plan,matrices);
      else
        DFT2I(matrices);
      int dKz = 0;
      const int skipx = matrices.rows()-result.rows();
      for( int l=0;l<matrices.matrices();++l ) {
        for( int j=0;j<matrices.cols();++j ) {
          int N = result.rows()/2;
          BLASCPFX(copy,BLASINT(N),BLASPCOMPLEX(matrices[l].data()[j]),BLASINT(mnlIntOne),BLASPCOMPLEX(result[l-dKz].data()[j]),BLASINT(mnlIntOne));
          BLASCPFX(copy,BLASINT(N),BLASPCOMPLEX(matrices[l].data()[j]+skipx+result.rows()/2),BLASINT(mnlIntOne),BLASPCOMPLEX(result[l-dKz].data()[j]+result.rows()/2),BLASINT(mnlIntOne));
        }
        if( l == result.matrices()/2 ) {
          dKz = matrices.matrices()-result.matrices();
          l += (matrices.matrices()-result.matrices());
        }
      }
      result *= static_cast<Real>(result.matrices())/static_cast<Real>(matrices.matrices())*static_cast<Real>(result.rows())/static_cast<Real>(matrices.rows());
    }

    /* TODO: plan_many! */
    basics::complexMatrix DFT(const basics::Matrix& matrix) // columnwise
    {
      LOG("Finding REAL DFT columns of "+matrix.name());

      basics::complexMatrix result("DFT of "+matrix.name(),matrix.rows()/2+1,matrix.cols(),false);
      for( int i=0;i<result.cols();i++ ) {
        FFT_PFX(plan) plan = FFT_PFX(plan_dft_r2c_1d)(result.rows(),const_cast<Real*>(matrix.data()[i]),const_cast<FFT_PFX(complex*)>(reinterpret_cast<const FFT_PFX(complex*)>(result.data()[i])),FFTW_ESTIMATE);
        FFT_PFX(execute)(plan);
        FFT_PFX(destroy_plan)(plan);
      }

      return( result );
    }

    /**
     *       Calculates the inverse DFT of an vector.
     * @param vector Input data
     * @return Inverse DFT of vector
     */
    basics::complexVector IDFT(const basics::complexVector& vector)
    {
      LOG("Finding IDFT of "+vector.name());

      basics::complexVector result("IDFT of "+vector.name(),vector.length(),false);
      FFT_PFX(complex*) out = (FFT_PFX(complex*))FFT_PFX(malloc)(sizeof(FFT_PFX(complex))*vector.length());

      FFT_PFX(plan) plan = FFT_PFX(plan_dft_1d)(vector.length(),reinterpret_cast<FFT_PFX(complex*)>(const_cast<std::complex<Real>*>(vector.data())),out,FFTW_BACKWARD,FFTW_ESTIMATE);
      assert( plan != NULL );

      FFT_PFX(execute)(plan);

      result = reinterpret_cast<std::complex<Real>*>(out);
      result /= static_cast<Real>(vector.length());

      FFT_PFX(destroy_plan)(plan);
      FFT_PFX(free)(out);

      return( result );
    }

    basics::complexMatrix IDFT(const basics::complexMatrix& matrix) // columnwise
    {
      LOG("Finding IDFT columns of "+matrix.name());

      basics::complexMatrix result("IDFT of "+matrix.name(),matrix.rows(),matrix.cols(),false);
      int *n = new int[result.rows()];

      /* setup length-of-transforms array */
      for( int i=0;i<result.rows();i++ ) 
        n[i] = result.rows();

      FFT_PFX(plan) plan = FFT_PFX(plan_many_dft)(1,n,result.cols(),const_cast<FFT_PFX(complex*)>(reinterpret_cast<const FFT_PFX(complex*)>(matrix.data()[0])),NULL,1,result.rows(),reinterpret_cast<FFT_PFX(complex*)>(result.data()[0]),NULL,1,result.rows(),FFTW_BACKWARD,FFTW_ESTIMATE);
      assert( plan != NULL );

      FFT_PFX(execute)(plan);

      FFT_PFX(destroy_plan)(plan);
      delete[] n;

      result *= static_cast<Real>(1.f)/static_cast<Real>(result.rows());

      return( result );
    }

    basics::complexMatrices IDFT2(const basics::complexMatrices& matrices)
    {
      LOG("Finding 2D IDFT of "+matrices.name());

      basics::complexMatrices result("IDFT of "+matrices.name(),matrices.rows(),matrices.cols(),matrices.matrices());

      FFT_PFX(iodim*) idim = new FFT_PFX(iodim)[2];
      idim[0].n = matrices.matrices();idim[0].is = idim[0].os = matrices.rows()*matrices.cols();
      idim[1].n = matrices.rows(); idim[1].is = idim[1].os = 1;
      FFT_PFX(iodim) hdim;
      hdim.n = matrices.cols(); hdim.is = hdim.os = matrices.rows();
      FFT_PFX(plan) plan = FFT_PFX(plan_guru_dft)(2,idim,1,&hdim,const_cast<FFT_PFX(complex*)>(reinterpret_cast<const FFT_PFX(complex*)>(matrices.data())),reinterpret_cast<FFT_PFX(complex*)>(result.data()),FFTW_BACKWARD,FFTW_ESTIMATE);

      assert( plan != NULL );

      FFT_PFX(execute)(plan);

      result *= static_cast<Real>(1.f)/static_cast<Real>(result.matrices()*result.rows());

      FFT_PFX(destroy_plan)(plan);
      delete[] idim;

      return( result );
    }

    void IDFTI(basics::complexMatrix& matrix) // columnwise, in-place
    {
      LOG("Finding IDFT columns of "+matrix.name()+" (inplace)");

      int* n = new int[matrix.cols()];

      /* setup length-of-transforms array */
      for( int i=0;i<matrix.cols();i++ )
        n[i] = matrix.rows();

      FFT_PFX(plan) plan = FFT_PFX(plan_many_dft)(1,n,matrix.cols(),reinterpret_cast<FFT_PFX(complex*)>(matrix.data()[0]),NULL,1,matrix.rows(),reinterpret_cast<FFT_PFX(complex*)>(matrix.data()[0]),NULL,1,matrix.rows(),FFTW_BACKWARD,FFTW_ESTIMATE);
      assert( plan != NULL );

      FFT_PFX(execute)(plan);

      matrix *= Real(1)/Real(matrix.rows());

      FFT_PFX(destroy_plan)(plan);
      delete[] n;
    }

    void IDFT2I(basics::complexMatrices& matrices)
    {
      LOG("Finding 2D IDFT of "+matrices.name()+" (inplace)");

      FFT_PFX(iodim*) idim = new FFT_PFX(iodim)[2];
      idim[0].n = matrices.matrices();idim[0].is = idim[0].os = matrices.rows()*matrices.cols();
      idim[1].n = matrices.rows(); idim[1].is = idim[1].os = 1;
      FFT_PFX(iodim) hdim;
      hdim.n = matrices.cols(); hdim.is = hdim.os = matrices.rows();
      FFT_PFX(plan) plan = FFT_PFX(plan_guru_dft)(2,idim,1,&hdim,reinterpret_cast<FFT_PFX(complex*)>(matrices.data()),reinterpret_cast<FFT_PFX(complex*)>(matrices.data()),FFTW_BACKWARD,FFTW_ESTIMATE);

      matrices *= static_cast<Real>(1.f)/static_cast<Real>(matrices.matrices()*matrices.rows());

      assert( plan != NULL );

      FFT_PFX(execute)(plan);

      FFT_PFX(destroy_plan)(plan);
      delete[] idim;
    }

    FFT_PFX(plan) planDFTI(basics::complexMatrix& matrix)
    {
      int* n = new int[matrix.cols()];

      /* setup length-of-transforms array */
      for( int i=0;i<matrix.cols();i++ )
        n[i] = matrix.rows();

      FFT_PFX(plan) result = FFT_PFX(plan_many_dft)(1,n,matrix.cols(),reinterpret_cast<FFT_PFX(complex*)>(matrix.data()[0]),NULL,1,matrix.rows(),reinterpret_cast<FFT_PFX(complex*)>(matrix.data()[0]),NULL,1,matrix.rows(),FFTW_FORWARD,FFTW_MEASURE);
      assert( result != NULL );

      delete[] n;

      return( result );
    }

    FFT_PFX(plan) planIDFTI(basics::complexMatrix& matrix) // remember: we won't be able to scale for you
    {
      int* n = new int[matrix.cols()];

      /* setup length-of-transforms array */
      for( int i=0;i<matrix.cols();i++ )
        n[i] = matrix.rows();

      FFT_PFX(plan) result = FFT_PFX(plan_many_dft)(1,n,matrix.cols(),reinterpret_cast<FFT_PFX(complex*)>(matrix.data()[0]),NULL,1,matrix.rows(),reinterpret_cast<FFT_PFX(complex*)>(matrix.data()[0]),NULL,1,matrix.rows(),FFTW_BACKWARD,FFTW_ESTIMATE);
      assert( result != NULL );

      delete[] n;

      return( result );
    }

    FFT_PFX(plan) planDFT2I(basics::complexMatrices& matrices)
    {
      FFT_PFX(iodim*) idim = new FFT_PFX(iodim)[2];
      idim[0].n = matrices.matrices();idim[0].is = idim[0].os = matrices.rows()*matrices.cols();
      idim[1].n = matrices.rows(); idim[1].is = idim[1].os = 1;
      FFT_PFX(iodim) hdim;
      hdim.n = matrices.cols(); hdim.is = hdim.os = matrices.rows();
      FFT_PFX(plan) result = FFT_PFX(plan_guru_dft)(2,idim,1,&hdim,reinterpret_cast<FFT_PFX(complex*)>(matrices.data()),reinterpret_cast<FFT_PFX(complex*)>(matrices.data()),FFTW_FORWARD,FFTW_ESTIMATE);
      assert( result != NULL );

      delete[] idim;

      return( result );
    }

    FFT_PFX(plan) planIDFT2I(basics::complexMatrices& matrices) // remember: we won't be able to scale for you
    {
      FFT_PFX(iodim*) idim = new FFT_PFX(iodim)[2];
      idim[0].n = matrices.matrices();idim[0].is = idim[0].os = matrices.rows()*matrices.cols();
      idim[1].n = matrices.rows(); idim[1].is = idim[1].os = 1;
      FFT_PFX(iodim) hdim;
      hdim.n = matrices.cols(); hdim.is = hdim.os = matrices.rows();
      FFT_PFX(plan) result = FFT_PFX(plan_guru_dft)(2,idim,1,&hdim,reinterpret_cast<FFT_PFX(complex*)>(matrices.data()),reinterpret_cast<FFT_PFX(complex*)>(matrices.data()),FFTW_BACKWARD,FFTW_ESTIMATE);
      assert( result != NULL );

      delete[] idim;

      return( result );
    }

    void executeDFTI(FFT_PFX(plan) plan, basics::complexMatrix& matrix)
    {
      FFT_PFX(execute_dft)(plan,reinterpret_cast<FFT_PFX(complex*)>(matrix.data()[0]),reinterpret_cast<FFT_PFX(complex*)>(matrix.data()[0]));
    }

    void executeDFT2I(FFT_PFX(plan) plan, basics::complexMatrices& matrices)
    {
      FFT_PFX(execute_dft)(plan,reinterpret_cast<FFT_PFX(complex*)>(matrices.data()),reinterpret_cast<FFT_PFX(complex*)>(matrices.data()));
    }

    void executeIDFTI(FFT_PFX(plan) plan, basics::complexMatrix& matrix)
    {
      FFT_PFX(execute_dft)(plan,reinterpret_cast<FFT_PFX(complex*)>(matrix.data()[0]),reinterpret_cast<FFT_PFX(complex*)>(matrix.data()[0]));
      matrix *= static_cast<Real>(1.f)/static_cast<Real>(matrix.rows());
    }

    void executeIDFT2I(FFT_PFX(plan) plan, basics::complexMatrices& matrices)
    {
      FFT_PFX(execute_dft)(plan,reinterpret_cast<FFT_PFX(complex*)>(matrices.data()),reinterpret_cast<FFT_PFX(complex*)>(matrices.data()));
      matrices *= static_cast<Real>(1.f)/static_cast<Real>(matrices.matrices()*matrices.rows());
    }
  }
}

