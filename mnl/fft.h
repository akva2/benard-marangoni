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
#ifndef MNL_FFTW_H_
#define MNL_FFTW_H_

#include <stdio.h>

#include "config.h"
#include "vector.h"
#include "matrix.h"
#include "matrices.h"
#include "field.h"
#include "function.h"

#ifdef MNL_USE_DOUBLE_
#define FFT_PFX(name) mnl::fftw::FFTW_MANGLE_DOUBLE(name)
#else
#define FFT_PFX(name) mnl::fftw::FFTW_MANGLE_FLOAT(name)
#endif

namespace mnl {
  /* Namespace for everything related to finding Fourier transforms. Uses fftw for the calculations. */
  namespace fftw {
#include <fftw3.h>

    /* Calculates the DFT of a vector. */
    basics::complexVector DFT(const basics::complexVector& vector);

    /* Calculates the DFT of a REAL vector. Returns transform in a packed format. */
    basics::complexVector DFT(const basics::Vector& vector);

    /* Calculates the DFT of the columns of a matrix. */
    basics::complexMatrix DFT(const basics::complexMatrix& matrix);

    /* Calculates the DFT of the columns of a REAL matrix. Returns transform in a packed format. */
    basics::complexMatrix DFT(const basics::Matrix& matrix);

    /* Calculates the 2D DFT along matrices and columns. */
    basics::complexMatrices DFT2(const basics::complexMatrices&);

    /* Calculates the inverse DFT of a vector. */
    basics::complexVector IDFT(const basics::complexVector& vector);

    /* Calculates the inverse DFT of the columns of a matrix. */
    basics::complexMatrix IDFT(const basics::complexMatrix& matrix);

    /* Calculates the inverse 2D DFT along matrices and columns. */
    basics::complexMatrices IDFT2(const basics::complexMatrices&);

    /* Calculates the DFT of the columns of a matrix in-place. */
    void DFTI(basics::complexMatrix& matrix);
    inline void DFTI(basics::complexField2D& field)
    {
      DFTI(field.X());
      DFTI(field.Y());
    }

    /* Calculates the DFT of matrix, truncating into result. */
    void DFTTrunc(basics::complexMatrix& result, basics::complexMatrix& matrices, FFT_PFX(plan) plan=NULL);
    inline void DFTTrunc(basics::complexField2D& result, basics::complexField2D& field, FFT_PFX(plan) plan=NULL)
    {
      DFTTrunc(result.X(),field.X(),plan);
      DFTTrunc(result.Y(),field.Y(),plan);
    }

    /* Calculates the 2D DFT of matrices and columns in-place. */
    void DFT2I(basics::complexMatrices& matrices);
    inline void DFT2I(basics::complexField3D& field)
    {
      DFT2I(field.X());
      DFT2I(field.Y());
      DFT2I(field.Z());
    }

    /* Calculates the 2D DFT of matrices, truncating into result. */
    void DFT2Trunc(basics::complexMatrices& result, basics::complexMatrices& matrices, FFT_PFX(plan) plan=NULL);
    inline void DFT2Trunc(basics::complexField3D& result, basics::complexField3D& field, FFT_PFX(plan) plan=NULL)
    {
      DFT2Trunc(result.X(),field.X(),plan);
      DFT2Trunc(result.Y(),field.Y(),plan);
      DFT2Trunc(result.Z(),field.Z(),plan);
    }

    /* Calculates the inverse DFT of the columns of a matrix in-place. */
    void IDFTI(basics::complexMatrix& matrix); // In place!
    inline void IDFTI(basics::complexField2D& field)
    {
      IDFTI(field.X());
      IDFTI(field.Y());
    }

    /* Calculates the inverse 2D DFT of matrices and columns in-place. */
    void IDFT2I(basics::complexMatrices& matrices); // In place!
    inline void IDFT2I(basics::complexField3D& field)
    {
      IDFT2I(field.X());
      IDFT2I(field.Y());
      IDFT2I(field.Z());
    }

    FFT_PFX(plan) planDFTI(basics::complexMatrix& matrix);
    FFT_PFX(plan) planIDFTI(basics::complexMatrix& matrix);

    FFT_PFX(plan) planDFT2I(basics::complexMatrices& matrix);
    FFT_PFX(plan) planIDFT2I(basics::complexMatrices& matrix);

    void executeDFTI(FFT_PFX(plan) plan, basics::complexMatrix& matrices);
    inline void executeDFTI(FFT_PFX(plan) plan, basics::complexField2D& result)
    {
      executeDFTI(plan,result.X());
      executeDFTI(plan,result.Y());
    }
    void executeDFT2I(FFT_PFX(plan) plan, basics::complexMatrices& matrices);
    void executeIDFTI(FFT_PFX(plan) plan, basics::complexMatrix& matrices);
    inline void executeDFT2I(FFT_PFX(plan) plan, basics::complexField3D& result)
    {
      executeDFT2I(plan,result.X());
      executeDFT2I(plan,result.Y());
      executeDFT2I(plan,result.Z());
    }
    inline void executeIDFTI(FFT_PFX(plan) plan, basics::complexField2D& result)
    {
      executeIDFTI(plan,result.X());
      executeIDFTI(plan,result.Y());
    }
    void executeIDFT2I(FFT_PFX(plan) plan, basics::complexMatrices& matrices);
    inline void executeIDFT2I(FFT_PFX(plan) plan, basics::complexField3D& result)
    {
      executeIDFT2I(plan,result.X());
      executeIDFT2I(plan,result.Y());
      executeIDFT2I(plan,result.Z());
    }

    class DoDFT : public basics::Functor1<mnl::basics::complexMatrix> {
      public:
        DoDFT()
        {
          m_plan = NULL;
        }

        virtual ~DoDFT()
        {
          FFT_PFX(destroy_plan)(m_plan);
        }

        void init(const mnl::basics::complexMatrix& mat)
        {
          m_plan = planDFTI(const_cast<mnl::basics::complexMatrix&>(mat));
        }

        virtual void operator()(mnl::basics::complexMatrix& mat)
        {
          if( m_plan )
            executeDFTI(m_plan,mat);
          else
            DFTI(mat);
        }
      protected:
        FFT_PFX(plan) m_plan;
    };

    class DoIDFT : public basics::Functor1<mnl::basics::complexMatrix> {
      public:
        DoIDFT()
        {
          m_plan = NULL;
        }

        virtual ~DoIDFT()
        {
          if( m_plan )
            FFT_PFX(destroy_plan)(m_plan);
        }

        virtual void init(const mnl::basics::complexMatrix& mat)
        {
          m_plan = planIDFTI(const_cast<mnl::basics::complexMatrix&>(mat));
        }

        virtual void operator()(mnl::basics::complexMatrix& mat)
        {
          if( m_plan )
            executeIDFTI(m_plan,mat);
          else
            IDFTI(mat);
        }
      protected:
        FFT_PFX(plan) m_plan;
    };
  }
}

#endif

