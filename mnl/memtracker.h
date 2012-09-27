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
#ifndef MEM_TRACKER_H_
#define MEM_TRACKER_H_

#include "config.h"

#include <vector>
#include <string>
#include <iostream>
#ifdef OPENMP
#include <omp.h>
#endif

namespace mnl {
  namespace basics {
    class Vector;
    class complexVector;
    class Matrix;
    class complexMatrix;
    class Matrices;
    class complexMatrices;
    template<class T> class Field2;
    template<class T> class Field3;
    typedef Field2<Matrix> Field2D;
    typedef Field2<complexMatrix> complexField2D;
    typedef Field3<Matrices> Field3D;
    typedef Field3<complexMatrices> complexField3D;
    class geometryStack;
    template<class T> class matrixStackT;
    typedef matrixStackT<basics::Matrix> matrixStack;
    class matricesStack;
  }
  namespace utilities {

#ifdef OPENMP
#define INIT_LOCK(l) omp_init_lock(l) 
#define DESTROY_LOCK(l) omp_destroy_lock(l)
#define SET_LOCK(l) omp_set_lock(l)
#define UNSET_LOCK(l) omp_unset_lock(l)
#else
#define INIT_LOCK(l)
#define DESTROY_LOCK(l)
#define SET_LOCK(l)
#define UNSET_LOCK(l)
#endif
#ifdef MNL_MEMORY_VERBOSE_
    class memTracker {
      public:
        memTracker();
        ~memTracker();

        enum chunkType {
          Vector 		= 1,
          Matrix 		= 2,
          Matrices 	= 3
        };


        int addChunk(chunkType type, unsigned long long size);
        void removeChunk(int id);

        std::string report();
      protected:
        struct memChunk {
          chunkType type;
          int id;
          unsigned long long size;
        };
        unsigned long long m_memory;
        std::vector<memChunk> m_chunk;
        int m_id;
#ifdef OPENMP
        omp_lock_t m_lock;
#endif
    };

    extern memTracker g_tracker;
#endif
    class bufferManager {
      public:
        bufferManager();
        ~bufferManager();

        basics::Vector* clone(const basics::Vector& input);
        basics::complexVector* clone(const basics::complexVector& input);
        basics::Matrix* clone(const basics::Matrix& input);
        basics::Field2D* clone(const basics::Field2D& input);
        basics::Field2D* cloneField(const basics::Matrix& input);
        basics::complexMatrix* clone(const basics::complexMatrix& input);
        basics::complexField2D* clone(const basics::complexField2D& input);
        basics::complexField2D* cloneField(const basics::complexMatrix& input);
        basics::Matrices* clone(const basics::Matrices& input);
        basics::Field3D* clone(const basics::Field3D& input);
        basics::Field3D* cloneField(const basics::Matrices& input);
        basics::complexMatrices* clone(const basics::complexMatrices& input);
        basics::complexField3D* clone(const basics::complexField3D& input);
        basics::complexField3D* cloneField(const basics::complexMatrices& input);
        basics::matrixStack* clone(const basics::matrixStack&);
        basics::Field2<basics::matrixStack>* clone(const basics::Field2<basics::matrixStack>&);
        basics::matricesStack* clone(const basics::matricesStack&);
        basics::Field3<basics::matricesStack>* clone(const basics::Field3<basics::matricesStack>&);

        basics::Vector* aquireVector(const std::string& name, int length);
        basics::complexVector* aquireComplexVector(const std::string& name, int length);
        basics::Matrix* aquireMatrix(const std::string& name, int rows, int cols);
        basics::Field2D* aquireField2D(const std::string& name, int rows, int cols);
        basics::Matrices* aquireMatrices(const std::string& name, int rows, int cols, int matrices);
        basics::Field3D* aquireField3D(const std::string& name, int rows, int cols, int matrices);

        basics::complexMatrix* aquireComplexMatrix(const std::string& name, int rows, int cols);
        basics::complexField2D* aquireComplexField2D(const std::string& name, int rows, int cols);
        basics::complexMatrices* aquireComplexMatrices(const std::string& name, int rows, int cols, int matrices);
        basics::complexField3D* aquireComplexField3D(const std::string& name, int rows, int cols, int matrices);

        basics::matrixStack* aquireMatrixStack(const std::string& name, const basics::geometryStack& geometry);
        basics::Field2<basics::matrixStack>* aquireMatrixStackField(const std::string& name, const basics::geometryStack& geometry);
        basics::matrixStack* aquireMatrixStack(const std::string& name, int rows, int cols, int size);
        basics::Field2<basics::matrixStack>* aquireMatrixStackField(const std::string& name, int rows, int cols, int size);
        basics::matricesStack* aquireMatricesStack(const std::string& name, const basics::geometryStack& geometry);
        basics::Field3<basics::matricesStack>* aquireMatricesStackField(const std::string& name, const basics::geometryStack& geometry);
        basics::matricesStack* aquireMatricesStack(const std::string& name, int rows, int cols, int matrices, int size);
        basics::Field3<basics::matricesStack>* aquireMatricesStackField(const std::string& name, int rows, int cols, int matrices, int size);

        void unlock(basics::Vector*&);
        void unlock(basics::complexVector*&);
        void unlock(basics::Matrix*&);
        void unlock(basics::Field2D*&);
        void unlock(basics::complexMatrix*&);
        void unlock(basics::complexField2D*&);
        void unlock(basics::Matrices*&);
        void unlock(basics::Field3D*&);
        void unlock(basics::complexMatrices*&);
        void unlock(basics::complexField3D*&);
        void unlock(basics::matrixStack*&);
        void unlock(basics::Field2<basics::matrixStack>*&);
        void unlock(basics::matricesStack*&);
        void unlock(basics::Field3<basics::matricesStack>*&);
      protected:
        struct bufferEntry {
          void* entry;
          bool locked;
        };
        std::vector<bufferEntry> m_vector;
        std::vector<bufferEntry> m_cvector;
        std::vector<bufferEntry> m_matrix;
        std::vector<bufferEntry> m_cmatrix;
        std::vector<bufferEntry> m_matrices;
        std::vector<bufferEntry> m_cmatrices;
        std::vector<bufferEntry> m_field2d;
        std::vector<bufferEntry> m_cfield2d;
        std::vector<bufferEntry> m_field3d;
        std::vector<bufferEntry> m_cfield3d;
        std::vector<bufferEntry> m_mstacks2d;
        std::vector<bufferEntry> m_mstacks3d;
#ifdef OPENMP
        omp_lock_t m_lock;
#endif
    };

    extern bufferManager g_manager;
  }
}

#endif

