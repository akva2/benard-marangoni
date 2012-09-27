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

#ifndef _MNL_HDF5_H_
#define _MNL_HDF5_H_

#ifndef __sun__

#include "config.h"
#include "vector.h"
#include "matrix.h"
#include "matrices.h"
#include "field.h"
#include "geometry.h"

#include "H5Cpp.h"

#include <string>

namespace mnl {
  namespace HDF5 {
    class HDF5Reader {
      public:
        HDF5Reader(const std::string& fName);
        ~HDF5Reader();

        bool read(void* where, const H5::CompType& type,
            const std::string& datasetname);
        bool read(basics::Vector& where,          const std::string& datasetName);
        bool read(basics::complexVector& where,   const std::string& datasetName);
        bool read(basics::Matrix& where, const std::string& datasetName);
        bool read(basics::Field2<basics::Matrix>& where,
            const std::string& datasetName);
        bool read(basics::complexMatrix& where,   const std::string& datasetName);
        bool read(basics::Field2<basics::complexMatrix>& where, 
            const std::string& datasetName);
        bool read(basics::complexMatrices& where, const std::string& datasetName);
        bool read(basics::Matrices& where,        const std::string& datasetName);
        bool read(basics::geometryStack& where,   const std::string& datasetName);
        bool read(basics::geometryStack3D& where, const std::string& datasetName);
        bool read(basics::matrixStack& where,     const std::string& datasetName);
        bool read(basics::matricesStack& where,   const std::string& datasetName);
      private:
        H5::H5File* h5File;
    };

    class HDF5Writer {
      public:
        HDF5Writer(const std::string& fName);
        ~HDF5Writer();

        void Close();

        H5::DataSpace createDataSpace(const int n_Nx, const int n_Ny=0, const int n_Nz=0, const int size=0);
        H5::DataSpace createDataSpace(const basics::Vector& vector);
        H5::DataSpace createDataSpace(const basics::complexVector& vector);
        H5::DataSpace createDataSpace(const basics::Matrix& matrix);
        H5::DataSpace createDataSpace(const basics::complexMatrix& matrix);
        H5::DataSpace createDataSpace(const basics::Matrices& matrices);
        H5::DataSpace createDataSpace(const basics::complexMatrices& matrices);

        /* these add to a dataspace return by createDataSpace()  */
        void add(const H5::DataSpace& dataSpace, H5::CompType set, void* data, const std::string& name = "");
        void add(const H5::DataSpace& dataSpace, const basics::Vector& vector, const std::string& name = "");
        void add(const H5::DataSpace& dataSpace, const basics::complexVector& vector, const std::string& name = "");
        void add(const H5::DataSpace& dataSpace, const basics::Matrix& matrix, const std::string& name = "");
        void add(const H5::DataSpace& dataSpace, const basics::complexMatrix& matrix, const std::string& name = "");
        void add(const H5::DataSpace& dataSpace, const basics::Matrices& matrices, const std::string& name = "");
        void add(const H5::DataSpace& dataSpace, const basics::complexMatrices& matrices, const std::string& name = "");
        void add(const H5::DataSpace& dataSpace, const basics::matricesStack& stack, const std::string& name = "");

        /* these add to its own dataspace */
        void add(const H5::CompType& set, void* data, const int N, const std::string& name = "");
        void add(const basics::Vector& vector, const std::string& name = "");
        void add(const basics::complexVector& vector, const std::string& name = "");
        void add(const basics::Matrix& matrix, const std::string& name = "");
        void add(const basics::Field2<basics::Matrix>& field, const std::string& name = "");
        void add(const basics::complexMatrix& matrix, const std::string& name = "");
        void add(const basics::Field2<basics::complexMatrix>& field, const std::string& name = "");
        void add(const basics::Matrices& matrices, const std::string& name = "");
        void add(const basics::complexMatrices& matrices, const std::string& name="");
        void add(const basics::matricesStack& stack, const std::string& name="");
      private:
        H5::H5File* h5File;
        unsigned int index;
    };
  }
}

#endif

#endif

