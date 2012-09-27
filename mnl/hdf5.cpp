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

#ifndef __sun__

#include "hdf5.h"

#include <sstream>

using namespace H5;
using namespace std;
namespace mnl {
  namespace HDF5 {
    HDF5Reader::HDF5Reader(const string& fName) 
    {
      H5::Exception::dontPrint();
      try {
        h5File = new H5File(fName,H5F_ACC_RDONLY);
      } catch( FileIException error ) {
        error.printError();
        /* do something drastic */
      }
    }

    HDF5Reader::~HDF5Reader()
    {
      delete h5File;
    }

    bool HDF5Reader::read(void* where, const CompType& type, const string& datasetName)
    {
      try {
        DataSet dataSet = h5File->openDataSet(datasetName);
        DataSpace dataSpace = dataSet.getSpace();
        dataSet.read(where,type,dataSpace);
        return( true );
      } catch( DataSetIException error ) {
        LOG("Failure during H5 dataset open");
        return( false );
      }
      catch( DataSpaceIException error ) {
        LOG("H5 Dataset has no dataspace");
        return( false );
      }
      catch( FileIException error ) {
        LOG("Failure during H5 I/O");
        return( false );
      }
    }

    bool HDF5Reader::read(basics::Vector& where, const string& datasetName)
    {
      try {
        DataSet dataSet = h5File->openDataSet(datasetName);
        DataSpace dataSpace = dataSet.getSpace();

#ifndef NDEBUG
        assert( dataSpace.getSimpleExtentNdims() == 1 );
        hsize_t dims;
        dataSpace.getSimpleExtentDims(&dims);
        assert( dims == where.length() );
#endif

#ifdef MNL_USE_DOUBLE_
        dataSet.read(where.data(),PredType::NATIVE_DOUBLE,dataSpace);
#else
        dataSet.read(where.data(),PredType::NATIVE_FLOAT,dataSpace);
#endif
        return( true );
      } catch( FileIException ) {
        LOG("Failure during H5 I/O");
        return( false );
      }
      catch( DataSetIException error ) {
        LOG("Failure during H5 dataset open");
        return( false );
      }
      catch( DataSpaceIException error ) {
        LOG("H5 Dataset has no dataspace");
        return( false );
      }
    }

    bool HDF5Reader::read(basics::complexVector& where, const string& datasetName)
    {
      try {
        DataSet dataSet = h5File->openDataSet(datasetName);
        DataSpace dataSpace = dataSet.getSpace();

#ifndef NDEBUG
        assert( dataSpace.getSimpleExtentNdims() == 1 );
        hsize_t dims;
        dataSpace.getSimpleExtentDims(&dims);
        assert( dims == 2*where.length() );
#endif

#ifdef MNL_USE_DOUBLE_
        dataSet.read(where.data(),PredType::NATIVE_DOUBLE,dataSpace);
#else
        dataSet.read(where.data(),PredType::NATIVE_FLOAT,dataSpace);
#endif
        return( true );
      } catch( FileIException ) {
        LOG("Failure during H5 I/O");
        return( false );
      }
      catch( DataSetIException error ) {
        LOG("Failure during H5 dataset open");
        return( false );
      }
      catch( DataSpaceIException error ) {
        LOG("H5 Dataset has no dataspace");
        return( false );
      }
    }

    bool HDF5Reader::read(basics::Matrix& where, const string& datasetName)
    {
      try {
        DataSet dataSet = h5File->openDataSet(datasetName);
        DataSpace dataSpace = dataSet.getSpace();

#ifndef NDEBUG
        assert( dataSpace.getSimpleExtentNdims() == 2 );
        hsize_t dims[2];
        dataSpace.getSimpleExtentDims(dims);
        assert( (dims[0] == where.rows()) && (dims[1] == where.cols()) );
#endif

#ifdef MNL_USE_DOUBLE_
        dataSet.read(where.data()[0],PredType::NATIVE_DOUBLE,dataSpace);
#else
        dataSet.read(where.data()[0],PredType::NATIVE_FLOAT,dataSpace);
#endif
        return( true );
      } catch( FileIException ) {
        LOG("Failure during H5 I/O");
        return( false );
      }
      catch( DataSetIException error ) {
        LOG("Failure during H5 dataset open");
        return( false );
      }
      catch( DataSpaceIException error ) {
        LOG("H5 Dataset has no dataspace");
        return( false );
      }
    }

    bool HDF5Reader::read(basics::Field2<basics::Matrix>& where, const string& datasetName)
    {
      bool bOkay;

      bOkay  = read(where.X(),datasetName+"-X");
      bOkay &= read(where.Y(),datasetName+"-Y");

      return( bOkay );
    }

    bool HDF5Reader::read(basics::complexMatrix& where, const string& datasetName)
    {
      try {
        DataSet dataSet = h5File->openDataSet(datasetName);
        DataSpace dataSpace = dataSet.getSpace();

#ifndef NDEBUG
        assert( dataSpace.getSimpleExtentNdims() == 2 );
        hsize_t dims[2];
        dataSpace.getSimpleExtentDims(dims);
        assert( (dims[0] == where.rows()) && (dims[1] == 2*where.cols()) );
#endif

#ifdef MNL_USE_DOUBLE_
        dataSet.read(where.data()[0],PredType::NATIVE_DOUBLE,dataSpace);
#else
        dataSet.read(where.data()[0],PredType::NATIVE_FLOAT,dataSpace);
#endif
        return( true );
      } catch( FileIException ) {
        LOG("Failure during H5 I/O");
        return( false );
      }
      catch( DataSetIException error ) {
        LOG("Failure during H5 dataset open");
        return( false );
      }
      catch( DataSpaceIException error ) {
        LOG("H5 Dataset has no dataspace");
        return( false );
      }
    }

    bool HDF5Reader::read(basics::Field2<basics::complexMatrix>& where,
        const string& datasetName)
    {
      bool bOkay;

      bOkay  = read(where.X(), datasetName+"-X");
      bOkay &= read(where.Y(), datasetName+"-Y");

      return( bOkay );
    }

    bool HDF5Reader::read(basics::Matrices& where, const string& datasetName)
    {
      try {
        DataSet dataSet = h5File->openDataSet(datasetName);
        DataSpace dataSpace = dataSet.getSpace();

#ifndef NDEBUG
        assert( dataSpace.getSimpleExtentNdims() == 3 );
        hsize_t dims[3];
        dataSpace.getSimpleExtentDims(dims);
        //assert( (dims[0] == where.matrices()) && (dims[1] == where.cols()) && (dims[2] == where.rows()) );
#endif

#ifdef MNL_USE_DOUBLE_
        dataSet.read(where.data(),PredType::NATIVE_DOUBLE,dataSpace);
#else
        dataSet.read(where.data(),PredType::NATIVE_FLOAT,dataSpace);
#endif
        return( true );
      } catch( FileIException ) {
        LOG("Failure during H5 I/O");
        return( false );
      }
      catch( DataSetIException error ) {
        LOG("Failure during H5 dataset open");
        return( false );
      }
      catch( DataSpaceIException error ) {
        LOG("H5 Dataset has no dataspace");
        return( false );
      }
    }

    bool HDF5Reader::read(basics::complexMatrices& where,
        const string& datasetName)
    {
      try {
        DataSet dataSet = h5File->openDataSet(datasetName);
        DataSpace dataSpace = dataSet.getSpace();

#ifndef NDEBUG
        assert( dataSpace.getSimpleExtentNdims() == 3 );
        hsize_t dims[3];
        dataSpace.getSimpleExtentDims(dims);
        assert( (dims[0] == where.matrices()) && (dims[1] == where.cols()*2) && (dims[2] == where.rows()) );
#endif

#ifdef MNL_USE_DOUBLE_
        dataSet.read(where.data(),PredType::NATIVE_DOUBLE,dataSpace);
#else
        dataSet.read(where.data(),PredType::NATIVE_FLOAT,dataSpace);
#endif
        return( true );
      } catch( FileIException ) {
        LOG("Failure during H5 I/O");
        return( false );
      }
      catch( DataSetIException error ) {
        LOG("Failure during H5 dataset open");
        return( false );
      }
      catch( DataSpaceIException error ) {
        LOG("H5 Dataset has no dataspace");
        return( false );
      }
    }

    bool HDF5Reader::read(basics::geometryStack& where, const string& datasetName)
    {
      try {
        DataSet dataSetX = h5File->openDataSet(datasetName+".x");
        DataSet dataSetY = h5File->openDataSet(datasetName+".y");
        DataSpace dataSpaceX = dataSetX.getSpace();
        DataSpace dataSpaceY = dataSetY.getSpace();

        assert( dataSpaceX.getSimpleExtentNdims() == 3 );
        hsize_t dims[3];
        dataSpaceX.getSimpleExtentDims(dims);
        DataSpace memspace(2, dims+1, NULL);
        for( int i=0;i<dims[0];++i ) {
          basics::spectralElement2D* foo = new basics::spectralElement2D(dims[1]-1,dims[2]-1);
          hsize_t sel[3];
          sel[0]=i; sel[1] = sel[2] = 0;
          hsize_t cnt[3];
          cnt[0] = 1; cnt[1] = dims[1]; cnt[2] = dims[2];
          dataSpaceX.selectHyperslab(H5S_SELECT_SET,cnt,sel);
          dataSetX.read(foo->getGH().getMapping().X().data()[0],PredType::NATIVE_DOUBLE,memspace,dataSpaceX);
          dataSpaceY.selectHyperslab(H5S_SELECT_SET,cnt,sel);
          dataSetY.read(foo->getGH().getMapping().Y().data()[0],PredType::NATIVE_DOUBLE,memspace,dataSpaceY);
          where.getGridVector().push_back(foo);
        }

        return( true );
      } catch( FileIException ) {
        LOG("Failure during H5 I/O");
        return( false );
      }
      catch( DataSetIException error ) {
        LOG("Failure during H5 dataset open");
        return( false );
      }
      catch( DataSpaceIException error ) {
        LOG("H5 Dataset has no dataspace");
        return( false );
      }
    }

    bool HDF5Reader::read(basics::geometryStack3D& where,
        const string& datasetName)
    {
      try {
        DataSet dataSetX = h5File->openDataSet(datasetName+"3.x");
        DataSet dataSetY = h5File->openDataSet(datasetName+"3.y");
        DataSet dataSetZ = h5File->openDataSet(datasetName+"3.z");
        DataSpace dataSpaceX = dataSetX.getSpace();
        DataSpace dataSpaceY = dataSetY.getSpace();
        DataSpace dataSpaceZ = dataSetZ.getSpace();

        assert( dataSpaceX.getSimpleExtentNdims() == 4 );
        hsize_t dims[4];
        dataSpaceX.getSimpleExtentDims(dims);
        DataSpace memspace(3, dims+1, NULL);
        for( int i=0;i<dims[0];++i ) {
          basics::spectralElement3D* foo = new basics::spectralElement3D(dims[1]-1,dims[2]-1,dims[3]-1);
          hsize_t sel[4];
          sel[0]=i; sel[1] = sel[2] = sel[3] = 0;
          hsize_t cnt[4];
          cnt[0] = 1; cnt[1] = dims[1]; cnt[2] = dims[2]; cnt[3] = dims[3];
          dataSpaceX.selectHyperslab(H5S_SELECT_SET,cnt,sel);
          dataSetX.read(foo->getGH().getMapping().X().data(),PredType::NATIVE_DOUBLE,memspace,dataSpaceX);
          dataSpaceY.selectHyperslab(H5S_SELECT_SET,cnt,sel);
          dataSetY.read(foo->getGH().getMapping().Y().data(),PredType::NATIVE_DOUBLE,memspace,dataSpaceY);
          dataSpaceZ.selectHyperslab(H5S_SELECT_SET,cnt,sel);
          dataSetZ.read(foo->getGH().getMapping().Z().data(),PredType::NATIVE_DOUBLE,memspace,dataSpaceZ);
          where.getGridVector().push_back(foo);
        }

        return( true );
      } catch( FileIException ) {
        LOG("Failure during H5 I/O");
        return( false );
      }
      catch( DataSetIException error ) {
        LOG("Failure during H5 dataset open");
        return( false );
      }
      catch( DataSpaceIException error ) {
        LOG("H5 Dataset has no dataspace");
        return( false );
      }
    }

    bool HDF5Reader::read(basics::matrixStack& where, const string& datasetName)
    {
      try {
        DataSet dataSet = h5File->openDataSet(datasetName);
        DataSpace dataSpace = dataSet.getSpace();

        assert( dataSpace.getSimpleExtentNdims() == 3 );
        hsize_t dims[3];
        dataSpace.getSimpleExtentDims(dims);
        DataSpace memspace(2, dims+1, NULL);
        for( int i=0;i<dims[0];++i ) {
          hsize_t sel[3];
          sel[0] = i; sel[1] = sel[2] = 0;
          hsize_t cnt[3];
          cnt[0] = 1; cnt[1] = dims[1]; cnt[2] = dims[2];
          dataSpace.selectHyperslab(H5S_SELECT_SET,cnt,sel);
          dataSet.read(where[i].data()[0],PredType::NATIVE_DOUBLE,memspace,dataSpace);
        }

        return( true );
      } catch( FileIException ) {
        LOG("Failure during H5 I/O");
        return( false );
      }
      catch( DataSetIException error ) {
        LOG("Failure during H5 dataset open");
        return( false );
      }
      catch( DataSpaceIException error ) {
        LOG("H5 Dataset has no dataspace");
        return( false );
      }
    }

    bool HDF5Reader::read(basics::matricesStack& where, const string& datasetName)
    {
      try {
        DataSet dataSet = h5File->openDataSet(datasetName);
        DataSpace dataSpace = dataSet.getSpace();

        assert( dataSpace.getSimpleExtentNdims() == 4 );
        hsize_t dims[4];
        dataSpace.getSimpleExtentDims(dims);
        DataSpace memspace(3, dims+1, NULL);
        for( int i=0;i<dims[0];++i ) {
          hsize_t sel[4];
          sel[0] = i; sel[1] = sel[2] = sel[3] = 0;
          hsize_t cnt[4];
          cnt[0] = 1; cnt[1] = dims[1]; cnt[2] = dims[2]; cnt[3] = dims[3];
          dataSpace.selectHyperslab(H5S_SELECT_SET,cnt,sel);
          dataSet.read(where[i].data(),PredType::NATIVE_DOUBLE,memspace,dataSpace);
        }
        return( true );
      } catch( FileIException ) {
        LOG("Failure during H5 I/O");
        return( false );
      }
      catch( DataSetIException error ) {
        LOG("Failure during H5 dataset open");
        return( false );
      }
      catch( DataSpaceIException error ) {
        LOG("H5 Dataset has no dataspace");
        return( false );
      }
    }

    HDF5Writer::HDF5Writer(const string& fName) : index(0)
    {
      H5::Exception::dontPrint();
      try {
        h5File = new H5File(fName,H5F_ACC_TRUNC);
      } catch( FileIException error ) {
        error.printError();
        /* do something drastic */
      }
    }

    HDF5Writer::~HDF5Writer()
    {
      Close();
    }

    void HDF5Writer::Close()
    {
      delete h5File;
      h5File = NULL;
    }

    DataSpace HDF5Writer::createDataSpace(const int n_Nx, const int n_Ny, const int n_Nz, int size)
    {
      int rank=1;
      if( n_Ny )
        rank = 2;
      if( n_Nz )
        rank = 3;
      if( size )
        rank = 4;

      hsize_t* dim = new hsize_t[rank];
      dim[0] = n_Nx;
      if( n_Ny )
        dim[1] = n_Ny;
      if( n_Nz )
        dim[2] = n_Nz;
      if( size )
        dim[3] = size;

      DataSpace space(rank,dim);

      delete[] dim;

      return( space );
    }

    DataSpace HDF5Writer::createDataSpace(const basics::Vector& vector)
    {
      return( createDataSpace(vector.length()) );
    }

    DataSpace HDF5Writer::createDataSpace(const basics::complexVector& vector)
    {
      return( createDataSpace(2*vector.length()) );
    }

    DataSpace HDF5Writer::createDataSpace(const basics::Matrix& matrix)
    {
      return( createDataSpace(matrix.rows(),matrix.cols()) );
    }

    DataSpace HDF5Writer::createDataSpace(const basics::complexMatrix& matrix)
    {
      return( createDataSpace(matrix.rows(),2*matrix.cols()) );
    }

    DataSpace HDF5Writer::createDataSpace(const basics::Matrices& matrices)
    {
      return( createDataSpace(matrices.matrices(),matrices.cols(),matrices.rows()) );
    }

    DataSpace HDF5Writer::createDataSpace(const basics::complexMatrices& matrices)
    {
      return( createDataSpace(matrices.matrices(),matrices.cols()*2,matrices.rows()) );
    }

    void HDF5Writer::add(const DataSpace& dataSpace, CompType set, void* data, const string& name)
    {
      DataSet dataSet = h5File->createDataSet(name,set,dataSpace);
      dataSet.write(data,set);
    }

    void HDF5Writer::add(const DataSpace& dataSpace, const basics::Vector& vector, const string& name)
    {
#ifndef NDEBUG
      assert( dataSpace.getSimpleExtentNdims() == 1 );
      hsize_t size;
      dataSpace.getSimpleExtentDims(&size);
      assert( size == vector.length() );
#endif

      string uname = name;
      if( name.empty() )
        uname = vector.name();

#ifdef MNL_USE_DOUBLE_
      DataSet dataSet = h5File->createDataSet(uname,PredType::NATIVE_DOUBLE,dataSpace);
      dataSet.write(vector.data(),PredType::NATIVE_DOUBLE);
#else
      DataSet dataSet = h5File->createDataSet(uname,PredType::NATIVE_DOUBLE,dataSpace);
      dataSet.write(vector.data(),PredType::NATIVE_DOUBLE);
#endif
    }

    void HDF5Writer::add(const DataSpace& dataSpace, const basics::complexVector& vector, const string& name)
    {
#ifndef NDEBUG
      assert( dataSpace.getSimpleExtentNdims() == 1 );
      hsize_t size;
      dataSpace.getSimpleExtentDims(&size);
      assert( size == 2*vector.length() );
#endif

      string uname = name;
      if( name.empty() )
        uname = vector.name();

#ifdef MNL_USE_DOUBLE_
      DataSet dataSet = h5File->createDataSet(uname,PredType::NATIVE_DOUBLE,dataSpace);
      dataSet.write(vector.data(),PredType::NATIVE_DOUBLE);
#else
      DataSet dataSet = h5File->createDataSet(uname,PredType::NATIVE_DOUBLE,dataSpace);
      dataSet.write(vector.data(),PredType::NATIVE_DOUBLE);
#endif
    }

    void HDF5Writer::add(const DataSpace& dataSpace, const basics::Matrix& matrix, const string& name)
    {
#ifndef NDEBUG
      assert( dataSpace.getSimpleExtentNdims() == 2 );
      hsize_t size[2];
      dataSpace.getSimpleExtentDims(size);
      assert( (size[0] == matrix.rows()) && (size[1] == matrix.cols()) );
#endif

      string uname = name;
      if( name.empty() )
        uname = matrix.name();

#ifdef MNL_USE_DOUBLE_
      DataSet dataSet = h5File->createDataSet(uname,PredType::NATIVE_DOUBLE,dataSpace);
      dataSet.write(matrix.data()[0],PredType::NATIVE_DOUBLE);
#else
      DataSet dataSet = h5File->createDataSet(uname,PredType::NATIVE_DOUBLE,dataSpace);
      dataSet.write(matrix.data()[0],PredType::NATIVE_DOUBLE);
#endif
    }

    void HDF5Writer::add(const DataSpace& dataSpace, const basics::complexMatrix& matrix, const string& name)
    {
#ifndef NDEBUG
      assert( dataSpace.getSimpleExtentNdims() == 2 );
      hsize_t size[2];
      dataSpace.getSimpleExtentDims(size);
      assert( (size[0] == matrix.rows()) && (size[1] == 2*matrix.cols()) );
#endif

      string uname = name;
      if( name.empty() )
        uname = matrix.name();

#ifdef MNL_USE_DOUBLE_
      DataSet dataSet = h5File->createDataSet(uname,PredType::NATIVE_DOUBLE,dataSpace);
      dataSet.write(matrix.data()[0],PredType::NATIVE_DOUBLE);
#else
      DataSet dataSet = h5File->createDataSet(uname,PredType::NATIVE_DOUBLE,dataSpace);
      dataSet.write(matrix.data()[0],PredType::NATIVE_DOUBLE);
#endif
    }

    void HDF5Writer::add(const DataSpace& dataSpace, const basics::Matrices& matrices, const string& name)
    {
#ifndef NDEBUG
      assert( dataSpace.getSimpleExtentNdims() == 3 );
      hsize_t size[3];
      dataSpace.getSimpleExtentDims(size);
      assert( (size[0] == matrices.matrices()) && (size[1] == matrices.cols()) && (size[2] == matrices.rows()) );
#endif

      string uname = name;
      if( name.empty() )
        uname = matrices.name();

#ifdef MNL_USE_DOUBLE_
      DataSet dataSet = h5File->createDataSet(uname,PredType::NATIVE_DOUBLE,dataSpace);
      dataSet.write(matrices.data(),PredType::NATIVE_DOUBLE);
#else
      DataSet dataSet = h5File->createDataSet(uname,PredType::NATIVE_DOUBLE,dataSpace);
      dataSet.write(matrices.data(),PredType::NATIVE_DOUBLE);
#endif
    }

    void HDF5Writer::add(const DataSpace& dataSpace, const basics::complexMatrices& matrices, const string& name)
    {
#ifndef NDEBUG
      assert( dataSpace.getSimpleExtentNdims() == 3 );
      hsize_t size[3];
      dataSpace.getSimpleExtentDims(size);
      assert( (size[0] == matrices.matrices()) && (size[1] == matrices.cols()*2) && (size[2] == matrices.rows()) );
#endif

      string uname = name;
      if( name.empty() )
        uname = matrices.name();

#ifdef MNL_USE_DOUBLE_
      DataSet dataSet = h5File->createDataSet(uname,PredType::NATIVE_DOUBLE,dataSpace);
      dataSet.write(matrices.data(),PredType::NATIVE_DOUBLE);
#else
      DataSet dataSet = h5File->createDataSet(uname,PredType::NATIVE_DOUBLE,dataSpace);
      dataSet.write(matrices.data(),PredType::NATIVE_DOUBLE);
#endif
    }

    void HDF5Writer::add(const DataSpace& dataSpace, const basics::matricesStack& stack, const string& name)
    {
#ifndef NDEBUG
      assert( dataSpace.getSimpleExtentNdims() == 4 );
      hsize_t size[4];
      dataSpace.getSimpleExtentDims(size);
      assert( (size[3] == stack[0].matrices()) && (size[1] == stack[0].cols()) && (size[2] == stack[0].rows()) && (size[0] == stack.size()) );
#endif

      string uname = name;
      if( name.empty() )
        uname = stack.name();

      Real* foo = new double[stack.size()*stack[0].length()];
      for( int i=0;i<stack.size();++i )
        memcpy(foo+i*stack[0].length(),stack[i].data(),stack[i].length()*sizeof(Real));

#ifdef MNL_USE_DOUBLE_
      DataSet dataSet = h5File->createDataSet(uname,PredType::NATIVE_DOUBLE,dataSpace);
      dataSet.write(foo,PredType::NATIVE_DOUBLE);
#else
      DataSet dataSet = h5File->createDataSet(uname,PredType::NATIVE_FLOAT,dataSpace);
      dataSet.write(foo,PredType::NATIVE_FLOAT);
#endif
      delete[] foo;
    }

    void HDF5Writer::add(const CompType& set, void* data, const int N, const string& name)
    {
      add(createDataSpace(1),set,data,name);
    }

    void HDF5Writer::add(const basics::Vector& vector, const string& name)
    {
      add(createDataSpace(vector),vector,name);
    }

    void HDF5Writer::add(const basics::complexVector& vector, const string& name)
    {
      add(createDataSpace(vector),vector,name);
    }

    void HDF5Writer::add(const basics::Matrix& matrix, const string& name)
    {
      add(createDataSpace(matrix),matrix,name);
    }

    void HDF5Writer::add(const basics::Field2<basics::Matrix>& field, const string& name)
    {
      add(createDataSpace(field.X()),field.X(),name+"-X");
      add(createDataSpace(field.Y()),field.Y(),name+"-Y");
    }

    void HDF5Writer::add(const basics::complexMatrix& matrix, const string& name)
    {
      add(createDataSpace(matrix),matrix,name);
    }

    void HDF5Writer::add(const basics::Field2<basics::complexMatrix>& field, const string& name)
    {
      add(createDataSpace(field.X()),field.X(),name+"-X");
      add(createDataSpace(field.Y()),field.Y(),name+"-Y");
    }

    void HDF5Writer::add(const basics::Matrices& matrices, const string& name)
    {
      add(createDataSpace(matrices),matrices,name);
    }

    void HDF5Writer::add(const basics::complexMatrices& matrices, const string& name)
    {
      add(createDataSpace(matrices),matrices,name);
    }

    void HDF5Writer::add(const basics::matricesStack& stack, const string& name)
    {
      add(createDataSpace(stack.size(),stack[0].rows(),
            stack[0].cols(),stack[0].matrices()),stack,name);
    }
  }
}
#endif

