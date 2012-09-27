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

#include "memtracker.h"
#include "vector.h"
#include "matrix.h"
#include "matrices.h"
#include "field.h"
#include "buffers.h"
#include "geometry.h"

#include <sstream>

using namespace std;

namespace mnl {
  namespace utilities {
    bufferManager g_manager;
#ifdef MNL_MEMORY_VERBOSE_
    memTracker g_tracker;

    memTracker::memTracker()
    {
      m_id = 0;
      m_memory = 0;
      cout << "init memtracker" << endl;
      INIT_LOCK(&m_lock);
    }

    memTracker::~memTracker()
    {
    }

    int memTracker::addChunk(chunkType type, unsigned long long size)
    {
      SET_LOCK(&m_lock);

      memChunk chunk;
      chunk.type = type;
      chunk.size = size;
      chunk.id = m_id++;
      m_chunk.push_back(chunk);
      m_memory += size;

      UNSET_LOCK(&m_lock);

      return chunk.id;
    }

    void memTracker::removeChunk(int id)
    {
      SET_LOCK(&m_lock);

      for( vector<memChunk>::iterator iter=m_chunk.begin();iter != m_chunk.end();++iter ) {
        if( iter->id  == id ) {
          m_memory -= iter->size;
          m_chunk.erase(iter);
          UNSET_LOCK(&m_lock);
          return;
        }	
      }

      UNSET_LOCK(&m_lock);
    }

    string memTracker::report()
    {
      stringstream stream;

      stream << "Total memory used: " << ((float)m_memory)/(1024.f*1024.f) << "MB" << endl;
      int vector=0, matrix=0, matrices=0;
      unsigned long long vectorSize=0, matrixSize=0, matricesSize=0;
      for( vector<memChunk>::iterator iter=m_chunk.begin();iter != m_chunk.end();++iter ) {
        if( iter->type == Vector ) {
          vector++;
          vectorSize += iter->size;
        }
        if( iter->type == Matrix ) {
          matrix++;
          matrixSize += iter->size;
        }
        if( iter->type == Matrices ) {
          matrices++;
          matricesSize += iter->size;
        }
      }
      stream << "	- " << vector << " vectors (" << ((float)(vectorSize)/(1024.f*1024.f)) << " MB)" << endl;
      stream << "	- " << matrix << " matrices (" << ((float)(matrixSize)/(1024.f*1024.f)) << " MB)" << endl;
      stream << "	- " << matrices << " matrices' (" << ((float)(matricesSize)/(1024.f*1024.f)) << " MB)" << endl;

      return stream.str();
    }
#endif

    bufferManager::bufferManager()
    {
      INIT_LOCK(&m_lock);
      m_matrices.clear();
    }

    bufferManager::~bufferManager()
    {
      SET_LOCK(&m_lock);
#ifdef DEBUG
      cout << "memtracker releasing " << m_vector.size() << " vectors at shutdown" << endl;
      cout << "memtracker releasing " << m_cvector.size() << " complex vectors at shutdown" << endl;
      cout << "memtracker releasing " << m_matrix.size() << " matrices at shutdown" << endl;
      cout << "memtracker releasing " << m_cmatrix.size() << " complex matrices at shutdown" << endl;
      cout << "memtracker releasing " << m_matrices.size() << " matrices' at shutdown" << endl;
      cout << "memtracker releasing " << m_cmatrices.size() << " complex matrices' at shutdown" << endl;
      cout << "memtracker releasing " << m_field2d.size() << " 2D fields at shutdown" << endl;
      cout << "memtracker releasing " << m_cfield2d.size() << " complex 2D fields at shutdown" << endl;
      cout << "memtracker releasing " << m_field3d.size() << " 3D fields at shutdown" << endl;
      cout << "memtracker releasing " << m_cfield3d.size() << " complex 3D fields at shutdown" << endl;
      cout << "memtracker releasing " << m_mstacks2d.size() << " 2D stacks at shutdown" << endl;
      cout << "memtracker releasing " << m_mstacks3d.size() << " 3D stacks at shutdown" << endl;
#endif

      for( vector<bufferEntry>::iterator it=m_vector.begin();it != m_vector.end();++it )
        delete (basics::Vector*)it->entry;

      for( vector<bufferEntry>::iterator it=m_cvector.begin();it != m_cvector.end();++it )
        delete (basics::complexVector*)it->entry;
      for( vector<bufferEntry>::iterator it=m_matrix.begin();it != m_matrix.end();++it )
        delete (basics::Matrix*)it->entry;
      for( vector<bufferEntry>::iterator it=m_cmatrix.begin();it != m_cmatrix.end();++it )
        delete (basics::complexMatrix*)it->entry;
      for( vector<bufferEntry>::iterator it=m_matrices.begin();it != m_matrices.end();++it )
        delete (basics::Matrices*)it->entry;
      for( vector<bufferEntry>::iterator it=m_cmatrices.begin();it != m_cmatrices.end();++it )
        delete (basics::complexMatrices*)it->entry;
      for( vector<bufferEntry>::iterator it=m_field2d.begin();it != m_field2d.end();++it )
        delete (basics::Field2D*)it->entry;
      for( vector<bufferEntry>::iterator it=m_cfield2d.begin();it != m_cfield2d.end();++it )
        delete (basics::complexField2D*)it->entry;
      for( vector<bufferEntry>::iterator it=m_field3d.begin();it != m_field3d.end();++it )
        delete (basics::Field3D*)it->entry;
      for( vector<bufferEntry>::iterator it=m_cfield3d.begin();it != m_cfield3d.end();++it )
        delete (basics::complexField3D*)it->entry;
      for( vector<bufferEntry>::iterator it=m_mstacks2d.begin();it != m_mstacks2d.end();++it )
        delete (basics::matrixStack*)it->entry;
      for( vector<bufferEntry>::iterator it=m_mstacks3d.begin();it != m_mstacks3d.end();++it )
        delete (basics::matricesStack*)it->entry;

      DESTROY_LOCK(&m_lock);
    }

    basics::Vector* bufferManager::clone(const basics::Vector& input)
    {
      return aquireVector("clone of "+input.name(),input.length());
    }

    basics::complexVector* bufferManager::clone(const basics::complexVector& input)
    {
      return aquireComplexVector("clone of "+input.name(),input.length());
    }

    basics::Matrix* bufferManager::clone(const basics::Matrix& input)
    {
      return aquireMatrix("clone of "+input.name(),input.rows(),input.cols());
    }

    basics::complexMatrix* bufferManager::clone(const basics::complexMatrix& input)
    {
      return aquireComplexMatrix("clone of "+input.name(),input.rows(),input.cols());
    }

    basics::Field2D* bufferManager::clone(const basics::Field2D& input)
    {
      return aquireField2D("clone of "+input.X().name().substr(0,input.X().name().size()-3),input.rows(),input.cols());
    }

    basics::Field2D* bufferManager::cloneField(const basics::Matrix& input)
    {
      return aquireField2D("clone of "+input.name().substr(0,input.name().size()-3),input.rows(),input.cols());
    }

    basics::complexField2D* bufferManager::clone(const basics::complexField2D& input)
    {
      return aquireComplexField2D("clone of "+input.X().name().substr(0,input.X().name().size()-3),input.rows(),input.cols());
    }

    basics::complexField2D* bufferManager::cloneField(const basics::complexMatrix& input)
    {
      return aquireComplexField2D("clone of "+input.name().substr(0,input.name().size()-3),input.rows(),input.cols());
    }

    basics::Matrices* bufferManager::clone(const basics::Matrices& input)
    {
      return aquireMatrices("clone of "+input.name(),input.rows(),input.cols(),input.matrices());
    }

    basics::complexMatrices* bufferManager::clone(const basics::complexMatrices& input)
    {
      return aquireComplexMatrices("clone of "+input.name(),input.rows(),input.cols(),input.matrices());
    }

    basics::Field3D* bufferManager::clone(const basics::Field3D& input)
    {
      return aquireField3D("clone of "+input.X().name().substr(0,input.X().name().size()-3),input.rows(),input.cols(),input.matrices());
    }

    basics::Field3D* bufferManager::cloneField(const basics::Matrices& input)
    {
      return aquireField3D("clone of "+input.name().substr(0,input.name().size()-3),input.rows(),input.cols(),input.matrices());
    }

    basics::complexField3D* bufferManager::clone(const basics::complexField3D& input)
    {
      return aquireComplexField3D("clone of "+input.X().name().substr(0,input.X().name().size()-3),input.rows(),input.cols(),input.matrices());
    }

    basics::complexField3D* bufferManager::cloneField(const basics::complexMatrices& input)
    {
      return aquireComplexField3D("clone of "+input.name().substr(0,input.name().size()-3),input.rows(),input.cols(),input.matrices());
    }

    basics::matrixStack* bufferManager::clone(const basics::matrixStack& input)
    {
      return aquireMatrixStack("clone of "+input.name(),input[0].rows(),
          input[0].cols(),input.size());
    }

    basics::Field2<basics::matrixStack>* bufferManager::clone(const basics::Field2<basics::matrixStack>& input)
    {
      return aquireMatrixStackField("clone of "+input.X().name(),input.X()[0].rows(),
          input.X()[0].cols(),input.X().size());
    }

    basics::matricesStack* bufferManager::clone(const basics::matricesStack& input)
    {
      return aquireMatricesStack("clone of "+input.name(),input[0].rows(),
          input[0].cols(),input[0].matrices(),input.size());
    }

    basics::Field3<basics::matricesStack>* bufferManager::clone(const basics::Field3<basics::matricesStack>& input)
    {
      return( aquireMatricesStackField("clone of "+input.X().name(),input.X()[0].rows(),
            input.X()[0].cols(),input.X()[0].matrices(),
            input.X().size()) );
    }

    basics::Vector* bufferManager::aquireVector(const string& name, int length)
    {
      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_vector.begin();it != m_vector.end();++it ) {
        if( !it->locked && ((basics::Vector*)it->entry)->length() == length ) {
          it->locked = true;
          basics::Vector* result = (basics::Vector*)it->entry;
          result->setName(name+" (managed)");
          UNSET_LOCK(&m_lock);
          return result;
        }
      }
      bufferEntry entry;
      entry.entry = (void*)new basics::Vector(name+" (managed)",length);
      entry.locked = true;
      m_vector.push_back(entry);

      UNSET_LOCK(&m_lock);

      return (basics::Vector*)entry.entry;
    }

    basics::complexVector* bufferManager::aquireComplexVector(const string& name, int length)
    {
      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_cvector.begin();it != m_cvector.end();++it ) {
        if( !it->locked && ((basics::complexVector*)it->entry)->length() == length ) {
          it->locked = true;
          basics::complexVector* result = (basics::complexVector*)it->entry;
          result->setName(name+" (managed)");
          UNSET_LOCK(&m_lock);
          return result;
        }
      }
      bufferEntry entry;
      entry.entry = (void*)new basics::complexVector(name+" (managed)",length);
      entry.locked = true;
      m_cvector.push_back(entry);

      UNSET_LOCK(&m_lock);

      return (basics::complexVector*)entry.entry;
    }

    basics::Matrix* bufferManager::aquireMatrix(const string& name, int rows, int cols)
    {
      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_matrix.begin();it != m_matrix.end();++it ) {
        if( !it->locked && ((basics::Matrix*)it->entry)->rows() == rows && ((basics::Matrix*)it->entry)->cols() == cols ) {
          it->locked = true;
          basics::Matrix* result = (basics::Matrix*)it->entry;
          result->setName(name+" (managed)");
          UNSET_LOCK(&m_lock);
          return result;
        }
      }
      bufferEntry entry;
      entry.entry = (void*)new basics::Matrix(name+" (managed)",rows,cols);
      entry.locked = true;
      m_matrix.push_back(entry);

      UNSET_LOCK(&m_lock);

      return (basics::Matrix*)entry.entry;
    }

    basics::Field2D* bufferManager::aquireField2D(const string& name, int rows, int cols)
    {
      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_field2d.begin();it != m_field2d.end();++it ) {
        if( !it->locked && ((basics::Field2D*)it->entry)->rows() == rows && ((basics::Field2D*)it->entry)->cols() == cols ) {
          it->locked = true;
          basics::Field2D* result = (basics::Field2D*)it->entry;
          result->setName(name+" (managed)");
          UNSET_LOCK(&m_lock);
          return result;
        }
      }
      bufferEntry entry;
      entry.entry = (void*)new basics::Field2D(name+" (managed)",rows,cols);
      entry.locked = true;
      m_field2d.push_back(entry);

      UNSET_LOCK(&m_lock);

      return (basics::Field2D*)entry.entry;
    }

    basics::Matrices* bufferManager::aquireMatrices(const string& name, int rows, int cols, int matrices)
    {
      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_matrices.begin();it != m_matrices.end();++it ) {
        if( !it->locked && ((basics::Matrices*)it->entry)->rows() == rows && ((basics::Matrices*)it->entry)->cols() == cols && ((basics::Matrices*)it->entry)->matrices() == matrices ) {
          it->locked = true;
          basics::Matrices* result = (basics::Matrices*)it->entry;
          result->setName(name+" (managed)");
          UNSET_LOCK(&m_lock);
          return result;
        }
      }
      bufferEntry entry;
      entry.entry = (void*)new basics::Matrices(name+" (managed)",rows,cols,matrices);
      entry.locked = true;
      m_matrices.push_back(entry);

      UNSET_LOCK(&m_lock);

      return (basics::Matrices*)entry.entry;
    }

    basics::Field3D* bufferManager::aquireField3D(const string& name, int rows, int cols, int matrices)
    {
      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_field3d.begin();it != m_field3d.end();++it ) {
        if( !it->locked && ((basics::Field3D*)it->entry)->rows() == rows && ((basics::Field3D*)it->entry)->cols() == cols ) {
          it->locked = true;
          basics::Field3D* result = (basics::Field3D*)it->entry;
          result->setName(name+" (managed)");
          UNSET_LOCK(&m_lock);
          return result;
        }
      }
      bufferEntry entry;
      entry.entry = (void*)new basics::Field3D(name+" (managed)",rows,cols,matrices);
      entry.locked = true;
      m_field3d.push_back(entry);

      UNSET_LOCK(&m_lock);

      return (basics::Field3D*)entry.entry;
    }

    basics::complexMatrix* bufferManager::aquireComplexMatrix(const string& name, int rows, int cols)
    {
      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_cmatrix.begin();it != m_cmatrix.end();++it ) {
        if( !it->locked && ((basics::complexMatrix*)it->entry)->rows() == rows && ((basics::complexMatrix*)it->entry)->cols() == cols ) {
          it->locked = true;
          basics::complexMatrix* result = (basics::complexMatrix*)it->entry;
          result->setName(name+" (managed)");
          UNSET_LOCK(&m_lock);
          return result;
        }
      }
      bufferEntry entry;
      entry.entry = (void*)new basics::complexMatrix(name+" (managed)",rows,cols);
      entry.locked = true;
      m_cmatrix.push_back(entry);

      UNSET_LOCK(&m_lock);

      return (basics::complexMatrix*)entry.entry;
    }

    basics::complexField2D* bufferManager::aquireComplexField2D(const string& name, int rows, int cols)
    {
      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_cfield2d.begin();it != m_cfield2d.end();++it ) {
        if( !it->locked && ((basics::complexField2D*)it->entry)->rows() == rows && ((basics::complexField2D*)it->entry)->cols() == cols ) {
          it->locked = true;
          basics::complexField2D* result = (basics::complexField2D*)it->entry;
          result->setName(name+" (managed)");
          UNSET_LOCK(&m_lock);
          return result;
        }
      }
      bufferEntry entry;
      entry.entry = (void*)new basics::complexField2D(name+" (managed)",rows,cols);
      entry.locked = true;
      m_cfield2d.push_back(entry);

      UNSET_LOCK(&m_lock);

      return (basics::complexField2D*)entry.entry;
    }

    basics::complexMatrices* bufferManager::aquireComplexMatrices(const string& name, int rows, int cols, int matrices)
    {
      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_cmatrices.begin();it != m_cmatrices.end();++it ) {
        if( !it->locked && ((basics::complexMatrices*)it->entry)->rows() == rows && ((basics::complexMatrices*)it->entry)->cols() == cols && ((basics::complexMatrices*)it->entry)->matrices() == matrices ) {
          it->locked = true;
          basics::complexMatrices* result = (basics::complexMatrices*)it->entry;
          result->setName(name+" (managed)");
          UNSET_LOCK(&m_lock);
          return result;
        }
      }
      bufferEntry entry;
      entry.entry = (void*)new basics::complexMatrices(name+" (managed)",rows,cols,matrices);
      entry.locked = true;
      m_cmatrices.push_back(entry);

      UNSET_LOCK(&m_lock);

      return (basics::complexMatrices*)entry.entry;
    }

    basics::complexField3D* bufferManager::aquireComplexField3D(const string& name, int rows, int cols, int matrices)
    {
      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_cfield3d.begin();it != m_cfield3d.end();++it ) {
        if( !it->locked && ((basics::complexField3D*)it->entry)->rows() == rows && ((basics::complexField3D*)it->entry)->cols() == cols ) {
          it->locked = true;
          basics::complexField3D* result = (basics::complexField3D*)it->entry;
          result->setName(name+" (managed)");
          UNSET_LOCK(&m_lock);
          return result;
        }
      }
      bufferEntry entry;
      entry.entry = (void*)new basics::complexField3D(name+" (managed)",rows,cols,matrices);
      entry.locked = true;
      m_cfield3d.push_back(entry);

      UNSET_LOCK(&m_lock);

      return (basics::complexField3D*)entry.entry;
    }

    basics::matrixStack* 
      bufferManager::aquireMatrixStack(const string& name, 
          const basics::geometryStack& geometry)
      {
        return aquireMatrixStack(name,geometry[0].getGH().getMapping().X().rows(),
            geometry[0].getGH().getMapping().X().cols(),
            geometry.size());
      }

    basics::matrixStack* 
      bufferManager::aquireMatrixStack(const string& name, 
          int rows, int cols, int size)
      {
        SET_LOCK(&m_lock);

        for( vector<bufferEntry>::iterator it=m_mstacks2d.begin();it != m_mstacks2d.end();++it ) {
          basics::matrixStack* foo = (basics::matrixStack*)it->entry;
          if( !it->locked && ((*foo)[0].rows() == rows && 
                (*foo)[0].cols() == cols &&
                foo->size()      == size   ) ) {
            it->locked = true;
            foo->setName(name);
            UNSET_LOCK(&m_lock);
            return foo;
          }
        }
        bufferEntry entry;
        basics::matrixStack* newop = new basics::matrixStack(name);
        basics::Matrix* foo2 = new basics::Matrix(name,rows,cols);
        newop->add(*foo2,size);
        entry.entry = (void*)newop;
        entry.locked = true;
        m_mstacks2d.push_back(entry);
        delete foo2;

        UNSET_LOCK(&m_lock);

        return newop;
      }

    basics::Field2<basics::matrixStack>* 
      bufferManager::aquireMatrixStackField(const string& name, 
          const basics::geometryStack& geometry)
      {
        return( new basics::Field2<basics::matrixStack>(aquireMatrixStack(name,geometry),
              aquireMatrixStack(name,geometry)) );
      }

    basics::Field2<basics::matrixStack>* 
      bufferManager::aquireMatrixStackField(const string& name, 
          int rows, int cols, int size)
      {
        return( new basics::Field2<basics::matrixStack>(aquireMatrixStack(name+" X",rows,cols,size),
              aquireMatrixStack(name+" Y",rows,cols,size)) );
      }

    basics::matricesStack* 
      bufferManager::aquireMatricesStack(const string& name, 
          const basics::geometryStack& geometry)
      {
        const basics::Matrix& X = geometry[0].getGH().getMapping().X();
        return aquireMatricesStack(name,X.rows(),X.cols(),
            X.rows(),geometry.size());
      }

    basics::matricesStack* 
      bufferManager::aquireMatricesStack(const string& name, int rows, int cols,
          int matrices, int size)
      {
        SET_LOCK(&m_lock);

        for( vector<bufferEntry>::iterator it=m_mstacks3d.begin();it != m_mstacks3d.end();++it ) {
          basics::matricesStack* foo = (basics::matricesStack*)it->entry;
          if( !it->locked && ((*foo)[0].rows()     == rows     && 
                (*foo)[0].cols()     == cols     &&
                (*foo)[0].matrices() == matrices &&
                foo->size()          == size )     ) {
            it->locked = true;
            foo->setName(name);
            UNSET_LOCK(&m_lock);
            return foo;
          }
        }
        bufferEntry entry;
        basics::matricesStack* newop = new basics::matricesStack(name,rows,cols,
            matrices,size);
        entry.entry = (void*)newop;
        entry.locked = true;
        m_mstacks3d.push_back(entry);

        UNSET_LOCK(&m_lock);

        return newop;
      }

    basics::Field3<basics::matricesStack>* 
      bufferManager::aquireMatricesStackField(const string& name, 
          const basics::geometryStack& geometry)
      {
        return( new basics::Field3<basics::matricesStack>(aquireMatricesStack(name+" X",geometry),
              aquireMatricesStack(name+" Z",geometry),
              aquireMatricesStack(name+" Y",geometry)) );
      }

    basics::Field3<basics::matricesStack>* 
      bufferManager::aquireMatricesStackField(const string& name, int rows, int cols,
          int matrices, int size)
      {
        return( new basics::Field3<basics::matricesStack>(aquireMatricesStack(name+" X",rows,cols,matrices,size),
              aquireMatricesStack(name+" Z",rows,cols,matrices,size),
              aquireMatricesStack(name+" Y",rows,cols,matrices,size)) );
      }

    void bufferManager::unlock(basics::Vector*& buffer)
    {
      if( buffer == NULL )
        return;

      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_vector.begin();it != m_vector.end();++it ) {
        if( it->entry == buffer ) {
          it->locked = false;
          buffer = NULL;
          UNSET_LOCK(&m_lock);
          return;
        }
      }

      LOG("unmanaged vector unlocked");
      assert( false ); // unmanaged matrix 
    }


    void bufferManager::unlock(basics::Matrix*& buffer)
    {
      if( buffer == NULL )
        return;

      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_matrix.begin();it != m_matrix.end();++it ) {
        if( it->entry == buffer ) {
          it->locked = false;
          buffer = NULL;
          UNSET_LOCK(&m_lock);
          return;
        }
      }

      LOG("unmanaged matrix unlocked");
      assert( false ); // unmanaged matrix 
    }

    void bufferManager::unlock(basics::complexMatrix*& buffer)
    {
      if( buffer == NULL )
        return;

      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_cmatrix.begin();it != m_cmatrix.end();++it ) {
        if( it->entry == buffer ) {
          it->locked = false;
          buffer = NULL;
          UNSET_LOCK(&m_lock);
          return;
        }
      }

      LOG("unmanaged matrix unlocked");
      assert( false ); // unmanaged matrix 
    }

    void bufferManager::unlock(basics::Matrices*& buffer)
    {
      if( buffer == NULL )
        return;

      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_matrices.begin();it != m_matrices.end();++it ) {
        if( it->entry == buffer ) {
          it->locked = false;
          buffer = NULL;
          UNSET_LOCK(&m_lock);
          return;
        }
      }

      LOG("unmanaged matrices unlocked");
      assert( false ); // unmanaged matrices 
    }

    void bufferManager::unlock(basics::complexMatrices*& buffer)
    {
      if( buffer == NULL )
        return;

      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_cmatrices.begin();it != m_cmatrices.end();++it ) {
        if( it->entry == buffer ) {
          it->locked = false;
          buffer = NULL;
          UNSET_LOCK(&m_lock);
          return;
        }
      }

      LOG("unmanaged matrices unlocked");
      assert( false ); // unmanaged matrices 
    }

    void bufferManager::unlock(basics::Field2D*& buffer)
    {
      if( buffer == NULL )
        return;

      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_field2d.begin();it != m_field2d.end();++it ) {
        if( it->entry == buffer ) {
          it->locked = false;
          buffer = NULL;
          UNSET_LOCK(&m_lock);
          return;
        }
      }

      LOG("unmanaged field2d unlocked");
      assert( false ); // unmanaged matrices 
    }

    void bufferManager::unlock(basics::complexField2D*& buffer)
    {
      if( buffer == NULL )
        return;

      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_cfield2d.begin();it != m_cfield2d.end();++it ) {
        if( it->entry == buffer ) {
          it->locked = false;
          buffer = NULL;
          UNSET_LOCK(&m_lock);
          return;
        }
      }

      LOG("unmanaged complexfield2d unlocked");
      assert( false ); // unmanaged matrices 
    }

    void bufferManager::unlock(basics::Field3D*& buffer)
    {
      if( buffer == NULL )
        return;

      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_field3d.begin();it != m_field3d.end();++it ) {
        if( it->entry == buffer ) {
          it->locked = false;
          buffer = NULL;
          UNSET_LOCK(&m_lock);
          return;
        }
      }

      LOG("unmanaged field3d unlocked");
      assert( false ); // unmanaged matrices 
    }

    void bufferManager::unlock(basics::complexField3D*& buffer)
    {
      if( buffer == NULL )
        return;

      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_cfield3d.begin();it != m_cfield3d.end();++it ) {
        if( it->entry == buffer ) {
          it->locked = false;
          buffer = NULL;
          UNSET_LOCK(&m_lock);
          return;
        }
      }

      LOG("unmanaged complexfield3d unlocked");
      assert( false );
    }

    void bufferManager::unlock(basics::matrixStack*& buffer)
    {
      if( buffer == NULL )
        return;

      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_mstacks2d.begin();it != m_mstacks2d.end();++it ) {
        if( it->entry == buffer ) {
          it->locked = false;
          buffer = NULL;
          UNSET_LOCK(&m_lock);
          return;
        }
      }

      LOG("unmanaged matrixstack unlocked");
      assert( false ); // unmanaged matrices 
    }

    void bufferManager::unlock(basics::Field2<basics::matrixStack>*& buffer)
    {
      if( buffer == NULL )
        return;

      basics::matrixStack* fx = &buffer->X();
      basics::matrixStack* fy = &buffer->Y();
      unlock(fx);
      unlock(fy);

      delete buffer;
      buffer = NULL;
    }

    void bufferManager::unlock(basics::matricesStack*& buffer)
    {
      if( buffer == NULL )
        return;

      SET_LOCK(&m_lock);

      for( vector<bufferEntry>::iterator it=m_mstacks3d.begin();it != m_mstacks3d.end();++it ) {
        if( it->entry == buffer ) {
          it->locked = false;
          buffer = NULL;
          UNSET_LOCK(&m_lock);
          return;
        }
      }

      LOG("unmanaged matricesstack unlocked");
      assert( false ); // unmanaged matrices 
    }

    void bufferManager::unlock(basics::Field3<basics::matricesStack>*& buffer)
    {
      if( !buffer )
        return;

      basics::matricesStack* fx = &buffer->X();
      basics::matricesStack* fy = &buffer->Y();
      basics::matricesStack* fz = &buffer->Z();
      unlock(fx);
      unlock(fy);
      unlock(fz);

      delete buffer;
      buffer = NULL;
    }
  }
}
