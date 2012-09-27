#include "quadratic.h"
#include "sem.h"

#include "mnl/gll.h"
#include "mnl/hdf5.h"

#ifdef HAS_MPI
#include <mpi.h>
#endif

using namespace mnl;
using namespace std;

quadraticGeometry::quadraticGeometry(int N, int M, const basics::Vector& ggrid, const basics::Vector& weight) :
  geometryStack(ggrid,weight)
{
  HDF5::HDF5Reader reader(getOperatorFile());
  reader.read(*this,"grid");
  int max=m_grid.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i )
    m_grid[i]->getGH().computeJacobian();

  basics::Vector temp("temp",2);
  reader.read(temp,"size");

  m_size1 = temp[0];
  m_size2 = temp[1];

  computeMultiplicities();
  computeMass();
}

quadraticGeometry::~quadraticGeometry()
{
}

void quadraticGeometry::periodicDssum(basics::matrixStack& op) const
{
  dssum(op);
  int N = op[0].rows();
  for( int i=0;i<m_size1;++i ) { // bottom and top
    op[                    i][0  ] += op[i+(m_size2-1)*m_size1][N-1];
    op[i+(m_size2-1)*m_size1][N-1]  = op[                    i][  0];
  }
  for( int i=0;i<m_size2;++i ) { // left and right
    for( int l=0;l<N;++l ) {
      op[i*m_size1          ][l][0  ] += op[i*m_size1+m_size1-1][l][N-1];
      op[i*m_size1+m_size1-1][l][N-1]  = op[i*m_size1][l][0];
    }
  }
}

void quadraticGeometry::dssum(basics::matrixStack& op, int rank,
    int size, int tag, void* group) const
{
  int N = op[0].rows();
  for( int m=1;m<m_size1;++m ) // horizontal sweep
    for( int n=0;n<m_size2;++n )
      for( int l=0;l<N;++l ) {
        op[m+n*m_size1  ][l][0  ] += op[m+n*m_size1-1][l][N-1];
        op[m+n*m_size1-1][l][N-1]  = op[m+n*m_size1  ][l][0  ];
      }
  for( int m=0;m<m_size1;++m )
    for( int n=1;n<m_size2;++n )
      for( int l=0;l<N;++l ) {
        op[m+n*m_size1    ][0   ][l] += op[m+(n-1)*m_size1][N-1][l];
        op[m+(n-1)*m_size1][N-1 ][l]  = op[m+n*m_size1    ][0  ][l];
      }
}

void quadraticGeometry::mask(basics::matrixStack& op, int rank, int size) const
{
  int N = op[0].rows();
  for( int i=0;i<m_size1;++i ) { // bottom and top
    op[i][0].clear();
    op[i+(m_size2-1)*m_size1][N-1].clear();
  }
  for( int i=0;i<m_size2;++i ) { // left and right
    op[i*m_size1].clearRow(0);
    op[i*m_size1+m_size1-1].clearRow(N-1);
  }
}

void quadraticGeometry::mask(basics::matricesStack& op, int rank, int size) const
{
  op.at(0).clear();
  op.at(op[0].matrices()-1).clear();
  for( int l=1;l<op[0].matrices()-1;++l )
    mask(op.at(l),rank,size);
}

