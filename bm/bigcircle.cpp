#include "bigcircle.h"
#include "sem.h"

#include "mnl/gll.h"

#ifdef HAS_MPI
#include <mpi.h>
#endif

#define TAG_INTERVAL	64

using namespace mnl;
using namespace std;

bigCircleGeometry::bigCircleGeometry(int N, int M, const basics::Vector& ggrid, const basics::Vector& weight) :
  geometryStack(ggrid,weight)
{
  HDF5::HDF5Reader reader("bigcircle.hdf5");
  reader.read(*this,"grid");
  int max=m_grid.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i )
    m_grid[i]->getGH().computeJacobian();

  m_division = (int)log2(size()/12); 
  m_size1 = (int)pow(2.f,m_division/2);
  m_size2 = m_size1*m_size1;
  /* g1 */
  if( m_division == 0 ) {
    m_size1 = m_size2 = 1;
    utilities::Range range;
    range.size(1);
    range[0] = 0;
    for( int i=0;i<4;++i ) {
      m_range.push_back(range);
    }
  } else {
    m_range.push_back(utilities::Range::colon(m_size1*(m_size1-1),m_size1*m_size1-1));
    /* g2 */
    m_range.push_back(utilities::Range::colon(m_size1-1,m_size1*m_size1-1,m_size1));
    /* g3 */
    m_range.push_back(utilities::Range::colon(0,m_size1-1));
    /* g4 */
    m_range.push_back(utilities::Range::colon(0,m_size1*(m_size1-1),m_size1));
  }

  computeMultiplicities();
  computeMass();
}

bigCircleGeometry::~bigCircleGeometry()
{
}

void bigCircleGeometry::innerDssum(basics::matrixStack& element, int group) const
{
  int N = element[0].rows();
  int skew = group*m_size2;
  for( int j=1;j<=m_size1-1;++j ) {
    for( int i=0;i<m_range[2].size();++i) {
      int elem11 = m_range[2][i]+skew+j*m_size1;
      int elem12 = m_range[2][i]+skew+(j-1)*m_size1;
      element[elem11][N-1] += element[elem12][0];
      element[elem12][0]    = element[elem11][N-1];
    }
  }
  for( int j=1;j<=m_size1-1;++j ) {
    for( int i=0;i<m_range[3].size();++i) {
      int elem21 = m_range[3][i]+skew+j;
      int elem22 = m_range[3][i]+skew+(j-1);
      for( int l=0;l<N;++l ) {
        element[elem21][l][0] += element[elem22][l][N-1];
        element[elem22][l][N-1] = element[elem21][l][0];
      }
    }
  }
}

void bigCircleGeometry::dssum_row_to_row(basics::matrixStack& element, int group1, int group2,
    int edge1, int edge2, int row1, int row2, int skip) const
{
  int skew1 = group1*m_size2;
  int skew2 = group2*m_size2;
  int i=0;
  int N = element[0].rows();
  if( skip == 1 ) {
    for( int l=0;l<element[0].rows()-1;++l ) {
      element[m_range[edge1][i]+skew1][l][row1] += element[m_range[edge2][i]+skew2][l][row2];
      element[m_range[edge2][i]+skew2][l][row2]  = element[m_range[edge1][i]+skew1][l][row1];
    }
    skip = 0;
    i++;
  }
  for( i;i<m_range[0].size()+skip;++i ) {
    for( int l=0;l<N;++l ) {
      element[m_range[edge1][i]+skew1][l][row1] += element[m_range[edge2][i]+skew2][l][row2];
      element[m_range[edge2][i]+skew2][l][row2]  = element[m_range[edge1][i]+skew1][l][row1];
    }
  }
  if( skip == -1 ) {
    for( int l=0;l<element[0].rows()-1;++l ) {
      element[m_range[edge1][i]+skew1][l+1][row1] += element[m_range[edge2][i]+skew2][l+1][row2];
      element[m_range[edge2][i]+skew2][l+1][row2]  = element[m_range[edge1][i]+skew1][l+1][row1];
    }
  }
}

void bigCircleGeometry::dssum_row_to_rowMPI(basics::matrixStack& element, 
    Real* temp, char* temp2,
    int group, int edge, int row,
    int& pos, bool dosum) const
{
  int skew = group*m_size2;
  int N = element[0].rows();
  int len = N*m_range[0].size();
  int N2 = N;

  if( dosum ) {
#ifdef HAS_MPI
    MPI_Unpack(temp2,8*len*sizeof(Real),&pos,temp,len,MPI_DOUBLE,MPI_COMM_WORLD);
#endif
    /* store data back in place */
    for( int i=0;i<m_range[0].size();++i )
      AXPY_ROW_FROM_COL(element[m_range[edge][i]+skew].data()[0]+row,temp+i*N);
  } else {
    for( int i=0;i<m_range[0].size();++i ) {
      COPY_COL_FROM_ROW(temp+i*N,element[m_range[edge][i]+skew].data()[0]+row);
    }
#ifdef HAS_MPI
    MPI_Pack(temp,len,MPI_DOUBLE,temp2,8*len*sizeof(Real),&pos,MPI_COMM_WORLD);
#endif
  }
}

void bigCircleGeometry::dssum_col_to_colMPI(basics::matrixStack& element, 
    Real* temp, char* temp2,
    int group, int edge, int col,
    int& pos, bool dosum) const
{
  int skew = group*m_size2;
  int N = element[0].rows();
  int len = N*m_range[0].size();

  int N2 = N;
  if( dosum ) {
#ifdef HAS_MPI
    MPI_Unpack(temp2,8*len*sizeof(Real),&pos,temp,len,MPI_DOUBLE,MPI_COMM_WORLD);
#endif
    /* store data back in place */
    for( int i=0;i<m_range[0].size();++i )
      AXPY_COL_FROM_COL(element[m_range[edge][i]+skew].data()[col],temp+i*N);
  } else {
    for( int i=0;i<m_range[0].size();++i )
      COPY_COL_FROM_COL(temp+i*N,element[m_range[edge][i]+skew].data()[col]);
#ifdef HAS_MPI
    MPI_Pack(temp,len,MPI_DOUBLE,temp2,8*len*sizeof(Real),&pos,MPI_COMM_WORLD);
#endif
  }
}

void bigCircleGeometry::dssum_col_to_col(basics::matrixStack& element, int group1, int group2,
    int edge1, int edge2, int col1, int col2, int skip) const
{
  int skew1 = group1*m_size2;
  int skew2 = group2*m_size2;
  int i=0;
  if( skip == 1 ) {
    int N2 = element[0].rows()-1;
    for( int l=0;l<N2;++l ) {
      element[m_range[edge1][i]+skew1][col1][l] += element[m_range[edge2][i]+skew2][col2][l];
      element[m_range[edge2][i]+skew2][col2][l] = element[m_range[edge1][i]+skew1][col1][l];
    }
    skip = 0;
    i++;
  }
  for( i;i<m_range[0].size()+skip;++i ) {
    element[m_range[edge1][i]+skew1][col1] +=  
      element[m_range[edge2][i]+skew2][col2];
    element[m_range[edge2][i]+skew2][col2] =
      element[m_range[edge1][i]+skew1][col1];
  }
  if( skip == -1 ) {
    int N2 = element[0].rows()-1;
    for( int l=0;l<N2;++l ) {
      element[m_range[edge1][i]+skew1][col1][l+1] += element[m_range[edge2][i]+skew2][col2][l+1];
      element[m_range[edge2][i]+skew2][col2][l+1] = element[m_range[edge1][i]+skew1][col1][l+1];
    }
  }
}

void bigCircleGeometry::dssum_col_to_row(basics::matrixStack& element, int group1, 
    int group2, int edge1, int edge2, int col1, int row2) const
{
  int skew1 = group1*m_size2;
  int skew2 = group2*m_size2;
  int N = element[0].rows();
  for( int i=0;i<m_range[0].size();++i ) {
    for( int l=0;l<N;++l ) {
      element[m_range[edge1][i]+skew1][col1][l] += element[m_range[edge2][m_size1-1-i]+skew2][l][row2];
      element[m_range[edge2][m_size1-1-i]+skew2][l][row2] = element[m_range[edge1][i]+skew1][col1][l];
    }
    element[m_range[edge2][m_size1-1-i]+skew2].copyRow(row2,
        element[m_range[edge1][i]+skew1][col1]);
  }
}

void bigCircleGeometry::dssum_row_to_colb(basics::matrixStack& element, int group1, 
    int group2, int edge1, int edge2, int row1, int col2) const
{
  int skew1 = group1*m_size2;
  int skew2 = group2*m_size2;
  int N = element[0].rows();
  for( int i=0;i<m_range[0].size();++i ) {
    int elem1 = m_range[edge1][i]+skew1;
    int elem2 = m_range[edge2][i]+skew2;
    for( int k=0;k<element[elem1].rows();++k ) {
      element[elem1][k][row1]    += element[elem2][col2][N-1-k];
      element[elem2][col2][N-1-k] = element[elem1][k][row1];
    }
  }
}

void bigCircleGeometry::periodicDssum(basics::matrixStack& element) const
{
  dssum(element);
  int N = element[0].rows();

  dssum_col_to_col(element, 4, 9,0,2,0,N-1,0); /* 5/g1 <-> 10/g3 */
  dssum_col_to_col(element, 5, 8,0,2,0,N-1,0); /* 6/g1 <->  9/g3 */

  dssum_row_to_row(element, 6,11,1,3,N-1,0,0); /*  7/g2 <-> 12/g4 */
  dssum_row_to_row(element, 7,10,1,3,N-1,0,0); /*  8/g2 <-> 11/g4 */
}

void bigCircleGeometry::dssum(basics::matrixStack& element, int rank,
    int size, int tag, void* group) const
{
  tag *= TAG_INTERVAL;

  int N = element[0].rows();
  /* sum within each subblock */
  if( m_division > 1 )
    for( int i=0;i<12/size;++i )
      innerDssum(element,i);

  Real* temp;  
  char* temp2;
  char* temp3;
  int pos=0;
  if( size > 1 ) {
    temp  = new Real[m_range[0].size()*element[0].rows()];
    temp2 = new char[m_range[0].size()*element[0].rows()*sizeof(Real)*8];
    temp3 = new char[m_range[0].size()*element[0].rows()*sizeof(Real)*8];
  }

  if( size == 1 ) {
    /* sum along the block boundaries */
    dssum_row_to_row (element, 4, 5,1,3,N-1,  0   ); //  5/g2 <->  6/g4
    dssum_row_to_colb(element, 5, 6,1,0,N-1,  0   ); //  6/g2 >-<  7/g1
    dssum_col_to_col (element, 6, 7,2,0,N-1,  0   ); //  7/g3 <->  8/g1
    dssum_col_to_row (element, 7, 8,2,1,N-1,N-1   ); //  8/g3 <->  9/g2
    dssum_row_to_row (element, 8, 9,3,1,  0,N-1   ); //  9/g4 <-> 10/g2
    dssum_row_to_colb(element, 9,10,3,2,  0,N-1   ); // 10/g4 >-< 11/g3
    dssum_col_to_col (element,10,11,0,2,  0,N-1   ); // 11/g1 <-> 12/g3
    dssum_col_to_row (element,11, 4,0,3,  0,  0   ); // 12/g1 <-> 5/g4
    dssum_row_to_row (element,0, 1,1,3,N-1,  0   ); //  1/g2 <->  2/g4
    dssum_row_to_row (element,2, 3,3,1,  0,N-1   ); //  3/g4 <->  4/g2
    dssum_col_to_col (element,0, 3,2,0,N-1,  0   ); //  1/g3 <->  4/g1
    dssum_col_to_col (element,1, 2,2,0,N-1,  0   ); //  2/g3 <->  3/g1

    /* lower left */
    dssum_col_to_col (element,0, 4,0,2,  0,N-1   ); //  1/g1 <->  5/g3
    dssum_row_to_row (element,0,11,3,1,  0,N-1,-1); //  1/g4 <-> 12/g2
    element[m_range[0][m_range[0].size()-1]+11*m_size2][  0][N-1] =
      element[m_range[0][                  0]           ][  0][  0];
    /* lower right */
    dssum_col_to_col (element,1, 5,0,2,  0,N-1   ); //  2/g1 <->  6/g3
    dssum_row_to_row (element,1, 6,1,3,N-1,  0,-1); //  2/g2 <->  7/g4
    element[m_range[0][                  0]+6*m_size2][  0][  0] = 
      element[m_range[0][m_range[0].size()-1]+  m_size2][  0][N-1];
    /* top right */
    dssum_col_to_col (element,2, 8,2,0,N-1,  0   ); //  3/g3 <->  9/g1
    dssum_row_to_row (element,2, 7,1,3,N-1,  0, 1); //  3/g2 <->  8/g4
    element[m_range[2][0]+7*m_size2][N-1][  0] =
      element[m_range[1][0]+2*m_size2][N-1][N-1];
    /* top left */
    dssum_col_to_col (element,3, 9,2,0,N-1,  0   ); //  4/g3 <-> 10/g1
    dssum_row_to_row (element,3,10,3,1,  0,N-1, 1); //  4/g4 <-> 11/g2
    element[m_range[1][  0]+10*m_size2][N-1][N-1] =
      element[m_range[2][  0]+ 3*m_size2][N-1][  0];
  } else if( size == 2 ) {
    /* sum within each subblock */
    if( rank == 0 ) { // 1 4 5 10 11 12
      /* sum along the block boundaries */
      int pos=0;
      dssum_row_to_rowMPI(element,temp,temp2,2,1,N-1,pos);	//  5/g2 <->  6/g4
      dssum_row_to_rowMPI(element,temp,temp2,3,1,N-1,pos);	//  9/g4 <-> 10/g2
      dssum_row_to_rowMPI(element,temp,temp2,0,1,N-1,pos);	//  1/g2 <->  2/g4
      dssum_row_to_rowMPI(element,temp,temp2,1,1,N-1,pos);	//  3/g4 <->  4/g2
#ifdef HAS_MPI
      MPI_Status res;
      MPI_Sendrecv(temp2,pos,MPI_PACKED,0,tag+2,
          temp3,pos,MPI_PACKED,0,tag+2,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssum_row_to_rowMPI(element,temp,temp3,2,1,N-1,pos,true);	//  5/g2 <-> 6/g4
      dssum_row_to_rowMPI(element,temp,temp3,3,1,N-1,pos,true);	//  9/g4 <-> 10/g2
      dssum_row_to_rowMPI(element,temp,temp3,0,1,N-1,pos,true);	//  1/g2 <->  2/g4
      dssum_row_to_rowMPI(element,temp,temp3,1,1,N-1,pos,true);	//  3/g4 <->  4/g2
      dssum_row_to_colb(element,3,4,3,2,  0,N-1   ); 				// 10/g4 >-< 11/g3
      dssum_col_to_col (element,4,5,0,2,  0,N-1   );				// 11/g1 <-> 12/g3
      dssum_col_to_row (element,5,2,0,3,  0,  0   );				// 12/g1 <-> 5/g4
      dssum_col_to_col (element,0,1,2,0,N-1,  0   );				//  1/g3 <->  4/g1

      /* lower left */
      dssum_col_to_col (element,0,2,0,2,  0,N-1   ); 	//  1/g1 <->  5/g3
      dssum_row_to_row (element,0,5,3,1,  0,N-1,-1); 	//  1/g4 <-> 12/g2
      element[m_range[0][m_range[0].size()-1]+5*m_size2][  0][N-1] =
        element[m_range[0][                  0]           ][  0][  0];
      /* top left */
      dssum_col_to_col (element,1,3,2,0,N-1,  0   );  //  4/g3 <-> 10/g1
      dssum_row_to_row (element,1,4,3,1,  0,N-1, 1);  //  4/g4 <-> 11/g2
      element[m_range[1][  0]+4*m_size2][N-1][N-1] =
        element[m_range[2][  0]+ 1*m_size2][N-1][  0];
    } else { // 2 3 6 7 8 9
      /* sum along the block boundaries */
      int pos=0;
      dssum_row_to_rowMPI(element,temp,temp2,2,3,0,pos);		//  5/g2 <->  6/g4
      dssum_row_to_rowMPI(element,temp,temp2,5,3,0,pos);		//  9/g4 <-> 10/g2
      dssum_row_to_rowMPI(element,temp,temp2,0,3,0,pos);		//  1/g2 <->  2/g4
      dssum_row_to_rowMPI(element,temp,temp2,1,3,0,pos);		//  3/g4 <->  4/g2
#ifdef HAS_MPI
      MPI_Status res;
      MPI_Sendrecv(temp2,pos,MPI_PACKED,1,tag+2,
          temp3,pos,MPI_PACKED,1,tag+2,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssum_row_to_rowMPI(element,temp,temp3,2,3,0,pos,true);		//  5/g2 <->  6/g4
      dssum_row_to_rowMPI(element,temp,temp3,5,3,0,pos,true);		//  9/g4 <-> 10/g2
      dssum_row_to_rowMPI(element,temp,temp3,0,3,0,pos,true);		//  1/g2 <->  2/g4
      dssum_row_to_rowMPI(element,temp,temp3,1,3,0,pos,true);		//  3/g4 <->  4/g2
      dssum_row_to_colb(element, 2, 3,1,0,N-1,  0); 	//  6/g2 >-<  7/g1
      dssum_col_to_col (element, 3, 4,2,0,N-1,  0); 	//  7/g3 <->  8/g1
      dssum_col_to_row (element, 4, 5,2,1,N-1,N-1); 	//  8/g3 <->  9/g2
      dssum_col_to_col (element, 0, 1,2,0,N-1,0);		//  2/g3 <->  3/g1

      /* lower right */
      dssum_col_to_col (element,0, 2,0,2,  0,N-1   ); //  2/g1 <->  6/g3
      dssum_row_to_row (element,0, 3,1,3,N-1,  0,-1); //  2/g2 <->  7/g4
      element[m_range[0][                  0]+3*m_size2][  0][  0] = 
        element[m_range[0][m_range[0].size()-1]][  0][N-1];

      /* top right */
      dssum_col_to_col (element,1, 5,2,0,N-1,  0   ); //  3/g3 <->  9/g1
      dssum_row_to_row (element,1, 4,1,3,N-1,  0, 1); //  3/g2 <->  8/g4
      element[m_range[2][0]+4*m_size2][N-1][  0] =
        element[m_range[1][0]+m_size2][N-1][N-1];
    }
  } else if( size == 4 ) {
    if( rank == 0 ) { // 1 5 12
      /* sum along the block boundaries */
      dssum_row_to_rowMPI(element,temp,temp2,0,1,N-1,pos);   	//  1/g2 <->  2/g4
      dssum_row_to_rowMPI(element,temp,temp2,1,1,N-1,pos); 	//  5/g2 <->  6/g4
#ifdef HAS_MPI
      MPI_Status res;
      MPI_Sendrecv(temp2,pos,MPI_PACKED,1,tag+1,
          temp3,pos,MPI_PACKED,1,tag+1,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssum_row_to_rowMPI(element,temp,temp3,0,1,N-1,pos,true);  	//  1/g2 <->  2/g4
      dssum_row_to_rowMPI(element,temp,temp3,1,1,N-1,pos,true); 	//  5/g2 <->  6/g4

      pos = 0;
      dssum_col_to_colMPI(element,temp,temp2,2,2,N-1,pos); 	// 11/g1 <-> 12/g3
      dssum_col_to_colMPI(element,temp,temp2,0,2,N-1,pos);   	//  1/g3 <->  4/g1
#ifdef HAS_MPI
      MPI_Sendrecv(temp2,pos,MPI_PACKED,3,tag+2,
          temp3,pos,MPI_PACKED,3,tag+2,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssum_col_to_colMPI(element,temp,temp3,2,2,N-1,pos,true); 	// 11/g1 <-> 12/g3
      dssum_col_to_colMPI(element,temp,temp3,0,2,N-1,pos,true);  	//  1/g3 <->  4/g1

      dssum_col_to_row (element,2,1,0,3,  0,  0   ); 			// 12/g1 <->  5/g4

      /* lower left */
      dssum_col_to_col (element,0,1,0,2,  0,N-1   ); //  1/g1 <->  5/g3
      dssum_row_to_row (element,0,2,3,1,  0,N-1,-1); //  1/g4 <-> 12/g2
      element[m_range[0][m_range[0].size()-1]+2*m_size2][  0][N-1] =
        element[m_range[0][                  0]           ][  0][  0];
    } else if( rank == 1 ) { // 2 6 7
      /* sum along the block boundaries */
      dssum_row_to_rowMPI(element,temp,temp2,0,3,0,pos);	//  1/g2 <->  2/g4
      dssum_row_to_rowMPI(element,temp,temp2,1,3,0,pos);	//  5/g2 <->  6/g4
#ifdef HAS_MPI
      MPI_Status res;
      MPI_Sendrecv(temp2,pos,MPI_PACKED,0,tag+1,
          temp3,pos,MPI_PACKED,0,tag+1,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssum_row_to_rowMPI(element,temp,temp3,0,3,0,pos,true);	//  1/g2 <->  2/g4
      dssum_row_to_rowMPI(element,temp,temp3,1,3,0,pos,true);	//  5/g2 <->  6/g4

      dssum_row_to_colb(element, 1, 2,1,0,N-1, 0);			//  6/g2 >-<  7/g1

      pos = 0;
      dssum_col_to_colMPI(element,temp,temp2,2,2,N-1,pos); 	//  7/g3 <->  8/g1
      dssum_col_to_colMPI(element,temp,temp2,0,2,N-1,pos);	//  2/g3 <->  3/g1
#ifdef HAS_MPI
      MPI_Sendrecv(temp2,pos,MPI_PACKED,2,tag+3,
          temp3,pos,MPI_PACKED,2,tag+3,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssum_col_to_colMPI(element,temp,temp3,2,2,N-1,pos,true); 	//  7/g3 <->  8/g1
      dssum_col_to_colMPI(element,temp,temp3,0,2,N-1,pos,true);	//  2/g3 <->  3/g1

      /* lower right */
      dssum_col_to_col (element,0, 1,0,2,  0,N-1   ); //  2/g1 <->  6/g3
      dssum_row_to_row (element,0, 2,1,3,N-1,  0,-1); //  2/g2 <->  7/g4
      element[m_range[0][                  0]+2*m_size2][  0][  0] = 
        element[m_range[0][m_range[0].size()-1]][  0][N-1];
    } else if( rank == 2 ) { // 3 8 9
      /* sum along the block boundaries */
      pos = 0;
      dssum_row_to_rowMPI(element,temp,temp2,2,3,0,pos);		//  9/g4 <-> 10/g2
      dssum_row_to_rowMPI(element,temp,temp2,0,3,0,pos);		//  3/g4 <->  4/g2
#ifdef HAS_MPI
      MPI_Status res;
      MPI_Sendrecv(temp2,pos,MPI_PACKED,3,tag+4,
          temp3,pos,MPI_PACKED,3,tag+4,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssum_row_to_rowMPI(element,temp,temp3,2,3,0,pos,true);		//  9/g4 <-> 10/g2
      dssum_row_to_rowMPI(element,temp,temp3,0,3,0,pos,true);		//  3/g4 <->  4/g2

      pos = 0;
      dssum_col_to_colMPI(element,temp,temp2,1,0,0,pos);	//  7/g3 <->  8/g1
      dssum_col_to_colMPI(element,temp,temp2,0,0,0,pos);	//  2/g3 <->  3/g1
#ifdef HAS_MPI
      MPI_Sendrecv(temp2,pos,MPI_PACKED,1,tag+3,
          temp3,pos,MPI_PACKED,1,tag+3,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssum_col_to_colMPI(element,temp,temp3,1,0,0,pos,true);	//  7/g3 <->  8/g1
      dssum_col_to_colMPI(element,temp,temp3,0,0,0,pos,true);	//  2/g3 <->  3/g1

      dssum_col_to_row   (element,1,2,2,1,N-1,N-1); 	//  8/g3 <->  9/g2

      /* top right */
      dssum_col_to_col (element,0, 2,2,0,N-1,  0   ); //  3/g3 <->  9/g1
      dssum_row_to_row (element,0, 1,1,3,N-1,  0, 1); //  3/g2 <->  8/g4
      element[m_range[2][0]+m_size2][N-1][  0] = element[m_range[1][0]][N-1][N-1];
    } else {
      /* sum along the block boundaries */
      pos = 0;
      dssum_row_to_rowMPI(element,temp,temp2,1,1,N-1,pos); 	//  9/g4 <-> 10/g2
      dssum_row_to_rowMPI(element,temp,temp2,0,1,N-1,pos); 	//  3/g4 <->  4/g2
#ifdef HAS_MPI
      MPI_Status res;
      MPI_Sendrecv(temp2,pos,MPI_PACKED,2,tag+4,
          temp3,pos,MPI_PACKED,2,tag+4,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssum_row_to_rowMPI(element,temp,temp3,1,1,N-1,pos,true); 	//  9/g4 <-> 10/g2
      dssum_row_to_rowMPI(element,temp,temp3,0,1,N-1,pos,true); 	//  3/g4 <->  4/g2

      pos = 0;
      dssum_col_to_colMPI(element,temp,temp2,2,0,0,pos);	// 11/g1 <-> 12/g3
      dssum_col_to_colMPI(element,temp,temp2,0,0,0,pos); 	//  1/g3 <->  4/g1
#ifdef HAS_MPI
      MPI_Sendrecv(temp2,pos,MPI_PACKED,0,tag+2,
          temp3,pos,MPI_PACKED,0,tag+2,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssum_col_to_colMPI(element,temp,temp3,2,0,0,pos,true);	// 11/g1 <-> 12/g3
      dssum_col_to_colMPI(element,temp,temp3,0,0,0,pos,true);	//  1/g3 <->  4/g1

      dssum_row_to_colb(element,1,2,3,2,  0,N-1   ); 	// 10/g4 <-> 11/g3

      /* top left */
      dssum_col_to_col (element,0,1,2,0,N-1,  0   );  //  4/g3 <-> 10/g1
      dssum_row_to_row (element,0,2,3,1,  0,N-1, 1);  //  4/g4 <-> 11/g2
      element[m_range[1][  0]+2*m_size2][N-1][N-1] =
        element[m_range[2][  0]][N-1][  0];
    }
  }
  if( size > 1 ) {
    delete[] temp;
    delete[] temp2;
    delete[] temp3;
  }
}

void bigCircleGeometry::mask(basics::matrixStack& op, int rank, int size) const
{
  int N = op[0].rows();
  for( int i=0;i<m_range[0].size();++i) {
    if( size == 1 ) {
      op[m_range[0][i]+4*m_size2][0].clear();    // 5/g1
      op[m_range[0][i]+5*m_size2][0].clear();    // 6/g1
      op[m_range[1][i]+6*m_size2].clearRow(N-1); // 7/g2
      op[m_range[1][i]+7*m_size2].clearRow(N-1); // 8/g2
      op[m_range[2][i]+8*m_size2][N-1].clear();  // 9/g3
      op[m_range[2][i]+9*m_size2][N-1].clear();  // 10/g3
      op[m_range[3][i]+10*m_size2].clearRow(0);  // 11/g4
      op[m_range[3][i]+11*m_size2].clearRow(0);  // 12/g4
    } else if( size == 2 ) {
      if( rank == 0 ) {
        op[m_range[0][i]+4*m_size2][0].clear();   // 5/g1
        op[m_range[0][i]+5*m_size2][0].clear();   // 6/g1
      } else {
        op[m_range[1][i]+0*m_size2].clearRow(N-1); // 7/g2
        op[m_range[1][i]+1*m_size2].clearRow(N-1); // 8/g2
        op[m_range[2][i]+2*m_size2][N-1].clear();  // 9/g3
        op[m_range[2][i]+3*m_size2][N-1].clear(); // 10/g3
        op[m_range[3][i]+4*m_size2].clearRow(0);  // 11/g4
        op[m_range[3][i]+5*m_size2].clearRow(0);  // 12/g4
      }
    } else if ( size == 3 ) {
      if( rank == 1 ) {
        op[m_range[0][i]+0*m_size2][0].clear();    // 5/g1
        op[m_range[0][i]+1*m_size2][0].clear();    // 6/g1
        op[m_range[1][i]+2*m_size2].clearRow(N-1); // 7/g2
        op[m_range[1][i]+3*m_size2].clearRow(N-1); // 8/g2
      } else if( rank == 2 ) {
        op[m_range[2][i]+0*m_size2][N-1].clear();  //  9/g3
        op[m_range[2][i]+1*m_size2][N-1].clear();  // 10/g3
        op[m_range[3][i]+2*m_size2].clearRow(0);   // 11/g4
        op[m_range[3][i]+3*m_size2].clearRow(0);   // 12/g4
      }
    } else if( size == 4 ) {
      if( rank == 1 ) {
        op[m_range[0][i]+1*m_size2][0].clear();    // 5/g1
        op[m_range[0][i]+2*m_size2][0].clear();    // 6/g1
      } else if( rank == 2 ) {
        op[m_range[1][i]+0*m_size2].clearRow(N-1); // 7/g2
        op[m_range[1][i]+1*m_size2].clearRow(N-1); // 8/g2
        op[m_range[2][i]+2*m_size2][N-1].clear();  // 9/g3
      } else if( rank == 3 ) {
        op[m_range[2][i]+0*m_size2][N-1].clear();  // 10/g3
        op[m_range[3][i]+1*m_size2].clearRow(0);   // 11/g4
        op[m_range[3][i]+2*m_size2].clearRow(0);   // 12/g4
      }
    } else if( size == 6 ) {
      if( rank == 2 ) {
        op[m_range[0][i]+0*m_size2][0].clear();    // 5/g1
        op[m_range[0][i]+1*m_size2][0].clear();    // 6/g1
      } else if( rank == 3 ) {
        op[m_range[1][i]+0*m_size2].clearRow(N-1); // 7/g2
        op[m_range[1][i]+1*m_size2].clearRow(N-1); // 8/g2
      } else if( rank == 4 ) {
        op[m_range[2][i]+0*m_size2][N-1].clear();  // 9/g3
        op[m_range[2][i]+1*m_size2][N-1].clear();  // 10/g3
      } else if( rank == 5 ) {
        op[m_range[3][i]+0*m_size2].clearRow(0);   // 11/g4
        op[m_range[3][i]+1*m_size2].clearRow(0);   // 12/g4
      }
    } else if( size == 12 ) {
      if( rank == 4 )
        op[m_range[0][i]+0*m_size2][0].clear();    // 5/g1
      else if( rank == 5 )
        op[m_range[0][i]+0*m_size2][0].clear();    // 6/g1
      else if( rank == 6 )
        op[m_range[1][i]+0*m_size2].clearRow(N-1); // 7/g2
      else if( rank == 7 )
        op[m_range[1][i]+0*m_size2].clearRow(N-1); // 8/g2
      else if( rank == 8 )
        op[m_range[2][i]+0*m_size2][N-1].clear();  // 9/g3
      else if( rank == 9 )
        op[m_range[2][i]+0*m_size2][N-1].clear();  // 10/g3
      else if( rank == 10 )
        op[m_range[3][i]+0*m_size2].clearRow(0);   // 11/g4
      else if( rank == 11 )
        op[m_range[3][i]+0*m_size2].clearRow(0);   // 12/g4
    }
  }
}

void bigCircleGeometry::mask(basics::matricesStack& op, int rank, int size) const
{
  op.at(0).clear();
  op.at(op[0].matrices()-1).clear();
  for( int l=1;l<op[0].matrices()-1;++l )
    mask(op.at(l),rank,size);
}

void bigCircleGeometry::addInnerGroup(basics::coarseGrid& result, int n, int offs) const
{
  for( int i=0;i<2;++i ) { // which quadrant of the group?
    for( int m=0;m<2;++m ) {
      basics::coarseDescriptor desc;
      desc.type  = basics::coarseDescriptor::FULL_OVERLAP;
      desc.type2 = basics::coarseDescriptor::NO_FLIP;
      desc.type3 = basics::coarseDescriptor::NORMAL;
      desc.size1 = m_size1/2;
      desc.size2 = m_size1/2;
      desc.size3 = result.size()+offs;
      for( int j=m_size1/2-1;j>=0;--j )
        for( int k=0;k<m_size1/2;++k )
          desc.elements.push_back(m_range[3][j+(1-m)*m_size1/2]+
              n*m_size2+k+i*m_size1/2);
      result.push_back(desc);
    }
  }
}

void bigCircleGeometry::addOuterGroup1(basics::coarseGrid& result, int n, int offs) const
{
  for( int i=0;i<2;++i ) { // which half of the group?
    basics::coarseDescriptor desc;
    desc.type = basics::coarseDescriptor::NO_BORDER_OVERLAP;
    desc.type2 = basics::coarseDescriptor::NO_FLIP;
    desc.type3 = basics::coarseDescriptor::NORMAL;
    desc.size1 = m_size1/2;
    desc.size2 = m_size1;
    desc.size3 = result.size()+offs;
    for( int j=m_size1-1;j>=0;--j )
      for( int k=0;k<m_size1/2;++k )
        desc.elements.push_back(m_range[3][j]+n*m_size2+i*m_size1/2+k);

    result.push_back(desc);
  }
}

void bigCircleGeometry::addOuterGroup2(basics::coarseGrid& result, int n, int offs) const
{
  for( int i=0;i<2;++i ) { // which half of the group?
    basics::coarseDescriptor desc;
    desc.type  = basics::coarseDescriptor::NO_BORDER_OVERLAP;
    desc.type2 = basics::coarseDescriptor::FLIP_XY;
    desc.type3 = basics::coarseDescriptor::BACKWARDS_ROW;
    desc.size1 = m_size1/2;
    desc.size2 = m_size1;
    desc.size3 = result.size()+offs;
    for( int j=m_size1-1;j>=0;--j )
      for( int k=0;k<m_size1/2;++k)
        desc.elements.push_back(m_range[0][j]+n*m_size2-k*m_size1-i*m_size2/2);

    result.push_back(desc);
  }
}

void bigCircleGeometry::addOuterGroup3(basics::coarseGrid& result, int n, int offs) const
{
  for( int i=0;i<2;++i ) { // which half of the group?
    basics::coarseDescriptor desc;
    desc.type = basics::coarseDescriptor::NO_BORDER_OVERLAP;
    desc.type2 = basics::coarseDescriptor::NO_FLIP;
    desc.type3 = basics::coarseDescriptor::BACKWARDS_COL_BACKWARDS;
    desc.size1 = m_size1/2;
    desc.size2 = m_size1;
    desc.size3 = result.size()+offs;
    for( int j=0;j<m_size1;++j )
      for( int k=0;k<m_size1/2;++k)
        desc.elements.push_back(m_range[1][j]+n*m_size2-i*m_size1/2-k);
    result.push_back(desc);
  }
}

void bigCircleGeometry::addOuterGroup4(basics::coarseGrid& result, int n, int offs) const
{
  for( int i=0;i<2;++i ) { // which half of the group?
    basics::coarseDescriptor desc;
    desc.type = basics::coarseDescriptor::NO_BORDER_OVERLAP;
    desc.type2 = basics::coarseDescriptor::FLIP_XY;
    desc.type3 = basics::coarseDescriptor::ROW_BACKWARDS;
    desc.size1 = m_size1/2;
    desc.size2 = m_size1;
    desc.size3 = result.size()+offs;
    for( int j=0;j<m_size1;++j ) 
      for( int k=0;k<m_size1/2;++k)
        desc.elements.push_back(m_range[2][j]+n*m_size2+i*m_size2/2+k*m_size1);
    result.push_back(desc);
  }
}

vector<basics::coarseGrid> bigCircleGeometry::getCoarseGroups(int size) const
{
  vector<basics::coarseGrid> result;
  assert( m_grid.size() >= 48 );

  if( size == 1 ) {
    result.resize(1);

    for( int n=0;n<4;++n )
      addInnerGroup(result[0],n);

    for( int n=4;n<6;++n )
      addOuterGroup1(result[0],n);
    for( int n=6;n<8;++n )
      addOuterGroup2(result[0],n);
    for( int n=8;n<10;++n )
      addOuterGroup3(result[0],n);
    for( int n=10;n<12;++n )
      addOuterGroup4(result[0],n);

  } else if( size == 2 ) {
    result.resize(2);
    int g = 0;
    /* rank = 0 */
    for( int n=0;n<4;++n )
      addInnerGroup(result[0],n);
    for( int n=4;n<6;++n )
      addOuterGroup1(result[0],n);

    /* rank = 1 */
    for( int n=0;n<2;++n )
      addOuterGroup2(result[1],n,result[0].size());
    for( int n=2;n<4;++n )
      addOuterGroup3(result[1],n,result[0].size());
    for( int n=4;n<6;++n )
      addOuterGroup4(result[1],n,result[0].size());
  } else if( size == 4 ) {
    result.resize(4);
    int curr=0;
    /* rank = 0 */
    for( int n=0;n<3;++n )
      addInnerGroup(result[0],n,curr);
    curr += 12;

    /* rank = 1 */
    addInnerGroup(result[1],0,curr);
    for( int n=1;n<3;++n )
      addOuterGroup1(result[1],n,curr);
    curr += 8;

    /* rank = 2 */
    for( int n=0;n<2;++n )
      addOuterGroup2(result[2],n,curr);
    addOuterGroup3(result[2],2,curr);
    curr += 6;

    /* rank = 3 */
    addOuterGroup3(result[3],0,curr);
    for( int n=1;n<3;++n )
      addOuterGroup4(result[3],n,curr);
  } else if( size == 12 ) {
    result.resize(12);

    for( int i=0;i<4;++i )
      addInnerGroup(result[i],0,4*i);
    for( int i=4;i<6;++i )
      addOuterGroup1(result[i],0,16+2*(i-4));
    for( int i=6;i<8;++i )
      addOuterGroup2(result[i],0,20+2*(i-6));
    for( int i=8;i<10;++i )
      addOuterGroup3(result[i],0,24+2*(i-8));
    for( int i=10;i<12;++i )
      addOuterGroup4(result[i],0,28+2*(i-10));
  }
  return( result );
}

void bigCircleGeometry::processFineToCoarse(basics::Matrix& foo,
    const basics::matrixStack& u,
    const basics::coarseDescriptor& desc,
    bool neumann) const
{
  int N = u[0].cols();
  int startcol=0;
  for( int k=0;k<desc.size2;++k ) {
    int startrow=0;
    int add = (k>0?0:1);
    if( neumann && desc.type == basics::coarseDescriptor::NO_BORDER_OVERLAP )
      add = 0;
    for( int i=0;i<desc.size1;++i ) {
      int element = desc.elements[i+k*desc.size1];
      for( int j=0;j<N;++j ) {
        if( desc.type == basics::coarseDescriptor::NO_BORDER_OVERLAP )
          if( k == 0 && j == N-1 )
            continue;
        if( desc.type3 == basics::coarseDescriptor::ROW_BACKWARDS ) {
          for( int i=0;i<N;++i )
            foo[j+startcol][startrow+i] = u[element][N-1-i][j+add];
        } else if( desc.type3 == basics::coarseDescriptor::BACKWARDS_COL_BACKWARDS ) {
          for( int i=0;i<N;++i )
            foo[j+startcol][i+startrow] = u[element][N-1-(j+add)][N-1-i];
        } else if( desc.type3 == basics::coarseDescriptor::BACKWARDS_ROW ) {
          for( int l=0;l<N;++l )
            foo[j+startcol][startrow+l] = u[element][l][N-1-(j+add)];
        } else if( desc.type == basics::coarseDescriptor::NO_BORDER_OVERLAP ) {
          for( int l=0;l<N;++l )
            foo[j+startcol][startrow+l] = u[element][j+add][l];
        } else {
          for( int l=0;l<N;++l )
            foo[j+startcol][startrow+l] = u[element][j][l];
        }
      }
      startrow += N-1;
    }
    startcol += N-1;
    if( !neumann && desc.type == basics::coarseDescriptor::NO_BORDER_OVERLAP && k == 0 )
      startcol--;
  }
}

void bigCircleGeometry::processFineToCoarseL2(basics::Matrix& foo,
    const basics::matrixStack& p, 
    const basics::coarseDescriptor& desc) const
{
  int N = p[0].cols();
  int startcol=0;
  for( int k=0;k<desc.size2;++k ) {
    int startrow=0;
    for( int i=0;i<desc.size1;++i ) {
      int element = desc.elements[i+k*desc.size1];
      for( int j=0;j<N;++j ) {
        if( desc.type3 ==  basics::coarseDescriptor::ROW_BACKWARDS ) {
          for( int i=0;i<N;++i )
            foo[j+startcol][startrow+i] = p[element][N-1-i][j];
        } else if( desc.type3 == basics::coarseDescriptor::BACKWARDS_COL_BACKWARDS ) {
          for( int i=0;i<N;++i )
            foo[j+startcol][i+startrow] = p[element][N-1-j][N-1-i];
        } else if( desc.type3 == basics::coarseDescriptor::BACKWARDS_ROW ) {
          int N2=N;
          for( int l=0;l<N;++l )
            foo[j+startcol][startrow+l] = p[element][l][N-1-j];
        } else if( desc.type == basics::coarseDescriptor::NO_BORDER_OVERLAP ) {
          for( int l=0;l<N;++l )
            foo[j+startcol][startrow+l] = p[element][j][l];
        } else {
          for( int l=0;l<N;++l )
            foo[j+startcol][startrow+l] = p[element][j][l];
        }
      }
      startrow += N;
    }
    startcol += N;
  }
}

basics::matrixStack* bigCircleGeometry::getCoarseBuffer(int rows, int cols,
    const basics::coarseGrid& desc,
    const std::string& name, bool L2,
    bool neumann) const
{
  basics::matrixStack* result = new basics::matrixStack(name);
  for( int n=0;n<desc.size();++n ) {
    int size1, size2;
    if( L2 ) {
      size1 = rows*desc[n].size1;
      size2 = cols*desc[n].size2;
    } else {
      size1 = rows*desc[n].size1-(desc[n].size1-1);
      if( desc[n].type == basics::coarseDescriptor::FULL_OVERLAP )
        size2 = rows*desc[n].size1-(desc[n].size1-1);
      else
        size2 = cols*desc[n].size2-(desc[n].size2-1)-(neumann?0:1);
    }

    basics::Matrix temp("coarse",size1,size2);
    result->add(temp,1);
  }

  return( result );
}

basics::matricesStack* bigCircleGeometry::getCoarseBuffer(int rows, int cols, int matrices,
    const basics::coarseGrid& desc,
    const std::string& name, bool L2,
    bool neumann) const
{
  basics::matricesStack* result = new basics::matricesStack(name);
  for( int n=0;n<desc.size();++n ) {
    int size1, size2;
    if( L2 ) {
      size1 = rows*desc[n].size1;
      size2 = cols*desc[n].size2;
    } else {
      size1 = rows*desc[n].size1-(desc[n].size1-1);
      if( desc[n].type == basics::coarseDescriptor::FULL_OVERLAP )
        size2 = rows*desc[n].size1-(desc[n].size1-1);
      else
        size2 = cols*desc[n].size2-(desc[n].size2-1)-(neumann?0:1);
    }

    basics::Matrices temp("coarse",size1,size2,matrices);
    result->add(temp,1);
  }
  /* setup the 'fake' matrix stack */
  result->m_stacks.clear();
  for( int l=0;l<matrices;++l ) {
    std::vector<basics::Matrix*> vec;
    for( int i=0;i<result->size();++i ) 
      vec.push_back(&((*result)[i][l]));
    basics::matrixStack mat(result->name()+" (fake matrix stack)");
    mat.add(vec);
    result->m_stacks.push_back(mat);
  }	

  return( result );
}

void bigCircleGeometry::fineToCoarse(basics::matrixStack& result,
    const basics::matrixStack& u,
    const basics::coarseGrid& desc,
    bool neumann) const
{
  int max=result.size();
  //#pragma omp parallel for schedule(static)
  for( int n=0;n<max;++n )
    processFineToCoarse(result[n],u,desc[n],neumann);
}

void bigCircleGeometry::fineToCoarse(basics::matricesStack& result,
    const basics::matricesStack& u,
    const basics::coarseGrid& desc,
    bool neumann) const
{
  int max=result[0].matrices();
  //#pragma omp parallel for schedule(static)
  for( int l=0;l<max;++l )
    fineToCoarse(result.at(l),u.at(l),desc,neumann);
}

void bigCircleGeometry::fineToCoarseL2(basics::matrixStack& result,
    const basics::matrixStack& p,
    const basics::coarseGrid& desc) const
{
  int max=result.size();
  //#pragma omp parallel for schedule(static)
  for( int n=0;n<max;++n )
    processFineToCoarseL2(result[n],p,desc[n]);
}

void bigCircleGeometry::fineToCoarseL2(basics::matricesStack& result,
    const basics::matricesStack& p,
    const basics::coarseGrid& desc) const
{
  int max=result[0].matrices();
  //#pragma omp parallel for schedule(static)
  for( int l=0;l<max;++l )
    for( int n=0;n<result.size();++n )
      processFineToCoarseL2(result.at(l)[n],p.at(l),desc[n]);
}

void bigCircleGeometry::processFineToCoarseRestrictionInner(basics::matrixStack& buffer,
    int n) const
{
  int skew = n*m_size2;
  int N = buffer[0].rows();
  for( int j=1;j<=m_size1/2-1;++j ) {
    for( int i=0;i<m_range[2].size();++i) {
      int elem11 = m_range[2][i]+skew+j*m_size1;
      int elem12 = m_range[2][i]+skew+(j-1)*m_size1;
      buffer[elem11][N-1] += buffer[elem12][0];
      buffer[elem12][0]    = buffer[elem11][N-1];
      elem11 += m_size2/2;
      elem12 += m_size2/2;
      buffer[elem11][N-1] += buffer[elem12][0];
      buffer[elem12][0]    = buffer[elem11][N-1];
    }
  }
  for( int j=1;j<=m_size1/2-1;++j ) {
    for( int i=0;i<m_range[3].size();++i) {
      int elem21 = m_range[3][i]+skew+j;
      int elem22 = m_range[3][i]+skew+(j-1);
      for( int l=0;l<N;++l ) {
        buffer[elem21][l][0] += buffer[elem22][l][N-1];
        buffer[elem22][l][N-1] = buffer[elem21][l][0];
      }
      elem21 += m_size1/2;
      elem22 += m_size1/2;
      for( int l=0;l<N;++l ) {
        buffer[elem21][l][0] += buffer[elem22][l][N-1];
        buffer[elem22][l][N-1] = buffer[elem21][l][0];
      }
    }
  }
}

void bigCircleGeometry::processFineToCoarseRestrictionOuter1(basics::matrixStack& buffer,
    int n) const
{
  int skew = n*m_size2;
  int N = buffer[0].rows();
  for( int j=1;j<=m_size1-1;++j ) {
    for( int i=0;i<m_range[2].size();++i) {
      int elem11 = m_range[2][i]+skew+j*m_size1;
      int elem12 = m_range[2][i]+skew+(j-1)*m_size1;
      buffer[elem11][N-1] += buffer[elem12][0];
      buffer[elem12][0]    = buffer[elem11][N-1];
    }
  }
  for( int j=1;j<=m_size1/2-1;++j ) {
    for( int i=0;i<m_range[3].size();++i) {
      int elem21 = m_range[3][i]+skew+j;
      int elem22 = m_range[3][i]+skew+(j-1);
      for( int l=0;l<N;++l ) {
        buffer[elem21][l][0] += buffer[elem22][l][N-1];
        buffer[elem22][l][N-1] = buffer[elem21][l][0];
      }
      elem21 += m_size1/2;
      elem22 += m_size1/2;
      for( int l=0;l<N;++l ) {
        buffer[elem21][l][0] += buffer[elem22][l][N-1];
        buffer[elem22][l][N-1] = buffer[elem21][l][0];
      }
    }
  }
}

void bigCircleGeometry::processFineToCoarseRestrictionOuter2(basics::matrixStack& buffer,
    int n) const
{
  int skew = n*m_size2;
  int N = buffer[0].rows();
  for( int j=1;j<=m_size1-1;++j ) {
    for( int i=0;i<m_range[3].size();++i) {
      int elem21 = m_range[3][i]+skew+j;
      int elem22 = m_range[3][i]+skew+(j-1);
      for( int l=0;l<N;++l ) {
        buffer[elem21][l][0] += buffer[elem22][l][N-1];
        buffer[elem22][l][N-1] = buffer[elem21][l][0];
      }
    }
  }
  for( int j=1;j<=m_size1/2-1;++j ) {
    for( int i=0;i<m_range[2].size();++i) {
      int elem11 = m_range[2][i]+skew+j*m_size1;
      int elem12 = m_range[2][i]+skew+(j-1)*m_size1;
      buffer[elem11][N-1] += buffer[elem12][0];
      buffer[elem12][0]    = buffer[elem11][N-1];
      elem11 += m_size2/2;
      elem12 += m_size2/2;
      buffer[elem11][N-1] += buffer[elem12][0];
      buffer[elem12][0]    = buffer[elem11][N-1];
    }
  }
}

void bigCircleGeometry::fineToCoarseRestriction(basics::matrixStack& result,
    basics::matrixStack& buffer,
    const basics::matrixStack& u,
    const basics::coarseGrid& desc,
    bool neumann, int rank, int size) const
{
  buffer = u;
  if( size == 1 ) {
    for( int n=0;n<4;++n )
      processFineToCoarseRestrictionInner(buffer,n);
    for( int n=4;n<6;++n ) 
      processFineToCoarseRestrictionOuter1(buffer,n);
    for( int n=6;n<8;++n )
      processFineToCoarseRestrictionOuter2(buffer,n);
    for( int n=8;n<10;++n )
      processFineToCoarseRestrictionOuter1(buffer,n);
    for( int n=10;n<12;++n )
      processFineToCoarseRestrictionOuter2(buffer,n);
  } else if( size == 2 ) {
    if( rank == 0 )  {
      for( int n=0;n<4;++n )
        processFineToCoarseRestrictionInner(buffer,n);
      for( int n=4;n<6;++n ) 
        processFineToCoarseRestrictionOuter1(buffer,n);
    } else {
      for( int n=0;n<2;++n )
        processFineToCoarseRestrictionOuter2(buffer,n);
      for( int n=2;n<4;++n )
        processFineToCoarseRestrictionOuter1(buffer,n);
      for( int n=4;n<6;++n )
        processFineToCoarseRestrictionOuter2(buffer,n);
    }
  } else if( size == 4 ) {
    if( rank == 0 ) {
      for( int n=0;n<3;++n )
        processFineToCoarseRestrictionInner(buffer,n);
    } else if( rank == 1 ) {
      processFineToCoarseRestrictionInner(buffer,0);
      for( int n=1;n<3;++n ) 
        processFineToCoarseRestrictionOuter1(buffer,n);
    } else if( rank == 2 ) {
      for( int n=0;n<2;++n )
        processFineToCoarseRestrictionOuter2(buffer,n);
      processFineToCoarseRestrictionOuter1(buffer,2);
    } else {
      processFineToCoarseRestrictionOuter1(buffer,0);
      for( int n=1;n<3;++n )
        processFineToCoarseRestrictionOuter2(buffer,n);
    }
  } else if( size == 12 ) {
    if( rank >= 0 && rank <= 3)
      processFineToCoarseRestrictionInner(buffer,0);
    if( rank >= 4 && rank < 6 ||
        rank >= 8 && rank < 10)
      processFineToCoarseRestrictionOuter1(buffer,0);
    if( rank >= 6 && rank < 8 ||
        rank >= 10 && rank < 12)
      processFineToCoarseRestrictionOuter2(buffer,0);
  }

  fineToCoarse(result,buffer,desc,neumann);
}

void bigCircleGeometry::fineToCoarseRestriction(basics::matricesStack& result,
    basics::matricesStack& buffer,
    const basics::matricesStack& u,
    const basics::coarseGrid& desc,
    bool neumann,
    int rank, int size) const
{
  int max=result[0].matrices();
  //#pragma omp parallel for schedule(static)
  for( int l=0;l<max;++l ) {
    fineToCoarseRestriction(result.at(l),buffer.at(l),u.at(l),desc,neumann,rank,size);
  }
}

void bigCircleGeometry::processCoarseToFine(basics::matrixStack& u,
    const basics::Matrix& foo,
    const basics::coarseDescriptor& desc,
    bool neumann) const
{
  int N = u[0].cols();
  int startcol=0;
  for( int k=0;k<desc.size2;++k ) {
    int add = (k>0?0:1);
    if( neumann && desc.type == basics::coarseDescriptor::NO_BORDER_OVERLAP )
      add = 0;
    int startrow=0;
    for( int i=0;i<desc.size1;++i ) {
      int element = desc.elements[i+k*desc.size1];
      for( int j=0;j<N;++j ) {
        if( desc.type == basics::coarseDescriptor::NO_BORDER_OVERLAP )
          if( k == 0 && j == N-1 && !neumann)
            continue;
        if( desc.type3 == basics::coarseDescriptor::ROW_BACKWARDS ) {
          for( int n=0;n<N;++n )
            u[element][N-1-n][j+add] = foo[j+startcol][startrow+n];
        } else if( desc.type3 == basics::coarseDescriptor::BACKWARDS_COL_BACKWARDS ) {
          for( int n=0;n<N;++n )
            u[element][N-1-(j+add)][N-1-n] = foo[j+startcol][n+startrow];
        } else if( desc.type3 == basics::coarseDescriptor::BACKWARDS_ROW ) {
          for( int l=0;l<N;++l )
            u[element][l][N-1-(j+add)] = foo[j+startcol][startrow+l];
        } else if( desc.type == basics::coarseDescriptor::NO_BORDER_OVERLAP ) {
          for( int l=0;l<N;++l )
            u[element][j+add][l] = foo[j+startcol][startrow+l];
        } else {
          for( int l=0;l<N;++l )
            u[element][j][l] = foo[j+startcol][startrow+l];
        }
      }
      startrow += N-1;
    }
    startcol += N-1;
    if( desc.type == basics::coarseDescriptor::NO_BORDER_OVERLAP && k == 0 && !neumann )
      startcol--;

  }
}

void bigCircleGeometry::processCoarseToFineL2(basics::matrixStack& p,
    const basics::Matrix& foo,
    const basics::coarseDescriptor& desc) const
{
  int N = p[0].cols();
  int startcol=0;
  for( int k=0;k<desc.size2;++k ) {
    int startrow=0;
    for( int i=0;i<desc.size1;++i ) {
      int element = desc.elements[i+k*desc.size1];
      for( int j=0;j<N;++j ) {
        if( desc.type3 == basics::coarseDescriptor::ROW_BACKWARDS ) {
          for( int i=0;i<N;++i )
            p[element][N-1-i][j] = foo[j+startcol][startrow+i];
        } else if( desc.type3 == basics::coarseDescriptor::BACKWARDS_COL_BACKWARDS ) {
          for( int i=0;i<N;++i )
            p[element][N-1-j][N-1-i] = foo[j+startcol][i+startrow];
        } else if( desc.type3 == basics::coarseDescriptor::BACKWARDS_ROW ) {
          for( int l=0;l<N;++l )
            p[element][l][N-1-j] = foo[j+startcol][startrow+l];
          //                    COPY_ROW_FROM_COL(p[element].data()[0]+N-1-j,foo.data()[j+startcol]+startrow);
        } else if( desc.type == basics::coarseDescriptor::NO_BORDER_OVERLAP ) {
          for( int l=0;l<N;++l )
            p[element][j][l] = foo[j+startcol][startrow+l];
          //                    COPY_COL_FROM_COL(p[element].data()[j],foo.data()[j+startcol]+startrow);
        } else {
          for( int l=0;l<N;++l )
            p[element][j][l] = foo[j+startcol][startrow+l];
          //                    COPY_COL_FROM_COL(p[element].data()[j],foo.data()[j+startcol]+startrow);
        }
      }
      startrow += N;
    }
    startcol += N;
  }
}

void bigCircleGeometry::coarseToFine(basics::matrixStack& u, 
    const basics::matrixStack& coarse,
    const basics::coarseGrid& desc,
    bool neumann) const
{
  int max=coarse.size();
  //#pragma omp parallel for schedule(static)
  for( int n=0;n<max;++n )
    processCoarseToFine(u,coarse[n],desc[n],neumann);
}

void bigCircleGeometry::coarseToFine(basics::matricesStack& u, 
    const basics::matricesStack& coarse,
    const basics::coarseGrid& desc,
    bool neumann) const
{
  int max=u[0].matrices();
  //#pragma omp parallel for schedule(static)
  for( int l=0;l<max;++l )
    coarseToFine(u.at(l),coarse.at(l),desc,neumann);
}

void bigCircleGeometry::coarseToFineL2(basics::matrixStack& p,
    const basics::matrixStack& coarse,
    const basics::coarseGrid& desc) const
{
  int max=desc.size();
  //#pragma omp parallel for schedule(static)
  for( int n=0;n<max;++n )
    processCoarseToFineL2(p,coarse[n],desc[n]);	
}

void bigCircleGeometry::coarseToFineL2(basics::matricesStack& p,
    const basics::matricesStack& coarse,
    const basics::coarseGrid& desc) const
{
  int max=p[0].matrices();
  //#pragma omp parallel for schedule(static)
  for( int l=0;l<max;++l )
    for( int n=0;n<coarse.size();++n )
      processCoarseToFineL2(p.at(l),coarse.at(l)[n],desc[n]);	
}

void bigCircleGeometry::dssumRowToRowMPI(basics::Matrix& result, char* temp, Real* temp2,
    int row, int& pos, bool sum) const
{
#ifdef HAS_MPI
  if( sum ) {
    MPI_Unpack(temp ,result.cols()*20*sizeof(Real),&pos,
        temp2,result.cols(),MPI_DOUBLE,MPI_COMM_WORLD);
    int N  = result.rows();
    int N2 = result.cols();
    AXPY_ROW_FROM_COL(result.data()[0]+row,temp2);
  } else  {
    MPI_Pack(result.row(row).data(),result.cols(),MPI_DOUBLE,temp,
        20*result.cols()*sizeof(Real),&pos,MPI_COMM_WORLD);
  }
#endif
}

void bigCircleGeometry::dssumColToColMPI(basics::Matrix& result, char* temp, Real* temp2,
    int col, int& pos, bool dosum) const
{
#ifdef HAS_MPI
  if( dosum ) {
    MPI_Unpack( temp,result.cols()*20*sizeof(Real),&pos,
        temp2,result.rows(),MPI_DOUBLE,MPI_COMM_WORLD);
    int N = result.rows();
    AXPY_COL_FROM_COL(result.data()[col],temp2)
  } else {
    MPI_Pack(result.data()[col],result.rows(),MPI_DOUBLE,temp,
        20*result.cols()*sizeof(Real),&pos,MPI_COMM_WORLD);
  }
#endif
}

void bigCircleGeometry::dssumCoarse(basics::matrixStack& result,
    const basics::coarseGrid& groups,
    int rank, int size, int mytag) const
{
  int N,N2;
#define DSSUM_COL_TO_COL(i1,i2,c1,c2) \
  result[i1][c1] += result[i2][c2]; \
  result[i2][c2]  = result[i1][c1];

#define DSSUM_COL_TO_COL_BACKWARDS(i1,i2,c1,c2) \
  for( int i=0;i<result[i1].rows();++i ) { \
    result[i1][c1][i] += result[i2][c2][result[i1].rows()-1-i]; \
    result[i2][c2][result[i1].rows()-1-i] = result[i1][c1][i];\
  }

#define DSSUM_COL_TO_COL_BACKWARDS_NO_EDGES(i1,i2,c1,c2,loffs) \
  for( int i=loffs;i<result[i1].rows()-(1-loffs);++i ) { \
    result[i1][c1][i] += result[i2][c2][result[i1].rows()-1-i]; \
    result[i2][c2][result[i1].rows()-1-i] = result[i1][c1][i]; \
  }

#define DSSUM_COL_TO_COL_NO_EDGES(i1,i2,c1,c2,edge) \
  AXPY_NO_EDGES(result[i1].data()[c1], \
      result[i2].data()[c2],edge) \
  COPY_NO_EDGES(result[i2].data()[c2], \
      result[i1].data()[c1],edge)

#define DSSUM_ROW_TO_ROW(i1,i2,r1,r2) \
  N  = result[i1].cols(); \
  N2 = result[i1].rows(); \
  BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealOne),BLASPREAL(result[i2].data()[0])+r2,BLASINT(N2),BLASPREAL(result[i1].data()[0])+r1,BLASINT(N2)); \
  BLASRPFX(copy,BLASINT(N),BLASPREAL(result[i1].data()[0])+r1,BLASINT(N2),BLASPREAL(result[i2].data()[0])+r2,BLASINT(N2));

#define DSSUM_COL_TO_ROW(i1,i2,r1,c2) \
  AXPY_ROW_FROM_COL(result[i1].data()[0]+r1,result[i2].data()[c2]); \
  COPY_COL_FROM_ROW(result[i2].data()[c2],result[i1].data()[0]+r1);

#define DSSUM_COL_TO_ROW_BACKWARDS(i1,i2,r1,c2) \
  for( int i=0;i<result[i1].cols();++i ) { \
    result[i1][i][r1] += result[i2][c2][result[i1].cols()-1-i]; \
    result[i2][c2][result[i1].cols()-1-i] = result[i1][i][r1]; \
  }
#define DSSUM_COL_TO_ROW_NO_EDGES(i1,i2,r1,c2,edge) \
  AXPY_ROW_FROM_COL(result[i1].data()[edge]+r1,result[i2].data()[c2]+edge); \
  COPY_COL_FROM_ROW(result[i2].data()[c2]+edge,result[i1].data()[edge]+r1);

#define DSSUM_COL_TO_ROW_BACKWARDS_NO_EDGES(i1,i2,r1,c2,edge) \
  for( int i=edge;i<result[i1].cols()-(1-edge);++i ) { \
    result[i1][i][r1] += result[i2][c2][result[i1].cols()-1-i]; \
    result[i2][c2][result[i1].cols()-1-i] = result[i1][i][r1]; \
  }

  N = result[0].cols();
  int tag = mytag*TAG_INTERVAL;
  tag += 13;
  int pos = 0;
  char* row;
  char* row3;
  Real* row2;
  if( size > 1 ) {
    row  = new char[result[0].cols()*20*sizeof(Real)];
    row3 = new char[result[0].cols()*20*sizeof(Real)];
    row2 = new Real[result[0].cols()*8];
  }
  if( size == 1 ) {
    /* 16 inner */
    DSSUM_COL_TO_COL( 0, 1,N-1,0)
      DSSUM_COL_TO_COL( 2, 3,N-1,0)
      DSSUM_COL_TO_COL( 4, 5,N-1,0)
      DSSUM_COL_TO_COL( 6, 7,N-1,0)
      DSSUM_COL_TO_COL( 1,12,N-1,0)
      DSSUM_COL_TO_COL( 3,14,N-1,0)
      DSSUM_COL_TO_COL( 5, 8,N-1,0)
      DSSUM_COL_TO_COL( 7,10,N-1,0)
      DSSUM_COL_TO_COL(12,13,N-1,0)
      DSSUM_COL_TO_COL(14,15,N-1,0)
      DSSUM_COL_TO_COL( 8, 9,N-1,0)
      DSSUM_COL_TO_COL(10,11,N-1,0)

      N = result[0].rows();
    DSSUM_ROW_TO_ROW( 0, 2,N-1,0)
      DSSUM_ROW_TO_ROW( 2, 4,N-1,0)
      DSSUM_ROW_TO_ROW( 4, 6,N-1,0)
      DSSUM_ROW_TO_ROW( 1, 3,N-1,0)
      DSSUM_ROW_TO_ROW( 3, 5,N-1,0)
      DSSUM_ROW_TO_ROW( 5, 7,N-1,0)
      DSSUM_ROW_TO_ROW(12,14,N-1,0)
      DSSUM_ROW_TO_ROW(14, 8,N-1,0)
      DSSUM_ROW_TO_ROW( 8,10,N-1,0)
      DSSUM_ROW_TO_ROW(13,15,N-1,0)
      DSSUM_ROW_TO_ROW(15, 9,N-1,0)
      DSSUM_ROW_TO_ROW( 9,11,N-1,0)

      int R1 = result[0].rows()-1;
    for( int i=16;i<31;++i ) {
      DSSUM_ROW_TO_ROW(i,i+1,R1,0)
    }
    DSSUM_ROW_TO_ROW(31,16,R1,0)

      N2 = result[0].rows()-1;
    N = result[16].cols();
    DSSUM_COL_TO_COL_NO_EDGES(0,16,0,N-1,1)
      N2 = result[2].rows();
    DSSUM_COL_TO_COL(2,17,0,N-1)
      DSSUM_COL_TO_COL(4,18,0,N-1)
      N2 = result[0].rows()-1;
    DSSUM_COL_TO_COL_NO_EDGES(6,19,0,N-1,0)

      N  = result[ 6].rows();
    N2 = result[20].rows();
    DSSUM_COL_TO_ROW( 6,20,result[ 6].rows()-1,result[20].cols()-1)
      DSSUM_COL_TO_ROW( 7,21,result[ 7].rows()-1,result[21].cols()-1)
      DSSUM_COL_TO_ROW(10,22,result[10].rows()-1,result[22].cols()-1)
      DSSUM_COL_TO_ROW(11,23,result[11].rows()-1,result[23].cols()-1)

      N2 = result[11].cols();
    N = result[24].cols();
    DSSUM_COL_TO_COL_BACKWARDS_NO_EDGES(11,24,N2-1,N-1,0)
      DSSUM_COL_TO_COL_BACKWARDS( 9,25,N2-1,N-1)
      DSSUM_COL_TO_COL_BACKWARDS(15,26,N2-1,N-1)
      DSSUM_COL_TO_COL_BACKWARDS_NO_EDGES(13,27,N2-1,N-1,1)

      N2 = result[28].cols();
    DSSUM_COL_TO_ROW_BACKWARDS(13,28,0,N2-1)
      DSSUM_COL_TO_ROW_BACKWARDS(12,29,0,N2-1)
      DSSUM_COL_TO_ROW_BACKWARDS( 1,30,0,N2-1)
      DSSUM_COL_TO_ROW_BACKWARDS( 0,31,0,N2-1)

      N  = result[16].cols();
    N2 = result[0].rows();
    result[16][N-1][   0] = result[31][N-1][N2-1] = result[ 0][0   ][   0];
    result[20][N-1][   0] = result[19][N-1][N2-1] = result[ 6][0   ][N2-1];
    result[23][N-1][N2-1] = result[24][N-1][   0] = result[11][N2-1][N2-1];
    result[27][N-1][N2-1] = result[28][N-1][   0] = result[13][N2-1][   0];
  } else if( size == 2 ) {
    if( rank == 0 ) {
      DSSUM_COL_TO_COL( 0, 1,N-1,0) 	//  0 <->  1
        DSSUM_COL_TO_COL( 2, 3,N-1,0) 	//  2 <->  3
        DSSUM_COL_TO_COL( 1, 4,N-1,0) 	//  1 <-> 12
        DSSUM_COL_TO_COL( 3, 6,N-1,0) 	//  3 <-> 14
        DSSUM_COL_TO_COL( 4, 5,N-1,0) 	// 12 <-> 13
        DSSUM_COL_TO_COL( 6, 7,N-1,0) 	// 14 <-> 15

        N = result[0].rows();
      DSSUM_ROW_TO_ROW( 0, 2,N-1,0) 	//  0 <->  2
        DSSUM_ROW_TO_ROW( 1, 3,N-1,0) 	//  1 <->  3
        DSSUM_ROW_TO_ROW( 4, 6,N-1,0) 	// 12 <-> 14
        DSSUM_ROW_TO_ROW( 5, 7,N-1,0) 	// 13 <-> 15

        pos=0;
      dssumRowToRowMPI(result[ 2],row,row2,N-1,pos);	//  2 <->  4
      dssumRowToRowMPI(result[ 3],row,row2,N-1,pos); 	//  3 <->  5
      dssumRowToRowMPI(result[ 6],row,row2,N-1,pos);	// 14 <->  8
      dssumRowToRowMPI(result[ 7],row,row2,N-1,pos);	// 15 <->  9
      N = result[9].rows();
      dssumRowToRowMPI(result[ 9],row,row2,N-1,pos);	// 17 <-> 18
      dssumRowToRowMPI(result[10],row,row2,  0,pos);	// 25 <-> 26 
#ifdef HAS_MPI
      MPI_Status res;
      MPI_Sendrecv( row,pos,MPI_PACKED,1,tag+2,
          row3,pos,MPI_PACKED,1,tag+2,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      N = result[0].rows();
      dssumRowToRowMPI(result[ 2],row3,row2,N-1,pos,true);	//  2 <->  4
      dssumRowToRowMPI(result[ 3],row3,row2,N-1,pos,true); 	//  3 <->  5
      dssumRowToRowMPI(result[ 6],row3,row2,N-1,pos,true);	// 14 <->  8
      dssumRowToRowMPI(result[ 7],row3,row2,N-1,pos,true);	// 15 <->  9
      N = result[9].rows();
      dssumRowToRowMPI(result[ 9],row3,row2,N-1,pos,true);	// 17 <-> 18
      dssumRowToRowMPI(result[10],row3,row2,  0,pos,true);	// 25 <-> 26 

      int R1 = result[0].rows()-1;
      DSSUM_ROW_TO_ROW(8,9,R1,0)
        N = result[9].rows();
      for( int i=10;i<15;++i ) {
        DSSUM_ROW_TO_ROW(i,i+1,R1,0)
      }
      DSSUM_ROW_TO_ROW(15,8,R1,0)

        N2 = result[0].rows()-1;
      N = result[8].cols();
      DSSUM_COL_TO_COL_NO_EDGES(0,8,0,N-1,1) 	// 0 <-> 16
        N2 = result[2].rows();
      DSSUM_COL_TO_COL(2,9,0,N-1)				// 2 <-> 17

        N2 = result[7].cols();
      N = result[10].cols();

      DSSUM_COL_TO_COL_BACKWARDS(7,10,N2-1,N-1)				// 15 <-> 26
        DSSUM_COL_TO_COL_BACKWARDS_NO_EDGES(5,11,N2-1,N-1,1) 	// 13 <-> 27

        N2 = result[12].cols();
      DSSUM_COL_TO_ROW_BACKWARDS(5,12,0,N2-1)				// 13 <-> 28
        DSSUM_COL_TO_ROW_BACKWARDS(4,13,0,N2-1)				// 12 <-> 29
        DSSUM_COL_TO_ROW_BACKWARDS(1,14,0,N2-1)				// 30 <->  1
        DSSUM_COL_TO_ROW_BACKWARDS(0,15,0,N2-1)				// 31 <->  0

        N  = result[8].cols();
      N2 = result[0].rows();
      result[ 8][N-1][   0] = result[15][N-1][N2-1] = result[0][0   ][   0]; // lower left
      result[11][N-1][N2-1] = result[12][N-1][   0] = result[5][N2-1][   0]; // upper left
    } else {
      DSSUM_COL_TO_COL( 0, 1,N-1,0) 			//  4 <->  5
        DSSUM_COL_TO_COL( 2, 3,N-1,0) 			//  6 <->  7
        DSSUM_COL_TO_COL( 1, 4,N-1,0) 			//  5 <->  8
        DSSUM_COL_TO_COL( 3, 6,N-1,0) 			//  7 <-> 10
        DSSUM_COL_TO_COL( 4, 5,N-1,0) 			//  8 <->  9
        DSSUM_COL_TO_COL( 6, 7,N-1,0) 			// 10 <-> 11
        N = result[0].rows();
      DSSUM_ROW_TO_ROW( 0, 2,N-1,0) 			//  4 <->  6
        DSSUM_ROW_TO_ROW( 1, 3,N-1,0) 			//  5 <->  7
        DSSUM_ROW_TO_ROW( 4, 6,N-1,0) 			//  8 <-> 10
        DSSUM_ROW_TO_ROW( 5, 7,N-1,0) 			//  9 <-> 11
        pos=0;
      dssumRowToRowMPI(result[ 0],row,row2,0,pos);	//  2 <->  4
      dssumRowToRowMPI(result[ 1],row,row2,0,pos);	//  3 <->  5
      dssumRowToRowMPI(result[ 4],row,row2,0,pos);	// 14 <->  8
      dssumRowToRowMPI(result[ 5],row,row2,0,pos);	// 15 <->  9
      dssumRowToRowMPI(result[ 8],row,row2,0,pos);	// 17 <-> 18
      N = result[15].rows();
      dssumRowToRowMPI(result[15],row,row2,N-1,pos);	// 25 <-> 26 
#ifdef HAS_MPI
      MPI_Status res;
      MPI_Sendrecv( row,pos,MPI_PACKED,0,tag+2,
          row3,pos,MPI_PACKED,0,tag+2,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      N = result[0].rows();
      dssumRowToRowMPI(result[ 0],row3,row2,0,pos,true);	//  2 <->  4
      dssumRowToRowMPI(result[ 1],row3,row2,0,pos,true);	//  3 <->  5
      dssumRowToRowMPI(result[ 4],row3,row2,0,pos,true);	// 14 <->  8
      dssumRowToRowMPI(result[ 5],row3,row2,0,pos,true);	// 15 <->  9
      dssumRowToRowMPI(result[ 8],row3,row2,0,pos,true);	// 17 <-> 18
      N = result[15].rows();
      dssumRowToRowMPI(result[15],row3,row2,N-1,pos,true);	// 25 <-> 26 
      N = result[0].rows();
      int R1 = result[0].rows()-1;
      for( int i=8;i<15;++i ) {
        DSSUM_ROW_TO_ROW(i,i+1,R1,0)
      }
      N = result[8].cols();
      N2 = result[2].rows();

      DSSUM_COL_TO_COL(0,8,0,N-1)					//  4 <-> 18
        N2 = result[0].rows()-1;
      DSSUM_COL_TO_COL_NO_EDGES(2,9,0,N-1,0)		//  6 <-> 19

        N  = result[ 2].rows();
      N2 = result[10].rows();
      int N3 = result[10].cols();
      DSSUM_COL_TO_ROW( 2,10,N-1,N3-1)		//  6 <-> 20
        DSSUM_COL_TO_ROW( 3,11,N-1,N3-1)		//  7 <-> 21
        DSSUM_COL_TO_ROW( 6,12,N-1,N3-1)		// 10 <-> 22
        DSSUM_COL_TO_ROW( 7,13,N-1,N3-1)		// 11 <-> 23
        N2 = result[7].cols();
      N = result[14].cols();
      DSSUM_COL_TO_COL_BACKWARDS_NO_EDGES(7,14,N2-1,N-1,0)	// 11 <-> 24
        DSSUM_COL_TO_COL_BACKWARDS( 5,15,N2-1,N-1)			 	//  9 <-> 25

        N  = result[8].cols();
      N2 = result[0].rows();
      result[10][N-1][   0] = result[ 9][N-1][N2-1] = result[2][0   ][N2-1]; // lower right
      result[13][N-1][N2-1] = result[14][N-1][   0] = result[7][N2-1][N2-1]; // upper right
    }
  } else if( size == 4 ) {
    if( rank == 0 ) {
      DSSUM_COL_TO_COL( 0, 1,N-1,0)		//  0 <->  1
        DSSUM_COL_TO_COL( 2, 3,N-1,0) 		//  2 <->  3
        N = result[0].rows();
      DSSUM_ROW_TO_ROW( 0, 2,N-1,0) 		//  0 <->  2
        DSSUM_ROW_TO_ROW( 1, 3,N-1,0) 		//  1 <->  3
        int R1 = result[0].rows()-1;
      DSSUM_ROW_TO_ROW(4,5,R1,0)			// 16 <-> 17
        DSSUM_ROW_TO_ROW(6,7,R1,0 )			// 30 <-> 31
        DSSUM_ROW_TO_ROW(7,4,R1,0)			// 31 <-> 16

        N2 = result[0].rows()-1;
      N = result[4].cols();
      DSSUM_COL_TO_COL_NO_EDGES(0,4,0,N-1,1) 		// 0 <-> 16
        N2 = result[2].rows();
      DSSUM_COL_TO_COL(2,5,0,N-1)					// 2 <-> 17

        N  = result[7].cols();
      N2 = result[7].cols();
      DSSUM_COL_TO_ROW_BACKWARDS(1,6,0,N2-1)		// 30 <->  1
        DSSUM_COL_TO_ROW_BACKWARDS(0,7,0,N2-1)		// 31 <->  0

        N  = result[4].cols();
      N2 = result[0].rows();
      result[4][N-1][0] = result[7][N-1][N2-1] = result[ 0][0][0]; // lower left

      N = result[0].cols();
      dssumColToColMPI(result[1],row,row2,N-1,pos);	//  1 <-> 12
      dssumColToColMPI(result[3],row,row2,N-1,pos);	//  3 <-> 14
      dssumRowToRowMPI(result[6],row,row2,  0,pos);	// 29 <-> 30
#ifdef HAS_MPI
      MPI_Status res;
      MPI_Sendrecv( row,pos,MPI_PACKED,3,tag+1,
          row3,pos,MPI_PACKED,3,tag+1,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssumColToColMPI(result[1],row3,row2,N-1,pos,true);	//  1 <-> 12
      dssumColToColMPI(result[3],row3,row2,N-1,pos,true);	//  3 <-> 14
      dssumRowToRowMPI(result[6],row3,row2,  0,pos,true);	// 29 <-> 30

      N = result[0].rows();
      pos = 0;
      dssumRowToRowMPI(result[2],row,row2,N-1,pos);	//  2 <->  4
      dssumRowToRowMPI(result[3],row,row2,N-1,pos);	//  3 <->  5
      dssumRowToRowMPI(result[5],row,row2, R1,pos);	// 17 <-> 18
#ifdef HAS_MPI
      MPI_Sendrecv( row,pos,MPI_PACKED,1,tag+2,
          row3,pos,MPI_PACKED,1,tag+2,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssumRowToRowMPI(result[2],row3,row2,N-1,pos,true);	//  2 <->  4
      dssumRowToRowMPI(result[3],row3,row2,N-1,pos,true);	//  3 <->  5
      dssumRowToRowMPI(result[5],row3,row2, R1,pos,true);	// 17 <-> 18
    } else if( rank == 1 ) {
      DSSUM_COL_TO_COL( 0, 1,N-1,0) 			 	//  4 <->  5
        DSSUM_COL_TO_COL( 2, 3,N-1,0) 			 	//  6 <->  7
        N = result[0].rows();
      DSSUM_ROW_TO_ROW( 0, 2,N-1,0) 				//  4 <->  6
        DSSUM_ROW_TO_ROW( 1, 3,N-1,0) 				//  5 <->  7
        int R1 = result[0].rows()-1;
      for( int i=4;i<7;++i ) {
        DSSUM_ROW_TO_ROW(i,i+1,R1,0)			// 18..20 <-> 19..21
      }
      N2 = result[2].rows();
      N  = result[4].cols();
      DSSUM_COL_TO_COL(0,4,0,N-1)					//  4 <-> 18
        N2 = result[0].rows()-1;
      DSSUM_COL_TO_COL_NO_EDGES(2,5,0,N-1,0)		//  6 <-> 19

        N  = result[2].rows();
      N2 = result[6].rows();
      int N3 = result[6].cols();
      DSSUM_COL_TO_ROW( 2, 6,N-1,N3-1)			//  6 <-> 20
        DSSUM_COL_TO_ROW( 3, 7,N-1,N3-1)			//  7 <-> 21

        N  = result[6].cols();
      N2 = result[0].rows();
      result[6][N-1][0] = result[5][N-1][N2-1] = result[2][0][N2-1]; // lower right

      pos = 0;
      N = result[0].cols();
      dssumColToColMPI(result[1],row,row2,N-1,pos);	//  5 <->  8
      dssumColToColMPI(result[3],row,row2,N-1,pos);	//  7 <-> 10
      dssumRowToRowMPI(result[7],row,row2, R1,pos);	// 21 <-> 22
#ifdef HAS_MPI
      MPI_Status res;
      MPI_Sendrecv( row,pos,MPI_PACKED,2,tag+3,
          row3,pos,MPI_PACKED,2,tag+3,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssumColToColMPI(result[1],row3,row2,N-1,pos,true);	//  5 <->  8
      dssumColToColMPI(result[3],row3,row2,N-1,pos,true);	//  7 <-> 10
      dssumRowToRowMPI(result[7],row3,row2, R1,pos,true);	// 21 <-> 22

      pos = 0;
      dssumRowToRowMPI(result[0],row,row2,0,pos);	//  2 <->  4
      dssumRowToRowMPI(result[1],row,row2,0,pos);	//  3 <->  5
      dssumRowToRowMPI(result[4],row,row2,0,pos);	// 17 <-> 18
#ifdef HAS_MPI
      MPI_Sendrecv( row,pos,MPI_PACKED,0,tag+2,
          row3,pos,MPI_PACKED,0,tag+2,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssumRowToRowMPI(result[0],row,row2,0,pos,true);	//  2 <->  4
      dssumRowToRowMPI(result[1],row,row2,0,pos,true);	//  3 <->  5
      dssumRowToRowMPI(result[4],row,row2,0,pos,true);	// 17 <-> 18
    } else if( rank == 2 ) {
      DSSUM_COL_TO_COL( 0, 1,N-1,0) 				//  8 <->  9
        DSSUM_COL_TO_COL( 2, 3,N-1,0) 				// 10 <-> 11
        N = result[0].rows();
      DSSUM_ROW_TO_ROW( 0, 2,N-1,0) 				//  8 <-> 10
        DSSUM_ROW_TO_ROW( 1, 3,N-1,0) 				//  9 <-> 11
        N = result[0].rows();
      int R1 = result[0].rows()-1;
      for( int i=4;i<7;++i ) {
        DSSUM_ROW_TO_ROW(i,i+1,R1,0)			// 22..24 <-> 23..25
      }
      N2 = result[4].rows();
      int N3 = result[4].cols();
      N = result[0].rows();
      DSSUM_COL_TO_ROW( 2, 4,N-1,N3-1)			// 10 <-> 22
        DSSUM_COL_TO_ROW( 3, 5,N-1,N3-1)			// 11 <-> 23
        N2 = result[3].cols();
      N = result[6].cols();
      DSSUM_COL_TO_COL_BACKWARDS_NO_EDGES(3,6,N2-1,N-1,0)	// 11 <-> 24
        DSSUM_COL_TO_COL_BACKWARDS( 1,7,N2-1,N-1)		 	//  9 <-> 25

        N  = result[4].cols();
      N2 = result[0].rows();
      result[5][N-1][N2-1] = result[6][N-1][0] = result[3][N2-1][N2-1]; // upper right

      dssumColToColMPI(result[0],row,row2,0,pos);	//  5 <->  8
      dssumColToColMPI(result[2],row,row2,0,pos);	//  7 <-> 10 
      dssumRowToRowMPI(result[4],row,row2,0,pos);	// 21 <-> 22
#ifdef HAS_MPI
      MPI_Status res;
      MPI_Sendrecv( row,pos,MPI_PACKED,1,tag+3,
          row3,pos,MPI_PACKED,1,tag+3,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssumColToColMPI(result[0],row3,row2,0,pos,true);	//  5 <->  8
      dssumColToColMPI(result[2],row3,row2,0,pos,true);	//  7 <-> 10 
      dssumRowToRowMPI(result[4],row3,row2,0,pos,true);	// 21 <-> 22

      pos = 0;
      N = result[7].rows();
      dssumRowToRowMPI(result[0],row,row2,  0,pos);	// 14 <->  8
      dssumRowToRowMPI(result[1],row,row2,  0,pos);	// 15 <->  9
      dssumRowToRowMPI(result[7],row,row2,N-1,pos);	// 25 <-> 26 
#ifdef HAS_MPI
      MPI_Sendrecv( row,pos,MPI_PACKED,3,tag+4,
          row3,pos,MPI_PACKED,3,tag+4,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssumRowToRowMPI(result[0],row3,row2,  0,pos,true);	// 14 <->  8
      dssumRowToRowMPI(result[1],row3,row2,  0,pos,true);	// 15 <->  9
      dssumRowToRowMPI(result[7],row3,row2,N-1,pos,true);	// 25 <-> 26 
    } else {
      DSSUM_COL_TO_COL( 0, 1,N-1,0) 				// 12 <-> 13
        DSSUM_COL_TO_COL( 2, 3,N-1,0) 				// 14 <-> 15

        N = result[0].rows();
      DSSUM_ROW_TO_ROW( 0, 2,N-1,0) 				// 12 <-> 14
        DSSUM_ROW_TO_ROW( 1, 3,N-1,0) 				// 13 <-> 15

        int R1 = result[0].rows()-1;
      N = result[4].rows();
      for( int i=4;i<7;++i ) {
        DSSUM_ROW_TO_ROW(i,i+1,R1,0)			// 26..28 <-> 27..29
      }

      N2 = result[3].cols();
      N = result[4].cols();
      DSSUM_COL_TO_COL_BACKWARDS(3,4,N2-1,N-1)			// 15 <-> 26
        DSSUM_COL_TO_COL_BACKWARDS_NO_EDGES(1,5,N2-1,N-1,1) // 13 <-> 27

        N2 = result[6].cols();
      DSSUM_COL_TO_ROW_BACKWARDS(1,6,0,N2-1)		// 13 <-> 28
        DSSUM_COL_TO_ROW_BACKWARDS(0,7,0,N2-1)		// 12 <-> 29

        N  = result[7].cols();
      N2 = result[0].rows();
      result[5][N-1][N2-1] = result[6][N-1][0] = result[1][N2-1][0]; // upper left

      dssumColToColMPI(result[0],row,row2, 0,pos);	//  1 <-> 12
      dssumColToColMPI(result[2],row,row2, 0,pos);	//  3 <-> 14
      dssumRowToRowMPI(result[7],row,row2,R1,pos); 	// 29 <-> 30
#ifdef HAS_MPI
      MPI_Status res;
      MPI_Sendrecv( row,pos,MPI_PACKED,0,tag+1,
          row3,pos,MPI_PACKED,0,tag+1,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssumColToColMPI(result[0],row,row2, 0,pos,true);	//  1 <-> 12
      dssumColToColMPI(result[2],row,row2, 0,pos,true);	//  3 <-> 14
      dssumRowToRowMPI(result[7],row,row2,R1,pos,true); 	// 29 <-> 30

      pos = 0;
      N = result[2].rows();
      dssumRowToRowMPI(result[2],row,row2,N-1,pos); 	// 14 <->  8
      dssumRowToRowMPI(result[3],row,row2,N-1,pos);	// 15 <->  9
      dssumRowToRowMPI(result[4],row,row2,  0,pos);	// 25 <-> 26 
#ifdef HAS_MPI
      MPI_Sendrecv( row,pos,MPI_PACKED,2,tag+4,
          row3,pos,MPI_PACKED,2,tag+4,MPI_COMM_WORLD,&res);
#endif
      pos = 0;
      dssumRowToRowMPI(result[2],row,row2,N-1,pos,true); 	// 14 <->  8
      dssumRowToRowMPI(result[3],row,row2,N-1,pos,true);	// 15 <->  9
      dssumRowToRowMPI(result[4],row,row2,  0,pos,true);	// 25 <-> 26 
    }
  }
  if( size > 1 )  {
    delete[] row;
    delete[] row2;
    delete[] row3;
  }
}

basics::coarseDescriptor bigCircleGeometry::getRestrictionGridInfo() const
{
  basics::coarseDescriptor result;
  result.size1 = 3;
  result.size2 = 113;
  result.size3 = 32;
  result.size4 = 145;

  return( result );
}

#define ADDELEMENTS(group,add) \
  for( int i=0;i<m_size2;++i ) \
desc.elements.push_back((group)*m_size2+i+add);

basics::coarseGrid bigCircleGeometry::getDivisionInfo(int size) const
{
  basics::coarseGrid result;
  basics::coarseDescriptor desc;
  if( 12 % size ) {
    cout << "invalid MPI comm size" << endl;
    exit(1);
  }
  for( int j=0;j<size;++j ) {
    for( int k=0;k<12/size;++k )
      ADDELEMENTS(k+j*12/size,0);
    result.push_back(desc);
    desc.elements.clear();
  }

  return( result );
}

void bigCircleGeometry::setupRestrictionOperators(vector<basics::Matrix*>& result,
    const basics::Vector& coarseGrid,
    const vector<const basics::Vector*> grids) const
{
  result.push_back(new basics::Matrix(utilities::GLL::interpolationMatrix(coarseGrid,*grids[0])));
  result.push_back(new basics::Matrix(utilities::GLL::interpolationMatrix(coarseGrid,*grids[1])));
  result.push_back(new basics::Matrix(utilities::GLL::interpolationMatrix(coarseGrid,*grids[2])));
  result.push_back(new basics::Matrix(utilities::GLL::interpolationMatrix(coarseGrid,*grids[3])));
  result.push_back(new basics::Matrix(utilities::GLL::interpolationMatrix(coarseGrid,*grids[4]))); 
}

void bigCircleGeometry::getRestriction(basics::Vector& result,
    basics::matrixStack& work,
    const basics::matrixStack& LG,
    const basics::matrixStack& input,
    basics::Matrix& temp,
    const basics::coarseGrid& group,
    const vector<basics::Matrix*>& RT)
{
  for( int n=0;n<16;++n ) {
    basics::multTranspose(temp,input[n], *RT[1], 'N','N');
    basics::multTranspose(work[n],*RT[0],temp,   'T','N');
  }
  for( int n=16;n<32;++n ) {
    basics::multTranspose(temp,input[n],*RT[3], 'N','N');
    basics::multTranspose(work[n],*RT[2],temp,  'T','N');
  }
  dssumCoarse(work,group);
  localToGlobal(result,work,LG);
}

void bigCircleGeometry::getRestriction(basics::Vector& result,
    basics::matricesStack& work,
    basics::matricesStack& work2,
    const basics::matricesStack& LG,
    const basics::matricesStack& input,
    basics::Matrix& temp,
    const basics::coarseGrid& group,
    const vector<basics::Matrix*>& RT)
{
  int max=input.size();
  //#pragma omp parallel for schedule(static)
  for( int k=0;k<max;++k )
    basics::applyLocalGlobal(work2[k],input[k],*RT[4],'N','N');
  for( int l=0;l<work2[0].matrices();++l )
    getRestriction(result,work.at(l),LG.at(l),work2.at(l),temp,group,RT);
}

void bigCircleGeometry::gatherAndRestrict(basics::Vector& result,
    basics::matricesStack& work,
    basics::matricesStack& work2,
    basics::matricesStack& work3,
    const basics::matricesStack& LG,
    const basics::matricesStack& input,
    basics::Matrix& temp,
    const vector<basics::Matrix*>& RT,
    int rank, int size)
{
  vector<basics::coarseGrid> group = getCoarseGroups(size);
  if( size == 1 ) {
    getRestriction(result,work,work2,LG,input,temp,group[0],RT);
    return;
  }
#ifdef HAS_MPI
  int max = work2.size();
  for( int i=0;i<max;++i )
    basics::applyLocalGlobal(work2[i],input[i],*RT[4],'N','N');
  if( rank == 0 ) {
    for( unsigned int l=0;l<group.size();++l ) {
      for( unsigned int k=0;k<group[l].size();++k ) {
        if( l == rank ) {
          work3[group[l][k].size3] = work2[k];
        }
        else {
          MPI_Recv(work3[group[l][k].size3].data(),
              work3[group[l][k].size3].length(),
              MPI_DOUBLE,l,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
      }
    }
    basics::coarseGrid group1 = getCoarseGroups(1)[0];
    for( int l=0;l<work3[0].matrices();++l )
      getRestriction(result,work.at(l),LG.at(l),work3.at(l),temp,group1,RT);
  } else {
    int sends=group[rank].size();
    for( int i=0;i<sends;++i )
      MPI_Send(work2[i].data(),work2[i].length(),MPI_DOUBLE,0,0,
          MPI_COMM_WORLD);
  }
#endif
}

void bigCircleGeometry::getProlongiation(basics::matrixStack& result,
    basics::matrixStack& work,
    const basics::matrixStack& LG,
    const basics::Vector& input,
    basics::Matrix& temp,
    const vector<basics::Matrix*>& RT)
{
  globalToLocal(work,input,LG);
  doProlongiation(result,work,temp,RT);
}

void bigCircleGeometry::doProlongiation(basics::matrixStack& result,
    basics::matrixStack& work,
    basics::Matrix& temp,
    const vector<basics::Matrix*>& RT)
{
  for( int n=0;n<16;++n ) {
    basics::multTranspose(temp,*RT[0],work[n], 'N','N');
    basics::multTranspose(result[n],temp,*RT[1], 'N','T');
  }
  for( int n=16;n<32;++n ) {
    basics::multTranspose(temp,*RT[2],work[n], 'N','N');
    basics::multTranspose(result[n],temp,*RT[3], 'N','T');
  }
}

void bigCircleGeometry::getProlongiation(basics::matricesStack& result,
    basics::matricesStack& work,
    basics::matricesStack& work2,
    const basics::matricesStack& LG,
    const basics::Vector& input,
    basics::Matrix& temp,
    const vector<basics::Matrix*>& RT)
{
  globalToLocal(work,input,LG);
  for( int l=0;l<work[0].matrices();++l )
    doProlongiation(work2.at(l),work.at(l),temp,RT);
  int max=work2.size();
  //#pragma omp parallel for schedule(static)
  for( int k=0;k<max;++k )
    applyLocalGlobal(result[k],work2[k],*RT[4],'N','T');
}

void bigCircleGeometry::prolongAndScatter(basics::matricesStack& result,
    basics::matricesStack& work,
    basics::matricesStack& work2,
    basics::matricesStack& work3,
    const basics::matricesStack& LG,
    const basics::Vector& input,
    basics::Matrix& temp,
    const vector<basics::Matrix*>& RT,
    int rank, int size)
{
  if( size == 1 ) {
    getProlongiation(result,work,work2,LG,input,temp,RT);
    return;
  }
#ifdef HAS_MPI
  vector<basics::coarseGrid> group = getCoarseGroups(size);
  if( rank == 0 ) {
    globalToLocal(work,input,LG);
    for( int l=0;l<work[0].matrices();++l )
      doProlongiation(work3.at(l),work.at(l),temp,RT);
    for( unsigned int l=0;l<group.size();++l) {
      for( unsigned int k=0;k<group[l].size();++k){
        if( l == rank) {
          work2[k] = work3[group[l][k].size3];
        }
        else {
          MPI_Send(work3[group[l][k].size3].data(),
              work3[group[l][k].size3].length(),
              MPI_DOUBLE,l,0,MPI_COMM_WORLD);
        }
      }
    }
  } else {
    int recv=group[rank].size();
    for( int i=0;i<recv;++i )
      MPI_Recv(work2[i].data(),work2[i].length(),MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }
  int max=work2.size();
#pragma omp parallel for schedule(static)
  for( int k=0;k<max;++k )
    applyLocalGlobal(result[k],work2[k],*RT[4],'N','T');
#endif
}

bigCircleGeometry3D::bigCircleGeometry3D(const basics::Vector& weight,
    bigCircleGeometry& geom2D) :
  geometryStack3D(weight), m_geometry2D(geom2D)
{
  HDF5::HDF5Reader reader("bigcircle.hdf5");
  reader.read(*this,"grid");
  int max=m_grid.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i )
    m_grid[i]->getGH().computeJacobian();
  m_level = max / 768;
  m_perlevel = 768;
  if( max % 768 ) {
    m_level = max / 192;
    m_perlevel = 192;
    if( max % 192 ) {
      m_level = max / 48;
      m_perlevel = 48;
    }
  }
  m_division = (int)log2(m_perlevel/12); 
  m_size1 = (int)pow(2.f,m_division/2);
  m_size2 = m_size1*m_size1;
  computeMultiplicities();
  computeMass();
  computeReferenceGeometryDerivatives();
}

bigCircleGeometry3D::~bigCircleGeometry3D()
{
}

void bigCircleGeometry3D::dssum(basics::matrixStack& op, int rank, int size) const
{
  m_geometry2D.dssum(op,rank,size);
}

void bigCircleGeometry3D::dssum(basics::matricesStack& op, int rank, int size) const
{
  std::vector<basics::matricesStack*> ops = setupFakeStacks(op);

  for( int i=0;i<ops.size();++i )
    m_geometry2D.dssum(*ops[i]);//,rank,size);

  if( m_level > 1 ) {
    int N = op[0].matrices();
    for( int i=0;i<ops.size()-1;++i ) {
      ops[i]->at(N-1) += ops[i+1]->at(0);
      ops[i+1]->at(0) = ops[i]->at(N-1);
    }
    killFakeStacks(ops);
  }
}

void bigCircleGeometry3D::mask(basics::matrixStack& op, int rank, int size) const
{
  m_geometry2D.mask(op,rank,size);
}

void bigCircleGeometry3D::mask(basics::matricesStack& op, int rank, int size) const
{
  std::vector<basics::matricesStack*> ops = setupFakeStacks(op,m_level,m_perlevel/size);

  int N = (*ops[0])[0].matrices();
  for( int i=0;i<ops.size();++i )
    for( int l=0;l<N;++l )
      m_geometry2D.mask((*ops[i]).at(l),rank,size);

  ops[0]->at(0).clear();
  ops[ops.size()-1]->at(N-1).clear();
  killFakeStacks(ops);
}

basics::coarseGrid bigCircleGeometry3D::getDivisionInfo(int size) const
{
  basics::coarseGrid result;
  basics::coarseDescriptor desc;
  if( 12 % size ) {
    cout << "invalid MPI comm size" << endl;
    exit(1);
  }
  for( int j=0;j<size;++j ) {
    for( int n=0;n<m_level;++n )
      for( int k=0;k<12/size;++k )
        ADDELEMENTS(k+j*12/size,n*m_perlevel);
    result.push_back(desc);
    desc.elements.clear();
  }

  return( result );
}

vector<basics::coarseGrid> bigCircleGeometry3D::getCoarseGroups(int size) const
{
  vector<basics::coarseGrid> result1 = m_geometry2D.getCoarseGroups(size);
  vector<basics::coarseGrid> result = result1;
  assert(result1.size());
  for( int j=1;j<m_level;++j ) {
    for( int l=0;l<result1[0].size();++l ) {
      basics::coarseDescriptor desc=result1[0][l];
      for( int i=0;i<desc.elements.size();++i )
        desc.elements[i] += j*m_perlevel;
      desc.size3 += 32*j;
      result[0].push_back(desc);
    }
  }

  return( result );
}

basics::matricesStack* bigCircleGeometry3D::getCoarseBuffer(int rows, int cols, int matrices,
    const basics::coarseGrid& desc,
    const string& name, bool L2) const
{
  return m_geometry2D.getCoarseBuffer(rows,cols,matrices,desc,name,L2);
}

void bigCircleGeometry3D::fineToCoarse(basics::matricesStack& result,
    const basics::matricesStack& u,
    const basics::coarseGrid& desc,
    bool neumann) const
{
  vector<basics::matricesStack*> ops = setupFakeStacks(result,m_level,32);
  vector<basics::matricesStack*> uops = setupFakeStacks(const_cast<basics::matricesStack&>(u));

  for( int i=0;i<ops.size();++i )
    m_geometry2D.fineToCoarse(*ops[i],*uops[i],desc,neumann);

  killFakeStacks(ops);
  killFakeStacks(uops);
}

void bigCircleGeometry3D::fineToCoarseL2(basics::matricesStack& result,
    const basics::matricesStack& p,
    const basics::coarseGrid& desc) const
{
  vector<basics::matricesStack*> ops = setupFakeStacks(result,m_level,32);
  vector<basics::matricesStack*> uops = setupFakeStacks(const_cast<basics::matricesStack&>(p));

  for( int i=0;i<ops.size();++i )
    m_geometry2D.fineToCoarseL2(*ops[i],*uops[i],desc);

  killFakeStacks(ops);
  killFakeStacks(uops);
}

void bigCircleGeometry3D::fineToCoarseRestriction(basics::matricesStack& result,
    basics::matricesStack& buffer,
    const basics::matricesStack& u,
    const basics::coarseGrid& desc,
    int rank, int size) const
{
  vector<basics::matricesStack*> ops = setupFakeStacks(result,m_level,32);
  vector<basics::matricesStack*> uops = setupFakeStacks(const_cast<basics::matricesStack&>(u));
  vector<basics::matricesStack*> bops = setupFakeStacks(buffer);
  for( int i=0;i<ops.size();++i )
    m_geometry2D.fineToCoarseRestriction(*ops[i],*bops[i],*uops[i],desc,rank,size);

  killFakeStacks(ops);
  killFakeStacks(uops);
  killFakeStacks(bops);
}

void bigCircleGeometry3D::coarseToFine(basics::matricesStack& result,
    const basics::matricesStack& u,
    const basics::coarseGrid& desc,
    bool neumann) const
{
  vector<basics::matricesStack*> ops = setupFakeStacks(const_cast<basics::matricesStack&>(u),
      m_level,32);
  vector<basics::matricesStack*> uops = setupFakeStacks(result);
  for( int i=0;i<ops.size();++i )
    m_geometry2D.coarseToFine(*uops[i],*ops[i],desc,neumann);

  killFakeStacks(ops);
  killFakeStacks(uops);
}

void bigCircleGeometry3D::coarseToFineL2(basics::matricesStack& result,
    const basics::matricesStack& p,
    const basics::coarseGrid& desc) const
{
  vector<basics::matricesStack*> ops = setupFakeStacks(const_cast<basics::matricesStack&>(p),
      m_level,32);
  vector<basics::matricesStack*> uops = setupFakeStacks(result);
  for( int i=0;i<ops.size();++i )
    m_geometry2D.coarseToFineL2(*uops[i],*ops[i],desc);

  killFakeStacks(ops);
  killFakeStacks(uops);
}

void bigCircleGeometry3D::dssumCoarse(basics::matricesStack& result,
    const basics::coarseGrid& groups,
    int rank, int size, int tag) const
{
  std::vector<basics::matricesStack*> ops = setupFakeStacks(result,m_level,32);
  for( int i=0;i<ops.size();++i )
    m_geometry2D.dssumCoarse(*ops[i],groups,rank,size,tag);	
  int N = (*ops[0])[0].matrices();
  for( int i=0;i<ops.size()-1;++i ) {
    ops[i]->at(N-1) += ops[i+1]->at(0);
    ops[i+1]->at(0)  = ops[i]->at(N-1);
  }
  killFakeStacks(ops);
}

basics::coarseDescriptor bigCircleGeometry3D::getRestrictionGridInfo() const
{
  basics::coarseDescriptor result;
  result.size1 = 3;
  result.size2 = 113*(2*m_level-1);
  result.size3 = 32*m_level;
  result.size4 = 145*(2*m_level+1);

  return( result );
}

void bigCircleGeometry3D::getRestriction(basics::Vector& result,
    basics::matricesStack& work,
    basics::matricesStack& work2,
    const basics::matricesStack& LG,
    const basics::matricesStack& input,
    basics::Matrix& temp,
    const basics::coarseGrid& group,
    const vector<basics::Matrix*>& RT)
{
  int max=input.size();
  for( int k=0;k<max;++k )
    basics::applyLocalGlobal(work2[k],input[k],*RT[4],'N','N');
  vector<basics::matricesStack*> ops = setupFakeStacks(work2,m_level,32);
  vector<basics::matricesStack*> wops = setupFakeStacks(work,m_level,32);
  vector<basics::matricesStack*> LGs = setupFakeStacks(const_cast<basics::matricesStack&>(LG),
      m_level,32);
  for( int i=0;i<ops.size()-1;++i ) {
    ops[i]->at(2)   += ops[i+1]->at(0);
    ops[i+1]->at(0)  = ops[i]->at(2);
  }

  for( int i=0;i<ops.size();++i )
    for( int l=0;l<(*ops[i])[0].matrices();++l )
      m_geometry2D.getRestriction(result,wops[i]->at(l),LGs[i]->at(l),
          ops[i]->at(l),temp,group,RT);

  killFakeStacks(ops);
  killFakeStacks(wops);
  killFakeStacks(LGs);
}

void bigCircleGeometry3D::getProlongiation(basics::matricesStack& result,
    basics::matricesStack& work,
    basics::matricesStack& work2,
    const basics::matricesStack& LG,
    const basics::Vector& input,
    basics::Matrix& temp,
    const vector<basics::Matrix*>& RT)
{
  vector<basics::matricesStack*> ops = setupFakeStacks(work2,m_level,32);
  vector<basics::matricesStack*> wops = setupFakeStacks(work,m_level,32);
  vector<basics::matricesStack*> rops = setupFakeStacks(result,m_level,32);
  vector<basics::matricesStack*> LGs = setupFakeStacks(const_cast<basics::matricesStack&>(LG),
      m_level,32);
  for( int i=0;i<rops.size();++i )
    m_geometry2D.getProlongiation(*rops[i],*wops[i],*ops[i],*LGs[i],input,temp,RT);

  killFakeStacks(ops);
  killFakeStacks(wops);
  killFakeStacks(rops);
  killFakeStacks(LGs);
}

