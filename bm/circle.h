#ifndef CIRCLE_H_
#define CIRCLE_H_

#include "mnl/geometry.h"

#define RADIUS Real(1)

#define XLEN RADIUS/Real(3)*cos(M_PI/4)
#define YLEN RADIUS/Real(3)*sin(M_PI/4)

class circleGeometry : public basics::geometryStack {
  public:
    class edge1 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(x*XLEN,-YLEN);
        }
    };

    class edge2 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(XLEN,x*YLEN);
        }
    };

    class edge3 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(x*XLEN,YLEN);
        }
    };

    class edge4 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(-XLEN,x*YLEN);
        }
    };

    class edge5 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair( XLEN + Real(2)*RADIUS/Real(3)*cos(M_PI/4)*(x+1)/Real(2),
              -YLEN - Real(2)/Real(3)*RADIUS*sin(M_PI/4)*(x+1)/Real(2));
        }
    };

    class edge52 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(RADIUS*cos(3*M_PI/2+M_PI/4)*(2-x)/3,
              RADIUS*sin(3*M_PI/2+M_PI/4)*(2-x)/3);
        }
    };

    class edge6 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(RADIUS*cos(M_PI/4*x),RADIUS*sin(M_PI/4*x));
        }
    };

    class edge7 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(XLEN + Real(2)/Real(3)*RADIUS*(x+1)/Real(2)*cos(M_PI/4),
              YLEN + Real(2)/Real(3)*RADIUS*(x+1)/Real(2)*sin(M_PI/4));
        }
    };

    class edge8 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(RADIUS*cos(M_PI/2-M_PI/4*x),RADIUS*sin(M_PI/2-M_PI/4*x));
        }
    };

    class edge9 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(-XLEN - Real(2)/Real(3)*RADIUS*(x+1)/Real(2)*cos(M_PI/4),
              YLEN + Real(2)/Real(3)*RADIUS*(x+1)/Real(2)*sin(M_PI/4));
        }
    };

    class edge92 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(RADIUS*cos(M_PI/2+M_PI/4)*(2-x)/3,
              RADIUS*sin(M_PI/2+M_PI/4)*(2-x)/3);
        }
    };

    class edge10 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(RADIUS*cos(M_PI-M_PI/4*x),RADIUS*sin(M_PI-M_PI/4*x));
        }
    };

    class edge11 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(RADIUS*cos(M_PI+M_PI/4)+(x+1)/2*cos(M_PI/4)*2*RADIUS/3,
              RADIUS*sin(M_PI+M_PI/4)+(x+1)/2*sin(M_PI/4)*2*RADIUS/3);
        }
    };

    class edge12 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(RADIUS*cos(3*M_PI/2+M_PI/4*x),RADIUS*sin(3*M_PI/2+M_PI/4*x));
        }
    };

    circleGeometry(int N, int M, const basics::Vector& weight) :
      geometryStack(weight)
  {
    m_grid.push_back(new basics::spectralElement2D(N,M,m_g1,m_g2,m_g3,m_g4));
    m_grid.push_back(new basics::spectralElement2D(N,M,m_g5,m_g6,m_g7,m_g2));
    m_grid.push_back(new basics::spectralElement2D(N,M,m_g3,m_g7,m_g8,m_g9));
    m_grid.push_back(new basics::spectralElement2D(N,M,m_g11,m_g4,m_g92,m_g10));
    m_grid.push_back(new basics::spectralElement2D(N,M,m_g12,m_g52,m_g1,m_g11));
    compute(0);
    computeMultiplicities();
  }

    virtual ~circleGeometry()
    {
    }

    void boundary(basics::matricesStack& U, const basics::function3D& source, Real t) const
    {
      int N = U[0].rows();
      for( int alpha=0;alpha<m_ggrid.length();++alpha ) {
        /* g6 */
        pair<Real,Real> point = m_grid[1]->getGH().evaluate(1,m_ggrid[alpha],m_ggrid);
        for( int beta=0;beta<m_ggrid.length();++beta)
          U[1][beta][alpha][N-1] = source.val(point.first,point.second,m_ggrid[beta],t);

        /* g8 */
        point = m_grid[2]->getGH().evaluate(m_ggrid[alpha],1,m_ggrid);
        for( int beta=0;beta<m_ggrid.length();++beta)
          U[2][beta][N-1][alpha] = source.val(point.first,point.second,m_ggrid[beta],t);

        /* g10 */
        point = m_grid[3]->getGH().evaluate(-1,m_ggrid[alpha],m_ggrid);
        for( int beta=0;beta<m_ggrid.length();++beta)
          U[3][beta][alpha][0] = source.val(point.first,point.second,m_ggrid[beta],t);

        /* g12 */
        point = m_grid[4]->getGH().evaluate(m_ggrid[alpha],-1,m_ggrid);
        for( int beta=0;beta<m_ggrid.length();++beta)
          U[4][beta][0][alpha] = source.val(point.first,point.second,m_ggrid[beta],t);
      }

      /* add top and bottom */
      for( int i=0;i<m_grid.size();++i ) {
        m_grid[i]->evaluate(U[i][0],m_ggrid,m_ggrid,source,-1,t);
        m_grid[i]->evaluate(U[i][N-1],m_ggrid,m_ggrid,source,1,t);
      }
    }

    void boundary(basics::matrixStack& U, const basics::function2D& source, Real t) const
    {
      int N = U[0].rows();
      for( int alpha=0;alpha<m_ggrid.length();++alpha ) {
        /* g6 */
        pair<Real,Real> point = m_grid[1]->getGH().evaluate(1,m_ggrid[alpha],m_ggrid);
        U[1][alpha][N-1] = source.val(point.first,point.second,t);

        /* g8 */
        point = m_grid[2]->getGH().evaluate(m_ggrid[alpha],1,m_ggrid);
        U[2][N-1][alpha] = source.val(point.first,point.second,t);

        /* g10 */
        point = m_grid[3]->getGH().evaluate(-1,m_ggrid[alpha],m_ggrid);
        U[3][alpha][0] = source.val(point.first,point.second,t);

        /* g12 */
        point = m_grid[4]->getGH().evaluate(m_ggrid[alpha],-1,m_ggrid);
        U[4][0][alpha] = source.val(point.first,point.second,t);
      }
    }

    virtual void dssum(basics::matrixStack& element) const
    {
      int N = element[0].rows();
      int N2 = N-2;

      /* element 1 and 2 - last row <-> first row */
      AXPY_ROW_FROM_ROW_NO_EDGES(element[0].data()[0]+N-1,element[1].data()[0],1)
        COPY_ROW_FROM_ROW_NO_EDGES(element[1].data()[0],element[0].data()[0]+N-1,1)

        /* element 1 and 3 - last column <-> first column */
        AXPY_NO_EDGES(element[0].data()[N-1],element[2].data()[0],1)
        COPY_NO_EDGES(element[2].data()[0],element[0].data()[N-1],1)

        /* element 1 and 4 - first row <-> last row */
        AXPY_ROW_FROM_ROW_NO_EDGES(element[0].data()[0],element[3].data()[0]+N-1,1)
        COPY_ROW_FROM_ROW_NO_EDGES(element[3].data()[0]+N-1,element[0].data()[0],1)

        /* element 1 and 5 - first column <-> last column */
        AXPY_NO_EDGES(element[0].data()[0],element[4].data()[N-1],1)
        COPY_NO_EDGES(element[4].data()[N-1],element[0].data()[0],1)

        /* element 2 and 3 - last column <-> last row */
        AXPY_ROW_TO_COLUMN(element[1].data()[N-1],element[2].data()[0]+N-1)
        element[2].copyRow(N-1,element[1][N-1]);

      /* element 2 and 5  - first column <-> last row (backwards) */
      for( int j=0;j<N;++j ) {
        element[1][0][N-j-1] += element[4][j][N-1];
        element[4][j][N-1] = element[1][0][N-j-1];
      }

      /* element 3 and 4  - first row <-> last column (backwards) */
      for( int j=0;j<N;++j ) {
        element[2][j][0] += element[3][N-1][N-1-j];
        element[3][N-1][N-1-j] = element[2][j][0];
      }

      /* element 4 and 5  - first column <-> first row */
      AXPY_ROW_TO_COLUMN(element[3].data()[0],element[4].data()[0])
        element[4].copyRow(0,element[3][0]);

      /* fix the corners between 1 and the others */

      /* lower left */
      element[0][0][0] = element[4][N-1][0] = element[3][0][N-1] = 	 element[0][0][0]+
        element[4][N-1][0];

      /* lower right */
      element[0][0][N-1] = element[1][0][0] = element[4][N-1][N-1] =	 element[0][0][N-1]+
        element[1][0][0];

      /* upper right */
      element[0][N-1][N-1] = element[1][N-1][0] = element[2][0][N-1] = element[0][N-1][N-1]+
        element[1][N-1][0];

      /* upper left */
      element[0][N-1][0] = element[2][0][0] =	 element[3][N-1][N-1] =  element[0][N-1][0]+
        element[2][0][0];
    }

    virtual void mask(basics::matrixStack& op) const
    {
      int N = op[1].rows();

      /* element 2 - g6 */
      op[1].clearRow(N-1);

      /* element 3 - g8 */
      op[2][N-1].clear();

      /* element 4 - g10 */
      op[3].clearRow(0);

      /* element 5 - g12 */
      op[4][0].clear();
    }

    virtual void mask(basics::matricesStack& op) const
    {
      for( int k=0;k<op.size();++k ) {
        op[k][0].clear();
        op[k][op[0].matrices()-1].clear();
      }   
      for( int l=1;l<op[0].matrices()-1;++l )
        mask(op.at(l));
    }
  protected:
    edge1   m_g1;
    edge2   m_g2;
    edge3   m_g3;
    edge4   m_g4;
    edge5   m_g5;
    edge52 m_g52;
    edge6   m_g6;
    edge7   m_g7;
    edge8   m_g8;
    edge9   m_g9;
    edge92 m_g92;
    edge10 m_g10;
    edge11 m_g11;
    edge12 m_g12;
};

#endif

