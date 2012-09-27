#ifndef TWO_ELEMENT_H_ 
#define TWO_ELEMENT_H_

#include "mnl/geometry.h"

class twoElement : public basics::geometryStack {
  public:
    class edge1 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair((x+1)/4,0);
        }
    };

    class edge2 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(.5f,(x+1)/2);
        }
    };

    class edge3 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair((x+1)/4,1);
        }
    };

    class edge4 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(0,(x+1)/2);
        }
    };

    class edge5 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(.5f+(x+1)/4*3,0);
        }
    };

    class edge6 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(2,(x+1)/2);
        }
    };

    class edge7 : public basics::function1Dto2D {
      public:
        inline pair<Real,Real> val(Real x, Real t) const
        {
          return make_pair(.5f+(x+1)/4*3,1);
        }
    };

    twoElement(int N, int M, const basics::Vector& ggrid, const basics::Vector& weight) :
      geometryStack(ggrid,weight)
  {
    m_grid.push_back(new basics::spectralElement2D(N,M,m_g1,m_g2,m_g3,m_g4));
    m_grid.push_back(new basics::spectralElement2D(N,M,m_g5,m_g6,m_g7,m_g2));
    compute(0);
    computeMultiplicities();
    computeMass();
  }

    void boundary(basics::matrixStack& U,
        const basics::function2D& source, Real t) const
    {
      int N = U[0].rows();
      for( int alpha=0;alpha<m_ggrid.length();++alpha ) {
        /* bottom */
        pair<Real,Real> point = m_grid[0]->getGH().evaluate(m_ggrid[alpha],-1,m_ggrid);
        U[0][0][alpha] = source.val(point.first,point.second,t);
        point = m_grid[1]->getGH().evaluate(m_ggrid[alpha],-1,m_ggrid);
        U[1][0][alpha] = source.val(point.first,point.second,t);

        /* right */
        point = m_grid[1]->getGH().evaluate(1,m_ggrid[alpha],m_ggrid);
        U[1][alpha][N-1] = source.val(point.first,point.second,t);

        /* top */
        point = m_grid[0]->getGH().evaluate(m_ggrid[alpha],1,m_ggrid);
        U[0][N-1][alpha] = source.val(point.first,point.second,t);
        point = m_grid[1]->getGH().evaluate(m_ggrid[alpha],1,m_ggrid);
        U[1][N-1][alpha] = source.val(point.first,point.second,t);

        /* left */
        point = m_grid[0]->getGH().evaluate(-1,m_ggrid[alpha],m_ggrid);
        U[0][alpha][0] = source.val(point.first,point.second,t);
      }
    }

    void boundary(basics::matricesStack& U, const basics::function3D& source, Real t) const
    {
      int N = U[0].rows();
      for( int alpha=0;alpha<m_ggrid.length();++alpha ) {
        /* bottom */
        pair<Real,Real> point = m_grid[0]->getGH().evaluate(m_ggrid[alpha],-1,m_ggrid);
        for( int beta=0;beta<m_ggrid.length();++beta )
          U[0][beta][0][alpha] = source.val(point.first,point.second,
              m_ggrid[beta],t);
        point = m_grid[1]->getGH().evaluate(m_ggrid[alpha],-1,m_ggrid);
        for( int beta=0;beta<m_ggrid.length();++beta )
          U[1][beta][0][alpha] = source.val(point.first,point.second,m_ggrid[beta],t);

        /* right */
        point = m_grid[1]->getGH().evaluate(1,m_ggrid[alpha],m_ggrid);
        for( int beta=0;beta<m_ggrid.length();++beta )
          U[1][beta][alpha][N-1] = source.val(point.first,point.second,m_ggrid[beta],t);

        /* top */
        point = m_grid[0]->getGH().evaluate(m_ggrid[alpha],1,m_ggrid);
        for( int beta=0;beta<m_ggrid.length();++beta )
          U[0][beta][N-1][alpha] = source.val(point.first,point.second,m_ggrid[beta],t);
        point = m_grid[1]->getGH().evaluate(m_ggrid[alpha],1,m_ggrid);
        for( int beta=0;beta<m_ggrid.length();++beta )
          U[1][beta][N-1][alpha] = source.val(point.first,point.second,m_ggrid[beta],t);

        /* left */
        point = m_grid[0]->getGH().evaluate(-1,m_ggrid[alpha],m_ggrid);
        for( int beta=0;beta<m_ggrid.length();++beta )
          U[0][beta][alpha][0] = source.val(point.first,point.second,m_ggrid[beta],t);
      }

      /* add top and bottom */
      for( int i=0;i<m_grid.size();++i ) {
        m_grid[i]->evaluate(U[i][0],  m_ggrid,m_ggrid,source,-1,t);
        m_grid[i]->evaluate(U[i][N-1],m_ggrid,m_ggrid,source,1,t);
      }
    }

    virtual void dssum(basics::matrixStack& element) const
    {
      element[0].axpyRow(element[0].rows()-1,Real(1),element[1].row(0));
      element[1].copyRow(0,element[0].row(element[0].rows()-1));
    }

    virtual void mask(basics::matrixStack& op) const
    {
      int N = op[0].rows();

      op[0][0].clear();
      op[1][0].clear();

      op[1].clearRow(N-1);

      op[0][N-1].clear();
      op[1][N-1].clear();

      op[0].clearRow(0);
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
    edge1 m_g1;
    edge2 m_g2;
    edge3 m_g3;
    edge4 m_g4;
    edge5 m_g5;
    edge6 m_g6;
    edge7 m_g7;
};

#endif

