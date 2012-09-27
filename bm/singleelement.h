#ifndef SINGLE_ELEMENT_H_ 
#define SINGLE_ELEMENT_H_

#include "mnl/geometry.h"

/* this is case 1, figure 6 in tormods master thesis */
class tormod1Bottom: public basics::function1Dto2D {
  public:
    inline pair<Real,Real> val(Real x, Real t) const
    {
      return make_pair((x+1)/2.f*2.f,0);
    }
};

class tormod1Right: public basics::function1Dto2D {
  public:
    inline pair<Real,Real> val(Real x, Real t) const
    {
      return make_pair(2.f,(x+1)/2.f*(3.f));
    }
};

class tormod1Top: public basics::function1Dto2D {
  public:
    inline pair<Real,Real> val(Real x, Real t) const
    {
      return make_pair((x+1)/2.f*2.f,(3.f)-cos(M_PI*(x+1)/4));
    }
};

class tormod1Left: public basics::function1Dto2D {
  public:
    inline pair<Real,Real> val(Real x, Real t) const
    {
      return make_pair(0.f,(x+1)/2.f*2.f);
    }
};

/* this is case 2, figure 9 in tormods master thesis */
class tormod2Bottom: public basics::function1Dto2D {
  public:
    inline pair<Real,Real> val(Real x, Real t) const
    {
      return make_pair((x+1)/2.f*2.f,0);
    }
};

class tormod2Right : public basics::function1Dto2D {
  public:
    inline pair<Real,Real> val(Real x, Real t) const
    {
      return make_pair(2.f+(sqrt(2.f)/2.f-2.f)*(x+1)/2.f,(2-sqrt(2.f)/2.f)*(x+1)/2.f);
    }
};

class tormod2Top : public basics::function1Dto2D {
  public:
    inline pair<Real,Real> val(Real x, Real t) const
    {
      return make_pair((x+1)/2.f*sqrt(2.f)/2.f,(2.f+sin(3*M_PI/2+((x+1)/2.f*M_PI/4))));
    }
};

class tormod2Left : public basics::function1Dto2D {
  public:
    inline pair<Real,Real> val(Real x, Real t) const
    {
      return make_pair(0.f,(x+1)/2.f);
    }
};

/* unit rectangle */
//#define LX 2*M_PI
//#define LY 2

#define LX -1 //2*M_PI
#define LY -1

class Rectangle : public basics::function1Dto2D {
  public:
    Rectangle(Real Lx, Real Ly) :
      m_Lx(Lx), m_Ly(Ly)
  {
  }
  protected:
    Real m_Lx;
    Real m_Ly;
};

class rectangleBottom : public Rectangle {
  public:
    rectangleBottom(Real Lx=1, Real Ly=1) :
      Rectangle(Lx,Ly)
  {
  }

    inline pair<Real,Real> val(Real x, Real t) const
    {
      return make_pair(x,-1.f);
      //            return make_pair((x+1)*M_PI,-1.f);
      //            if( m_Lx != -1 )
      //                return make_pair((x+Real(1))/Real(2)*m_Lx,0);
      //            else
      //                return make_pair(x,Real(-1));
    }
};

class rectangleRight : public Rectangle {
  public:
    rectangleRight(Real Lx=1, Real Ly=1) :
      Rectangle(Lx,Ly)
  {
  }

    inline pair<Real,Real> val(Real x, Real t) const
    {
      return make_pair(1.f,x);
      //            return make_pair(2*M_PI,x);
      //            if( m_Lx != -1 )
      //                return make_pair(m_Lx,(x+1)/2.f*m_Ly);
      //            else
      //                return make_pair(Real(1),x);
    }
};

class rectangleTop : public Rectangle {
  public:
    rectangleTop(Real Lx=1, Real Ly=1) :
      Rectangle(Lx,Ly)
  {
  }

    inline pair<Real,Real> val(Real x, Real t) const
    {
      return make_pair(x,1.f);
      //            return make_pair((x+1)*M_PI,1.f);
      //            if( m_Lx != -1 )
      //                return make_pair((x+1)/2.f*(m_Lx),m_Ly);
      //            else
      //                return make_pair(x,Real(1));
    }
};

class rectangleLeft : public Rectangle {
  public:
    rectangleLeft(Real Lx=1, Real Ly=1) :
      Rectangle(Lx,Ly)
  {
  }

    inline pair<Real,Real> val(Real x, Real t) const
    {
      return make_pair(-1.f,x);
      //            return make_pair(0.f,x);
      //            if( m_Lx != -1 )
      //                return make_pair(0.f,(x+1)/2.f*m_Ly);
      //            else
      //                return make_pair(Real(-1),x);
    }
};

/* unit trapezoidal */
class trapezoidBottom: public Rectangle {
  public:
    trapezoidBottom(Real Lx=1, Real Ly=1) :
      Rectangle(Lx,Ly)
  {
  }

    inline pair<Real,Real> val(Real x, Real t) const
    {
      return make_pair((x+1)/2.f*m_Lx,0.f);
    }
};

class trapezoidRight : public Rectangle {
  public:
    trapezoidRight(Real Lx=1, Real Ly=1) :
      Rectangle(Lx,Ly)
  {
  }

    inline pair<Real,Real> val(Real x, Real t) const
    {
      return make_pair(m_Lx,(x+1)/2.f*(m_Ly+t)*2);
    }
};

class trapezoidTop : public Rectangle {
  public:
    trapezoidTop(Real Lx=1, Real Ly=1) :
      Rectangle(Lx,Ly)
  {
  }

    inline pair<Real,Real> val(Real x, Real t) const
    {
      return make_pair((x+1)/2.f*(m_Lx),m_Ly+(x+1)/2.f*(m_Ly+t));
    }
};

class trapezoidLeft : public Rectangle {
  public:
    trapezoidLeft(Real Lx=1, Real Ly=1) :
      Rectangle(Lx,Ly)
  {
  }

    inline pair<Real,Real> val(Real x, Real t) const
    {
      return make_pair(0.f,(x+1)/2.f*m_Ly);
    }
};

class singleElement : public basics::geometryStack {
  public:
    singleElement(int N, int M, const basics::Vector& ggrid) : 
      geometryStack(ggrid),
      m_g1(LX,LY), m_g2(LX,LY), m_g3(LX,LY), m_g4(LX,LY)
  {
    m_grid.push_back(new basics::spectralElement2D(N,M,m_g1,m_g2,m_g3,m_g4));
    compute(0);
    computeMultiplicities();
  }

    virtual ~singleElement()
    {
    }

    void boundary(basics::matricesStack& U, const basics::function3D& source, Real t)
    {
      int N = U[0].rows();
      for( int alpha=0;alpha<m_ggrid.length();++alpha ) {
        /* bottom */
        pair<Real,Real> point = m_grid[0]->getGH().evaluate(m_ggrid[alpha],-1,m_ggrid);
        for( int beta=0;beta<m_ggrid.length();++beta)
          U[0][beta][0][alpha] = source.val(point.first,point.second,
              m_ggrid[beta],t);

        /* right */
        point = m_grid[0]->getGH().evaluate(1,m_ggrid[alpha],m_ggrid);
        for( int beta=0;beta<m_ggrid.length();++beta)
          U[0][beta][alpha][N-1] = source.val(point.first,point.second,
              m_ggrid[beta],t);

        /* top */
        point = m_grid[0]->getGH().evaluate(m_ggrid[alpha],1,m_ggrid);
        for( int beta=0;beta<m_ggrid.length();++beta)
          U[0][beta][N-1][alpha] = source.val(point.first,point.second,
              m_ggrid[beta],t);

        /* left */
        point = m_grid[0]->getGH().evaluate(-1,m_ggrid[alpha],m_ggrid);
        for( int beta=0;beta<m_ggrid.length();++beta)
          U[0][beta][alpha][0] = source.val(point.first,point.second,
              m_ggrid[beta],t);
      }

      /* add top and bottom */
      for( int i=0;i<m_grid.size();++i ) {
        m_grid[i]->evaluate(U[i][0],  m_ggrid,m_ggrid,source,m_ggrid[0],  t);
        m_grid[i]->evaluate(U[i][N-1],m_ggrid,m_ggrid,source,m_ggrid[N-1],t);
      }
    }

    virtual void dssum(basics::matrixStack& element) const
    {
    }

    virtual void mask(basics::matrixStack& op) const
    {
      int N = op[0].rows();

      op[0].clearRow(N-1);
      op[0].clearRow(0);
      op[0][0].clear();
      op[0][N-1].clear();
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
    rectangleBottom m_g1;
    rectangleRight m_g2;
    rectangleTop m_g3;
    rectangleLeft m_g4;
    //        tormod1Bottom m_g1;
    //        tormod1Right m_g2;
    //        tormod1Top m_g3;
    //        tormod1Left m_g4;
};

#endif

