#ifndef FUNCTION_2D_H_
#define FUNCTION_2D_H_

class shenTest2DX : public basics::function2D {
  public:
    shenTest2DX()
    {
      m_id = 2;
    }

    Real val(Real x, Real y, Real t) const
    {
      return( sin(2*M_PI*y)*sin(x)*sin(x)*sin(t) );
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return( sin(2*x)*sin(2*M_PI*y)*sin(t) );
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( 2*M_PI*cos(2*M_PI*y)*sin(t)*sin(x)*sin(x) );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return( sin(2*M_PI*y)*sin(x)*sin(x)*cos(t) );
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return( 2*cos(x)*cos(x)*sin(t)*sin(2*M_PI*y) - 2*sin(x)*sin(x)*sin(t)*sin(2*M_PI*y) );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return( -4*M_PI*M_PI*sin(2*M_PI*y)*sin(t)*sin(x)*sin(x) );
    }
};

class shenStationaryTest2DX : public basics::function2D {
  public:
    shenStationaryTest2DX()
    {
      m_id = 3;
    }

    Real val(Real x, Real y, Real t) const
    {
      return( sin(2*M_PI*y)*sin(x)*sin(x) );
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return( sin(2*x)*sin(2*M_PI*y) );
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( 2*M_PI*cos(2*M_PI*y)*sin(x)*sin(x) );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return 0;
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return( 2*cos(x)*cos(x)*sin(2*M_PI*y) - 2*sin(x)*sin(x)*sin(2*M_PI*y) );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return( -4*M_PI*M_PI*sin(2*M_PI*y)*sin(x)*sin(x) );
    }
};

class shenTest2DY : public basics::function2D {
  public:
    shenTest2DY()
    {
      m_id = 2;
    }

    Real val(Real x, Real y, Real t) const
    {
      return( -sin(Real(2)*x)*sin(M_PI*y)*sin(M_PI*y)*sin(t)/M_PI );
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return( Real(-2)*cos(Real(2)*x)*sin(t)*sin(M_PI*y)*sin(M_PI*y)/M_PI );
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( -sin(Real(2)*M_PI*y)*sin(t)*sin(Real(2)*x) );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return( -cos(t)*sin(Real(2)*x)*sin(M_PI*y)*sin(M_PI*y)/M_PI );
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return( Real(4)*sin(Real(2)*x)*sin(t)*sin(M_PI*y)*sin(M_PI*y)/M_PI );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return( Real(2)*M_PI*sin(M_PI*y)*sin(M_PI*y)*sin(t)*sin(Real(2)*x) - Real(2)*M_PI*cos(M_PI*y)*cos(M_PI*y)*sin(t)*sin(Real(2)*x) );
    }
};

class shenStationaryTest2DY : public basics::function2D {
  public:
    shenStationaryTest2DY()
    {
      m_id = 3;
    }

    Real val(Real x, Real y, Real t) const
    {
      return( -sin(Real(2)*x)*sin(M_PI*y)*sin(M_PI*y)/M_PI );
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return( Real(-2)*cos(Real(2)*x)*sin(M_PI*y)*sin(M_PI*y)/M_PI );
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( -sin(Real(2)*M_PI*y)*sin(Real(2)*x) );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return 0;
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return( Real(4)*sin(Real(2)*x)*sin(M_PI*y)*sin(M_PI*y)/M_PI );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return( Real(2)*M_PI*sin(M_PI*y)*sin(M_PI*y)*sin(Real(2)*x) - Real(2)*M_PI*cos(M_PI*y)*cos(M_PI*y)*sin(Real(2)*x) );
    }
};

class shenTest2DP : public basics::function2D {
  public:
    shenTest2DP()
    {
      m_id = 2;
    }

    Real val(Real x, Real y, Real t) const
    {
      return( cos(x)*cos(M_PI*y)*sin(t) );
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return( -sin(x)*cos(M_PI*y)*sin(t) );
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( -M_PI*cos(x)*sin(M_PI*y)*sin(t) );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return( cos(x)*cos(M_PI*y)*cos(t) );
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return( -cos(x)*cos(M_PI*y)*sin(t) );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return( -M_PI*M_PI*cos(x)*cos(M_PI*y)*sin(t) );
    }
};

class shenTest2DPD : public basics::function2D {
  public:
    shenTest2DPD()
    {
      m_id = 2;
    }

    Real val(Real x, Real y, Real t) const
    {
      return( cos(x)*sin(M_PI*y)*sin(t) );
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return( -sin(x)*sin(M_PI*y)*sin(t) );
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( M_PI*cos(x)*cos(M_PI*y)*sin(t) );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return( cos(x)*sin(M_PI*y)*cos(t) );
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return( -cos(x)*sin(M_PI*y)*sin(t) );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return( -M_PI*M_PI*cos(x)*sin(M_PI*y)*sin(t) );
    }
};

class shenStationaryTest2DP : public basics::function2D {
  public:
    shenStationaryTest2DP()
    {
      m_id = 3;
    }

    Real val(Real x, Real y, Real t) const
    {
      return( cos(x)*cos(M_PI*y) );
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return( -sin(x)*cos(M_PI*y) );
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( -M_PI*cos(x)*sin(M_PI*y) );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return( 0 );
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return( -cos(x)*cos(M_PI*y) );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return( -M_PI*M_PI*cos(x)*cos(M_PI*y) );
    }
};

class easyTest2DX : public basics::function2D {
  public:
    easyTest2DX()
    {
      m_id = 4;
    }

    Real val(Real x, Real y, Real t) const
    {
      return( cos(x)*sin(M_PI*y)*sin(t) );
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return( -sin(x)*sin(M_PI*y)*sin(t) );
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( M_PI*cos(x)*cos(M_PI*y)*sin(t) );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return( cos(x)*sin(M_PI*y)*cos(t) );
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return( -cos(x)*sin(M_PI*y)*sin(t) );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return( -M_PI*M_PI*cos(x)*sin(M_PI*y)*sin(t) );
    }
};

class easyTest2DY : public basics::function2D {
  public:
    easyTest2DY()
    {
      m_id = 4;
    }

    Real val(Real x, Real y, Real t) const
    {
      return( -sin(x)*(Real(1)+cos(M_PI*y))*sin(t)/M_PI );
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return( -cos(x)*(Real(1)+cos(M_PI*y))*sin(t)/M_PI );
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( sin(x)*sin(M_PI*y)*sin(t) );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return( -sin(x)*(Real(1)+cos(M_PI*y))*cos(t)/M_PI );
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return( sin(x)*(Real(1)+cos(M_PI*y))*sin(t)/M_PI );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return( M_PI*sin(x)*cos(M_PI*y)*sin(t) );
    }
};

class easyTest2DP : public basics::function2D {
  public:
    easyTest2DP()
    {
      m_id = 4;
    }

    Real val(Real x, Real y, Real t) const
    {
      return( cos(x)*cos(M_PI*y)*sin(t) );
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return( -sin(x)*cos(M_PI*y)*sin(t) );
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( -M_PI*cos(x)*sin(M_PI*y)*sin(t) );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return( cos(x)*cos(M_PI*y)*cos(t) );
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return( -cos(x)*cos(M_PI*y)*sin(t) );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return( -M_PI*M_PI*cos(x)*cos(M_PI*y)*sin(t) );
    }
};

//class einarTest2DX : public basics::function2D {
//public:
//    Real val(Real x, Real y, Real t) const
//    {
//        return( sin(x-t)*sin(x-t)*sin(Real(2)*M_PI*y) );
//    }

//    Real diffx(Real x, Real y, Real t) const
//    {
//        return( sin(Real(2)*(x-t))*sin(Real(2)*M_PI*y) );
//    }
//    
//    Real diffy(Real x, Real y, Real t) const
//    {
//        return( Real(2)*M_PI*sin(x-t)*sin(x-t)*cos(Real(2)*M_PI*y) );
//    }

//    Real difft(Real x, Real y, Real t) const
//    {
//        return( sin(2*(t-x))*sin(Real(2)*M_PI*y) );
//    }

//    Real diff2x(Real x, Real y, Real t) const
//    {
//        return( Real(2)*sin(Real(2)*M_PI*y)*(cos(x-t)*cos(x-t)-sin(x-t)*sin(x-t)) );
//    }

//    Real diff2y(Real x, Real y, Real t) const
//    {
//        return( -Real(4)*M_PI*M_PI*sin(Real(2)*M_PI*y)*sin(t-x)*sin(t-x) );
//    }
//};

//class einarTest2DY : public basics::function2D {
//    public:
//    Real val(Real x, Real y, Real t) const
//    {
//        return( -sin(Real(2)*(x-t))*sin(M_PI*y)*sin(M_PI*y)/M_PI );
//    }

//    Real diffx(Real x, Real y, Real t) const
//    {
//        return( -Real(2)/M_PI*cos(Real(2)*(x-t))*sin(M_PI*y)*sin(M_PI*y) );
//    }

//    Real diffy(Real x, Real y, Real t) const
//    {
//        return( Real(2)*sin(M_PI*y)*cos(M_PI*y)*sin(2*(t-x)) );
//    }

//    Real difft(Real x, Real y, Real t) const
//    {
//        return( Real(2)/M_PI*cos(Real(2)*(t-x))*sin(M_PI*y)*sin(M_PI*y) );
//    }

//    Real diff2x(Real x, Real y, Real t) const
//    {
//        return( Real(4)/M_PI*sin(Real(2)*(x-t))*sin(M_PI*y)*sin(M_PI*y) );
//    }

//    Real diff2y(Real x, Real y, Real t) const
//    {
//        return( Real(2)*M_PI*sin(Real(2)*(t-x))*(cos(M_PI*y)*cos(M_PI*y)-sin(M_PI*y)*sin(M_PI*y)) );
//    }
//};

//class einarTest2DP : public basics::function2D {
//    public:
//    Real val(Real x, Real y, Real t) const
//    {
//        return( 0 );
//    }
//    
//    Real diffx(Real x, Real y, Real t) const
//    {
//        return( 0 );
//    }
//    
//    Real diffy(Real x, Real y, Real t) const
//    {
//        return( 0 );
//    }
//    
//    Real difft(Real x, Real y, Real t) const
//    {
//        return( 0 );
//    }
//    
//    Real diff2x(Real x, Real y, Real t) const
//    {
//        return( 0 );
//    }
//    
//    Real diff2y(Real x, Real y, Real t) const
//    {
//        return( 0 );
//    }
//};

//class rappTest2DP : public basics::function2D {
//    public:  
//    Real val(Real x, Real y, Real t) const
//    {
//        return( sin(x)*cos(M_PI*y) );
//    }

//    Real diffx(Real x, Real y, Real t) const
//    {
//        return( cos(x)*cos(M_PI*y) );
//    }

//    Real diffy(Real x, Real y, Real t) const
//    {
//        return( -M_PI*sin(x)*sin(M_PI*y) );
//    }

//    Real difft(Real x, Real y, Real t) const
//    {
//        return 0;
//    }

//    Real diff2x(Real x, Real y, Real t) const
//    {
//        return( -sin(x)*cos(M_PI*y) );
//    }
//    
//    Real diff2y(Real x, Real y, Real t) const
//    {
//        return( -M_PI*M_PI*sin(x)*cos(M_PI*y) );
//    }
//};

class channelTest2DX : public basics::function2D {
  public:
    channelTest2DX()
    {
      m_id = 5;
    }

    Real val(Real x, Real y, Real t) const
    {
      return( (1.f-y*y)*sin(t) );
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( -2.f*y*sin(t) );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return( (1.f-y*y)*cos(t) );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return -2*sin(t);
    }
};

class channelTest2DY : public basics::function2D {
  public:
    channelTest2DY()
    {
      m_id = -5;
    }

    Real val(Real x, Real y, Real t) const
    {
      return( 0 );
    }
};


class channelTest2DP : public basics::function2D {
  public:
    channelTest2DP(Real n_nu) :
      nu(n_nu)
  {
    m_id = 5;
  }

    Real val(Real x, Real y, Real t) const
    {
      return 2*nu*y;
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( 2*nu );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return( -2*nu*y );
    }
  protected:
    Real nu;
};

#define AMPLITUDE 1.e-8

class channelPerturbedTest2DX : public basics::function2D {
  public:
    Real val(Real x, Real y, Real t) const
    {
      //sin(x-t)*sin(x-t)*sin(Real(2)*M_PI*y) // 1.avi
      //sin(x-t)*cos(Real(2)*M_PI*y)*t*t*t // 2.avi
      //sin(Real(2)*M_PI*y)*sin(x)*sin(x)*sin(t) // 3.avi
      //sin(x)*sin(x)*sin(t) // 4.avi
      //sin(Real(2)*M_PI*y)*sin(t) // 5.avi
      //sin(x)*sin)(x)
      return( 1-y*y+AMPLITUDE*sin(x-t)*sin(x-t)*sin(Real(2)*M_PI*y) );
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return( 0 );
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( -2.f*y );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return( 0.f ); 
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return( 0 );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return 0.f;
    }

    static int id() 
    {
      return( 2 );
    }
};

class channelPerturbedTest2DY : public basics::function2D {
  public:
    Real val(Real x, Real y, Real t) const
    {
      return( AMPLITUDE*(-sin(Real(2)*(x-t))*sin(M_PI*y)*sin(M_PI*y)/M_PI) );
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return( 0.f );
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( 0.f );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return( 0.f ); 
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return( 0.f );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return( 0.f );
    }

    static int id() 
    {
      return( 2 );
    }
};

class extremelySimpleTestX : public basics::function2D {
  public:
    Real val(Real x, Real y, Real t) const
    {
      return( 1.f-y*y );
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return( 0 );
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( -2*y );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return( 0 );
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return( 0 );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return( -2 );
    }
};

class extremelySimpleTestY : public basics::function2D {
  public:
    Real val(Real x, Real y, Real t) const
    {
      return( 0 );
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return( 0 );
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( 0 );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return( 0 );
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return( 0 );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return( 0 );
    }
};

//class extremelySimpleTestP : public basics::function2D {
//    public:
//        Real val(Real x, Real y, Real t) const
//        {
//            return( cos(x) );
//        }
//        
//        Real diffx(Real x, Real y, Real t) const
//        {
//            return( -sin(x) );
//        }
//        
//        Real diffy(Real x, Real y, Real t) const
//        {
//            return( 0 );
//        }
//        
//        Real difft(Real x, Real y, Real t) const
//        {
//            return( 0 );
//        }
//        
//        Real diff2x(Real x, Real y, Real t) const
//        {
//            return( -cos(x) );
//        }
//        
//        Real diff2y(Real x, Real y, Real t) const
//        {
//            return( 0 );
//        }
//};
class extremelySimpleTestP : public basics::function2D {
  public:
    Real val(Real x, Real y, Real t) const
    {
      return( sin(x)*cos(M_PI*y) );
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return( cos(x)*cos(M_PI*y) );
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( -sin(x)*M_PI*sin(M_PI*y) );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return( 0 );
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return( -sin(x)*cos(M_PI*y) );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return( -sin(x)*M_PI*M_PI*cos(M_PI*y) );
    }
};

class articleTest1 : public basics::function2D {
  public:
    Real val(Real x, Real y, Real t) const
    {
      return( sin((Real(7)/8)*M_PI*(pow(x-2,2)+pow(y-2,2)+(Real(1)/2)*M_PI))*sin(t) );
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return( (Real(7)/8)*cos((Real(7)/8)*M_PI*(pow(x-2,2)+pow(y-2,2))+
            (Real(1)/2)*M_PI)*M_PI*(2*x-4)*sin(t) );
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return( (Real(7)/8)*cos((Real(7)/8)*M_PI*(pow(x-2,2)+pow(y-2,2))
            +(Real(1)/2)*M_PI)*M_PI*(2*y-4)*sin(t) );
    }

    Real difft(Real x, Real y, Real t) const
    {
      return( sin((Real(7)/8)*M_PI*(pow(x-2,2)+pow(y-2,2)+(Real(1)/2)*M_PI))*cos(t) );
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return( -(Real(49)/64)*sin((Real(7)/8)*M_PI*(pow(x-2,2)+pow(y-2,2))+
            (Real(1)/2)*M_PI)*pow(M_PI,2)*pow(2*x-4,2)*sin(t)+
          (Real(7)/4)*cos((Real(7)/8)*M_PI*(pow(x-2,2)+pow(y-2,2))+
            (Real(1)/2)*M_PI)*M_PI*sin(t) );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return( -(Real(49)/64)*sin((Real(7)/8)*M_PI*(pow(x-2,2)+pow(y-2,2))+
            (Real(1)/2)*M_PI)*pow(M_PI,2)*pow(2*y-4,2)*sin(t)+
          (Real(7)/4)*cos((Real(7)/8)*M_PI*(pow(x-2,2)+pow(y-2,2))+
            (Real(1)/2)*M_PI)*M_PI*sin(t) );
    }
};

class articleTest1X : public basics::function2D {
  public:
    Real val(Real x, Real y, Real t) const
    {
      if( x == 0 && y == 0 )
        return 0;

      Real r            = sqrt(x*x+y*y);
      Real r2  = x*x+y*y;
      Real cospr= cos(M_PI*r);
      Real sinpr= sin(M_PI*r);
      Real st  = sin(t);

      return ((Real(1)/5)*sinpr*sinpr*y*x/r2 -
          (Real(1)/5)*sinpr*(2*M_PI*r*cospr+sinpr)*x*y/r2)*st;
    }

    Real diffx(Real x, Real y, Real t) const
    {
      Real r= sqrt(x*x+y*y);
      Real r2= x*x+y*y;
      Real cospr= cos(M_PI*r);
      Real cospr2= cospr*cospr;
      Real sinpr= sin(M_PI*r);
      Real x2= x*x;
      Real x4= pow(x,4);
      Real y2= y*y;
      Real st= sin(t);

      return (-(Real(2)/5)*(-x2*M_PI*r+2*x2*cospr2*M_PI*r+
            sinpr*cospr*y*y)*M_PI*y/pow(r2,Real(3)/2))*st;
    }

    Real diffy(Real x, Real y, Real t) const
    {
      Real r= sqrt(x*x+y*y);
      Real r2= x*x+y*y;
      Real cospr= cos(M_PI*r);
      Real cospr2= cospr*cospr;
      Real sinpr= sin(M_PI*r);
      Real x2= x*x;
      Real x4= pow(x,4);
      Real y2= y*y;
      Real st= sin(t);

      return (-(Real(2)/5)*(sinpr*cospr*x2-M_PI*y2*r+2*cospr2*M_PI*y2*r)*
          M_PI*x/pow(r2,Real(3)/2))*st;
    }

    Real difft(Real x, Real y, Real t) const
    {
      if( x == 0 && y == 0 )
        return 0;

      Real r= sqrt(x*x+y*y);
      Real r2= x*x+y*y;
      Real cospr= cos(M_PI*r);
      Real sinpr= sin(M_PI*r);
      Real ct= cos(t);

      return ((Real(1)/5)*sinpr*sinpr*y*x/r2 -
          (Real(1)/5)*sinpr*(2*M_PI*r*cospr+sinpr)*x*y/r2)*ct;
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      Real r= sqrt(x*x+y*y);
      Real r2= x*x+y*y;
      Real cospr= cos(M_PI*r);
      Real cospr2= cospr*cospr;
      Real sinpr= sin(M_PI*r);
      Real x2= x*x;
      Real x4= pow(x,4);
      Real y2= y*y;
      Real p2    = M_PI*M_PI;
      Real st    = sin(t);

      return (-(Real(2)/5)*(-4*cospr*p2*sinpr*x4-
            4*cospr*p2*sinpr*y2*x2-3*M_PI*y2*r+
            6*cospr2*M_PI*y2*r-3*sinpr*cospr*y2)*M_PI*x*y/pow(r2,Real(5)/2))*st;
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      Real r    = sqrt(x*x+y*y);
      Real r2    = x*x+y*y;
      Real cospr  = cos(M_PI*r);
      Real cospr2  = cospr*cospr;
      Real sinpr  = sin(M_PI*r);
      Real x2    = x*x;
      Real y2    = y*y;
      Real y4    = pow(y,4);
      Real p2    = M_PI*M_PI;
      Real st    = sin(t);

      return (-(Real(2)/5)*(-4*cospr*p2*sinpr*y2*x2+
            6*x2*M_PI*r*cospr2-3*sinpr*cospr*x2-
            3*x2*M_PI*r-4*cospr*p2*y4*sinpr)*M_PI*y*x/pow(r2,Real(5)/2))*st;
    }
};

class articleTest1Y : public basics::function2D {
  public:
    Real val(Real x, Real y, Real t) const
    {
      if( x == 0 && y == 0 )
        return 0;

      Real r    = sqrt(x*x+y*y);
      Real r2    = x*x+y*y;
      Real cospr  = cos(M_PI*r);
      Real sinpr  = sin(M_PI*r);
      Real x2    = x*x;
      Real y2    = y*y;
      Real st    = sin(t);

      return (Real(1)/5*sinpr*sinpr*y2/r2+
          Real(1)/5*sinpr*(2*M_PI*r*cospr+sinpr)*x2/r2)*st;
    }

    Real diffx(Real x, Real y, Real t) const
    {
      Real r    = sqrt(x*x+y*y);
      Real r2    = x*x+y*y;
      Real cospr  = cos(M_PI*r);
      Real cospr2  = cospr*cospr;
      Real sinpr  = sin(M_PI*r);
      Real x2    = x*x;
      Real x4    = pow(x,4);
      Real y2    = y*y;
      Real st    = sin(t);

      return ((Real(2)/5)*(2*x2*M_PI*r*cospr2+2*sinpr*cospr*x2-
            x2*M_PI*r+3*sinpr*y2*cospr)*M_PI*x/pow(r2,Real(3)/2))*st;
    }

    Real diffy(Real x, Real y, Real t) const
    {
      Real r    = sqrt(x*x+y*y);
      Real r2    = x*x+y*y;
      Real cospr  = cos(M_PI*r);
      Real cospr2  = cospr*cospr;
      Real sinpr  = sin(M_PI*r);
      Real x2    = x*x;
      Real y2    = y*y;
      Real st    = sin(t);

      return ((Real(2)/5)*(-x2*M_PI*r+2*x2*M_PI*r*cospr2+
            sinpr*y2*cospr)*M_PI*y/pow(r2,Real(3)/2))*st;
    }

    Real difft(Real x, Real y, Real t) const
    {
      if( x == 0 && y == 0 )
        return 0;

      Real r    = sqrt(x*x+y*y);
      Real r2    = x*x+y*y;
      Real cospr  = cos(M_PI*r);
      Real sinpr  = sin(M_PI*r);
      Real x2    = x*x;
      Real y2    = y*y;
      Real ct    = cos(t);

      return( (Real(1)/5*sinpr*sinpr*y2/r2+
            Real(1)/5*sinpr*(2*M_PI*r*cospr+sinpr)*x2/r2)*ct );
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      Real r    = sqrt(x*x+y*y);
      Real r2    = x*x+y*y;
      Real cospr  = cos(M_PI*r);
      Real cospr2  = cospr*cospr;
      Real sinpr  = sin(M_PI*r);
      Real x2    = x*x;
      Real x4    = pow(x,4);
      Real y2    = y*y;
      Real y4    = pow(y,4);
      Real p2    = M_PI*M_PI;
      Real st    = sin(t);

      return (-(Real(2)/5)*(4*sinpr*p2*cospr*pow(x,6)-6*x4*cospr2*M_PI*r+
            4*x4*sinpr*p2*cospr*y2+3*x4*M_PI*r+
            6*x2*y2*M_PI*r-12*x2*cospr2*M_PI*r*y2-
            3*sinpr*y4*cospr)*M_PI/pow(r2,Real(5)/2))*st;
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      Real r    = sqrt(x*x+y*y);
      Real r2    = x*x+y*y;
      Real cospr  = cos(M_PI*r);
      Real cospr2 = cospr*cospr;
      Real sinpr  = sin(M_PI*r);
      Real x2     = x*x;
      Real x4     = pow(x,4);
      Real y2     = y*y;
      Real y4    = pow(y,4);
      Real p2    = M_PI*M_PI;
      Real st    = sin(t);

      return (-(Real(2)/5)*(4*x4*sinpr*p2*cospr*y2-
            2*x4*cospr2*M_PI*r+x4*M_PI*r-
            3*x2*sinpr*cospr*y2+4*x2*sinpr*p2*cospr*y4+
            2*x2*cospr2*M_PI*r*y2-x2*y2*M_PI*r-2*cospr2*M_PI*r*y4+
            M_PI*r*y4)*M_PI/pow(r2,Real(5)/2))*st;
    }
};

class articleTest1P : public basics::function2D {
  public:
    Real val(Real x, Real y, Real t) const
    {
      Real r      = sqrt(x*x+y*y);
      Real sinpr  = sin(M_PI*r);
      Real st     = sin(t);

      return sinpr*sinpr*st;
    }

    Real diffx(Real x, Real y, Real t) const
    {
      Real r      = sqrt(x*x+y*y);
      Real cospr  = cos(M_PI*r);
      Real sinpr  = sin(M_PI*r);
      Real st     = sin(t);

      return 2*cospr*sinpr*M_PI*x/r*st;
    }

    Real diffy(Real x, Real y, Real t) const
    {
      Real r      = sqrt(x*x+y*y);
      Real cospr  = cos(M_PI*r);
      Real sinpr  = sin(M_PI*r);
      Real st     = sin(t);

      return 2*cospr*sinpr*M_PI*y/r*st;
    }

    Real difft(Real x, Real y, Real t) const
    {
      Real r      = sqrt(x*x+y*y);
      Real sinpr  = sin(M_PI*r);
      Real ct     = cos(t);

      return sinpr*sinpr*ct;
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      Real r    = sqrt(x*x+y*y);
      Real r2    = x*x+y*y;
      Real cospr  = cos(M_PI*r);
      Real cospr2  = cospr*cospr;
      Real sinpr  = sin(M_PI*r);
      Real x2    = x*x;
      Real x4    = pow(x,4);
      Real y2    = y*y;
      Real y4    = pow(y,4);
      Real p2    = M_PI*M_PI;
      Real st    = sin(t);

      return (2*(2*cospr2*M_PI*x2*r-M_PI*x2*r+cospr*sinpr*y2)*M_PI*st)/pow(r2,Real(3)/2);
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      Real r    = sqrt(x*x+y*y);
      Real r2    = x*x+y*y;
      Real cospr  = cos(M_PI*r);
      Real cospr2  = cospr*cospr;
      Real sinpr  = sin(M_PI*r);
      Real x2    = x*x;
      Real x4    = pow(x,4);
      Real y2    = y*y;
      Real y4    = pow(y,4);
      Real p2    = M_PI*M_PI;
      Real st    = sin(t);

      return ((2*(cospr*sinpr*x2-M_PI*y2*r+2*cospr2*M_PI*y2*r))*M_PI*st)/pow(r2,Real(3)/2);
    }
};

#endif

