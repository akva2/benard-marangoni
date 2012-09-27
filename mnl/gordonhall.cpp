#include "gordonhall.h"
#include "gll.h"

namespace mnl {
  namespace utilities {
    GordonHall::GordonHall(int N, int M) :
      m_bottom(NULL), m_right(NULL), m_top(NULL), m_left(NULL),
      m_J("jacobian",N+1,M+1), m_Dx("geometry derivative X",N+1,N+1), 
      m_Dy("geometry derivative Y",M+1,M+1), m_map("mapping",N+1,M+1),
      m_D("geometry derivatives")
    {
      m_Dx = utilities::GLL::LagrangeDerivativeMatrix(N+1);
      m_Dy = utilities::GLL::LagrangeDerivativeMatrix(M+1);
      m_D.add(m_Dx,4);
    }

    GordonHall::GordonHall(int N, int M, basics::function1Dto2D& bottom,
        basics::function1Dto2D& right, basics::function1Dto2D& top, basics::function1Dto2D& left) : 
      m_bottom(&bottom), m_right(&right), m_top(&top), m_left(&left),
      m_J("jacobian",N+1,M+1), m_Dx("geometry derivative X",N+1,N+1), 
      m_Dy("geometry derivative Y",M+1,M+1), m_map("mapping",N+1,M+1),
      m_D("geometry derivatives")
    {
      m_Dx = utilities::GLL::LagrangeDerivativeMatrix(N+1);
      m_Dy = utilities::GLL::LagrangeDerivativeMatrix(M+1);
      m_D.add(m_Dx,4);
    }

    void GordonHall::compute(Real t)
    {
      compute(m_map,*m_bottom,*m_right,*m_top,*m_left,t);
      computeJacobian(m_J,m_map);
    }

    void GordonHall::compute(basics::Field2D& map, const basics::function1Dto2D& g1, const basics::function1Dto2D& g2, const basics::function1Dto2D& g3, const basics::function1Dto2D& g4, Real t)
    {
      basics::Vector grid("GLL grid xi",map.rows());
      GLL::GaussLobattoLegendreGrid(grid);

      int N = map.rows();

      /* distribute points along the boundaries */
      for( int i=0;i<N;++i ) { // bottom
        /* bottom */
        std::pair<Real,Real> val = g1(grid[i],t);
        map.X()[0][i] = val.first;
        map.Y()[0][i] = val.second;

        /* right */
        val = g2(grid[i],t);
        map.X()[i][N-1] = val.first;
        map.Y()[i][N-1] = val.second;

        /* top */
        val = g3(grid[i],t);
        map.X()[N-1][i] = val.first;
        map.Y()[N-1][i] = val.second;

        /* left */
        val = g4(grid[i],t);
        map.X()[i][0] = val.first;
        map.Y()[i][0] = val.second;
      }

      compute(map.X(),grid);
      compute(map.Y(),grid);
    }

    void GordonHall::compute(basics::Matrix& x, const basics::Vector& grid)
    {
      basics::Matrix* xa = utilities::g_manager.clone(x);
      basics::Matrix* xb = utilities::g_manager.clone(x);
      basics::Matrix* xc = utilities::g_manager.clone(x);

      int N=x.rows()-1;

      (*xa) = x;
      (*xb) = x;

      /* prefers memory over speed */
#pragma omp parallel for schedule(static)
      for( int n=0;n<N+1;++n )
        for( int m=1;m<N;++m )
          (*xa)[n][m] = (Real(1)-grid[m])/Real(2)*x[n][0] + 
            (Real(1)+grid[m])/Real(2)*x[n][N];

#pragma omp parallel for schedule(static)
      for( int n=1;n<N;++n )
        for( int m=0;m<N+1;++m )
          (*xb)[n][m] = (Real(1)-grid[n])/Real(2)*x[0][m] + 
            (Real(1)+grid[n])/Real(2)*x[N][m];

#pragma omp parallel for schedule(static)
      for( int n=0;n<N+1;++n )
        for( int m=0;m<N+1;++m )
          (*xc)[n][m] = (Real(1)-grid[m])/Real(2)*(Real(1)-grid[n])/Real(2)*x[0][0] +
            (Real(1)+grid[m])/Real(2)*(Real(1)-grid[n])/Real(2)*x[0][N] +
            (Real(1)-grid[m])/Real(2)*(Real(1)+grid[n])/Real(2)*x[N][0] +
            (Real(1)+grid[m])/Real(2)*(Real(1)+grid[n])/Real(2)*x[N][N];

      x = *xa;	
      x += *xb;
      x -= *xc;

      utilities::g_manager.unlock(xa);
      utilities::g_manager.unlock(xb);
      utilities::g_manager.unlock(xc);
    }

    std::pair<Real,Real> GordonHall::evaluate(Real x, Real y, const basics::Matrix& X,
        const basics::Matrix& Y, const basics::Vector& grid)
    {
      std::pair<Real,Real> result(0,0);
      for (int m=0;m<X.rows();++m)
        for (int n=0;n<X.cols();++n) {
          Real l = GLL::Lagrange(x,n,grid)*GLL::Lagrange(y,m,grid);
          result.first += X[m][n]*l;
          result.second += Y[m][n]*l;
        }

      return( result );
    }

    void GordonHall::computeJacobian(basics::Matrix& J, const basics::Field2D& map)
    {
      basics::multTranspose(m_D[0],m_Dx,map.X(),'N','N');
      basics::multTranspose(m_D[1],map.X(),m_Dy,'N','T');
      basics::multTranspose(m_D[2],m_Dx,map.Y(),'N','N');
      basics::multTranspose(m_D[3],map.Y(),m_Dy,'N','T');

#pragma omp parallel for schedule(static)
      for( int beta=0;beta<map.cols();++beta) {
        for( int alpha=0;alpha<map.rows();++alpha ) {
          Real xxi  = m_D[0][beta][alpha];
          Real xeta = m_D[1][beta][alpha];
          Real yxi  = m_D[2][beta][alpha];
          Real yeta = m_D[3][beta][alpha];
          J[beta][alpha] = xxi*yeta-xeta*yxi;
        }
      }
    }

    GordonHall3D::GordonHall3D(int N, int M, int P) :
      m_J("jacobian",N+1,M+1,P+1), m_Dx("geometry derivative X",N+1,N+1),
      m_map("mapping",N+1,M+1,P+1), m_D("geometry derivatives")
    {
      m_Dx = utilities::GLL::LagrangeDerivativeMatrix(N+1);
      m_D.add(m_J,9);
    }

    void GordonHall3D::computeJacobian(basics::Matrices& J, const basics::Field3D& map)
    {
#pragma omp parallel for schedule(static)
      for( int l=0;l<map.X().matrices();++l ) {
        basics::multTranspose(m_D[0][l],m_Dx,map.X()[l],'N','N');
        basics::multTranspose(m_D[1][l],map.X()[l],m_Dx,'N','T');
        basics::multTranspose(m_D[3][l],m_Dx,map.Y()[l],'N','N');
        basics::multTranspose(m_D[4][l],map.Y()[l],m_Dx,'N','T');
        basics::multTranspose(m_D[6][l],m_Dx,map.Z()[l],'N','N');
        basics::multTranspose(m_D[7][l],map.Z()[l],m_Dx,'N','T');
      }
      applyLocalGlobal(m_D[2],map.X(),m_Dx,'N','T');
      applyLocalGlobal(m_D[5],map.Y(),m_Dx,'N','T');
      applyLocalGlobal(m_D[8],map.Z(),m_Dx,'N','T');

#pragma omp parallel for schedule(static)
      for( int gamma=0;gamma<map.matrices();++gamma ) {
        for( int beta=0;beta<map.cols();++beta) {
          for( int alpha=0;alpha<map.rows();++alpha ) {
            Real xxi    = m_D[0][gamma][beta][alpha];
            Real xeta   = m_D[1][gamma][beta][alpha];
            Real xgamma = m_D[2][gamma][beta][alpha];
            Real yxi    = m_D[3][gamma][beta][alpha];
            Real yeta   = m_D[4][gamma][beta][alpha];
            Real ygamma = m_D[5][gamma][beta][alpha];
            Real zxi    = m_D[6][gamma][beta][alpha];
            Real zeta   = m_D[7][gamma][beta][alpha];
            Real zgamma = m_D[8][gamma][beta][alpha];
            J[gamma][beta][alpha] =   xxi*(yeta*zgamma-zeta*ygamma)
              -xeta*(yxi*zgamma-zxi*ygamma)
              +xgamma*(yxi*zeta-zxi*yeta);
          }
        }
      }
    }
  }
}

