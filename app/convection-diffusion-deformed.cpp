#include "mnl/vector.h"
#include "mnl/gll.h"
#include "mnl/matrix.h"
#include "mnl/gordonhall.h"
#include "mnl/function.h"
#include "mnl/cgsolver.h"
#include "mnl/field.h"
#include "mnl/timer.h"

#include "bm/poissonsolver-lgw.h"

using namespace mnl;
using namespace std;
using namespace legendreLegendreW;

#include "bm/functions3d.h"
#include "bm/bigcircle.h"
#include "bm/laplacian.h"
#include "bm/sem.h"
#include "bm/convection.h"

#ifdef HAS_MPI
#include <mpi.h>
#endif

class source3 : public basics::function3D {
  public:
    source3(basics::function3D& n_X, basics::function3D& n_Y,
        basics::function3D& n_Z, basics::function3D& n_phi) : 
      X(n_X), Y(n_Y), Z(n_Z), phi(n_phi) 
  {
  }

    Real val(Real x, Real y, Real z, Real t) const
    {
      //            if( abs(x) < 1.e-14 && abs(y) < 1.e-14)
      //                return 0;

      return(  phi.difft(x,y,z,t) - (phi.diff2x(x,y,z,t)+phi.diff2y(x,y,z,t)+phi.diff2z(x,y,z,t)) 
          +X.val(x,y,z,t)*phi.diffx(x,y,z,t)
          +Y.val(x,y,z,t)*phi.diffy(x,y,z,t)
          +Z.val(x,y,z,t)*phi.diffz(x,y,z,t) );
    }

    basics::function3D& X;
    basics::function3D& Y;
    basics::function3D& Z;
    basics::function3D& phi;
};

class articleTest1XD : public basics::function3D {
  public:
    articleTest1XD(Real R) : 
      m_R(R)
  {
  }

    Real val(Real x, Real y, Real z, Real t) const
    {
      return sin(t)*(pow(m_R,2)-x*x-y*y)*sin(M_PI*z);
    }

    Real difft(Real x, Real y, Real z, Real t) const
    {
      return cos(t)*(pow(m_R,2)-x*x-y*y)*sin(M_PI*z);
    }

    Real diffx(Real x, Real y, Real z, Real t) const
    {
      return -2*x*sin(t)*sin(M_PI*z);
    }

    Real diffy(Real x, Real y, Real z, Real t) const
    {
      return -2*y*sin(t)*sin(M_PI*z);
    }

    Real diffz(Real x, Real y, Real z, Real t) const
    {
      return sin(t)*(pow(m_R,2)-x*x-y*y)*M_PI*cos(M_PI*z);
    }

    Real diff2x(Real x, Real y, Real z, Real t) const
    {
      return -2*sin(t)*sin(M_PI*z);
    }

    Real diff2y(Real x, Real y, Real z, Real t) const
    {
      return -2*sin(t)*sin(M_PI*z);
    }

    Real diff2z(Real x, Real y, Real z, Real t) const
    {
      return -sin(t)*(pow(m_R,2)-x*x-y*y)*M_PI*M_PI*sin(M_PI*z);
    }

    Real m_R;
};

class articleTest1XV : public basics::function3D {
  public:
    articleTest1XV(Real R) : 
      m_R(R)
  {
  }

    Real val(Real x, Real y, Real z, Real t) const
    {
      return sin(t)*(pow(m_R,2)-x*x-y*y)*z*(2-z);
    }

    Real difft(Real x, Real y, Real z, Real t) const
    {
      return cos(t)*(pow(m_R,2)-x*x-y*y)*z*(2-z);
    }

    Real diffx(Real x, Real y, Real z, Real t) const
    {
      return -2*x*z*(2-z)*sin(t);
    }

    Real diffy(Real x, Real y, Real z, Real t) const
    {
      return -2*y*z*(2-z)*sin(t);
    }

    Real diffz(Real x, Real y, Real z, Real t) const
    {
      return -2*sin(t)*(pow(m_R,2)-x*x-y*y)*(-1+z);
    }

    Real diff2x(Real x, Real y, Real z, Real t) const
    {
      return -2*z*(2-z)*sin(t);
    }

    Real diff2y(Real x, Real y, Real z, Real t) const
    {
      return -2*z*(2-z)*sin(t);
    }

    Real diff2z(Real x, Real y, Real z, Real t) const
    {
      return -2*sin(t)*(pow(m_R,2)-x*x-y*y);
    }

    Real m_R;
};

//class articleTest1XN : public basics::function3D {
//public:
//    articleTest1XN(Real R) : 
//        m_R(R)
//    {
//    }

//    Real val(Real x, Real y, Real z, Real t) const
//    {
//        Real r = sqrt(x*x+y*y);

//        return sin(t)*cos(M_PI*r/m_R)*z*(2-z);
//    }

//    Real difft(Real x, Real y, Real z, Real t) const
//    {
//        Real r = sqrt(x*x+y*y);

//        return cos(t)*cos(M_PI*r/m_R)*z*(2-z);
//    }

//    Real diffx(Real x, Real y, Real z, Real t) const
//    {
//        Real r = sqrt(x*x+y*y);

//        return -2*sin(M_PI*r/m_R)*M_PI/m_R*x*z*(2-z);
//    }
//    
//    Real diffy(Real x, Real y, Real z, Real t) const
//    {
//        Real r = sqrt(x*x+y*y);

//        return -2*sin(M_PI*r/m_R)*M_PI/m_R*y*z*(2-z);
//    }

//    Real diffz(Real x, Real y, Real z, Real t) const
//    {
//        Real r = sqrt(x*x+y*y);

//        return -2*cos(M_PI*r/m_R)*(-1+z);
//    }
//    
//    Real diff2x(Real x, Real y, Real z, Real t) const
//    {
//        Real r = sqrt(x*x+y*y);

//        return -4*cos(M_PI*r/m_R)*pow(M_PI*x/m_R,2)*z*(2-z)-2*sin(M_PI*r/m_R)*M_PI*z*(2-z)/m_R;
//    }
//    
//    Real diff2y(Real x, Real y, Real z, Real t) const
//    {
//        Real r = sqrt(x*x+y*y);

//        return -4*cos(M_PI*r/m_R)*pow(M_PI*y/m_R,2)*z*(2-z)-2*sin(M_PI*r/m_R)*M_PI*z*(2-z)/m_R;
//    }

//    Real diff2z(Real x, Real y, Real z, Real t) const
//    {
//        Real r = sqrt(x*x+y*y);

//        return -2*cos(M_PI*r/m_R);
//    }

//    Real m_R;
//};

class articleTest1XN : public basics::function3D {
  public:
    articleTest1XN(Real R) : 
      m_R(R)
  {
  }

    Real val(Real x, Real y, Real z, Real t) const
    {
      return sin(t)*z*(2-z);
    }

    Real difft(Real x, Real y, Real z, Real t) const
    {
      return cos(t)*z*(2-z);
    }

    Real diffz(Real x, Real y, Real z, Real t) const
    {
      return (2-2*z)*sin(t);
    }

    Real diff2z(Real x, Real y, Real z, Real t) const
    {
      return -2*sin(t);
    }

    Real m_R;
};

int main(int argc, char** argv)
{
  int rank=0, size=1;

#ifdef HAS_MPI
  int aquired;
  MPI_Init_thread(&argc, &argv,MPI_THREAD_SERIALIZED,&aquired);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if (aquired != MPI_THREAD_SERIALIZED ) {
    cout << "threading unavailable (got " << aquired << "!) mpi/openmp combo will not work." << endl;
  } else
    cout << "init threaded mpi successfull " << MPI_THREAD_MULTIPLE << " " << aquired << endl;
#endif
  int N;
  if( argc < 2)
    N=12;
  else
    N = atoi(argv[1]);

  Real Dt = 1.e-2f;
  Real T = 1.f;
  if( argc > 3 )
    Dt = Real(atoi(argv[2]))/atoi(argv[3]);
  if( argc > 5 )
    T = Real(atoi(argv[4]))/atoi(argv[5]);
  bool preconditioned=true;
  if( argc > 6 )
    preconditioned=(atoi(argv[6])==1?false:true);

  Real Lz = 1;

  basics::Vector grid("GLL grid",N+1);
  basics::Vector weight("GLL grid",N+1);
  utilities::GLL::GaussLobattoLegendreGrid(grid);
  utilities::GLL::GaussLobattoLegendreWeights(weight,grid);

  bigCircleGeometry circle(N,N,grid,weight);
  bigCircleGeometry3D circle3D(weight,circle);
  circle.setLz(Lz);

  Real max = circle3D.max();
  articleTest1XD xD(max);
  articleTest1XD yD(max);
  articleTest1XD zD(max);
  articleTest1XD pD(max);
  articleTest1XV xV(max);
  articleTest1XV yV(max);
  articleTest1XV zV(max);
  articleTest1XV pV(max);
  articleTest1XN xN(max);
  articleTest1XN yN(max);
  articleTest1XN zN(max);
  articleTest1XN pN(max);
  source3* sourceX;

  SEMLaplacianEvaluator::BC bc3d = SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM;
  SEMLaplacianEvaluator::BC bc2d = SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET;
  poissonSolver::BC bcp=poissonSolver::Homogenous;
  if( bc3d == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET ) {
    sourceX = new source3(xD,yD,zD,pD);
  } else if( bc3d == SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM ) {
    bcp = poissonSolver::HomogenousLeft;
    bc2d = SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET;
    sourceX = new source3(xV,yV,zV,pV);
  } else if( bc3d == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM ) {
    bcp = poissonSolver::HomogenousLeft;
    bc2d = SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN;
    sourceX = new source3(xN,yN,zN,pN);
  }

  poissonSolver  SP(N,N,0,bcp,true,Lz,-1,-1,circle3D.m_level);
  poissonSolver SP2(N,N,0,poissonSolver::Nonhomogenous,true,Lz,-1,-1,circle3D.m_level);

  convectionEvaluator conv(circle,SP.Dx(),SP.weightX(),rank,size,&circle3D);
  conv.m_order = 1;
  conv.m_bc = bc3d;

  Real nu; 
  if( conv.m_order == 1 )
    nu = Real(3)/(2*Dt);
  else
    nu = Real(1)/Dt;
  SEMLaplacianEvaluator eval2(circle,SP.Dx(),SP.weightX(),NULL,0,1,bc2d);
  extrudedLaplacianEvaluator eval(eval2,SP,nu,&SP2.Ax(),preconditioned,rank,size,bc3d);
  extrudedLaplacianPreconditioner velocity3Dpre(eval,circle3D);
  SEM3DLaplacianEvaluator velocity3D(circle3D,SP.weightX(),SP.Dx(),nu,rank,size,bc3d);

  basics::coarseGrid desc = circle3D.getDivisionInfo(size);
  utilities::ringBuffer<basics::Field3<basics::matricesStack> > u;
  u.add(*utilities::g_manager.aquireMatricesStackField("velocity 1",N+1,N+1,N+1,desc[rank].elements.size()));
  u.add(*utilities::g_manager.aquireMatricesStackField("velocity 2",N+1,N+1,N+1,desc[rank].elements.size()));
  utilities::ringBuffer<basics::matricesStack> phi;
  phi.add(*utilities::g_manager.aquireMatricesStack("phi 1",N+1,N+1,N+1,desc[rank].elements.size()));
  phi.add(*utilities::g_manager.aquireMatricesStack("phi 2",N+1,N+1,N+1,desc[rank].elements.size()));
  phi.add(*utilities::g_manager.aquireMatricesStack("phi 3",N+1,N+1,N+1,desc[rank].elements.size()));
  basics::matricesStack* F   = utilities::g_manager.aquireMatricesStack("source",N+1,N+1,N+1,desc[rank].elements.size());

  /* set up initial conditions */
  Real t0=0;
  circle3D.evaluate(  u.get(0,-1),sourceX->X,sourceX->Y,sourceX->Z,t0-Dt,desc[rank].elements);
  circle3D.evaluate(phi.get(0,-1),sourceX->phi,t0-Dt,desc[rank].elements);
  circle3D.evaluate(phi.get(0, 0),sourceX->phi,t0,desc[rank].elements);

  int convection = utilities::g_profiler.add("convection");
  int eliptic = utilities::g_profiler.add("eliptic");

  int n;
  for( n=0;n<(T-t0)/Dt;++n ) {
    cout << "Time: " << t0+(n+1)*Dt << endl;

    /* set up references */
    basics::Field3<basics::matricesStack>& unm1 = u.get(n,-1);
    basics::Field3<basics::matricesStack>& un  	= u.get(n, 0);
    basics::matricesStack& pnm1 				= phi.get(n,-1);
    basics::matricesStack& pn  					= phi.get(n, 0);
    basics::matricesStack& pnp1 				= phi.get(n, 1);

    /* source */
    circle3D.evaluate(  un,sourceX->X,sourceX->Y,sourceX->Z,t0+n*Dt,desc[rank].elements);

    /* solve convection problems */
    utilities::g_profiler.start(convection);
    conv.solve(pnp1,pn,u,n,Dt,Dt,0);
    if( conv.m_order )
      conv.solve(*F,pnm1,u,n,Dt,Dt,-Dt);
    utilities::g_profiler.pause(convection);

    /* BDF2 */
    if( conv.m_order ) {
      pnp1 *= Real(2)/Dt;
      pnp1.axpy(Real(-1)/(Real(2)*Dt),*F);
    } else
      pnp1 *= Real(1)/Dt;

    circle3D.evaluate(*F,*sourceX,t0+(n+1)*Dt,desc[rank].elements);
    pnp1 += *F;

    circle3D.mass(pnp1,desc[rank].elements);

    /* delta value */
    velocity3D.evaluate(pnm1,pn,false,false);
    pnp1 -= pnm1;

    mask(pnp1,circle3D,rank,size,bc3d);
    *F = pnp1;
    dssum(pnp1,circle3D,rank,size);

    utilities::g_profiler.start(eliptic);
    cout << "start velocity solve" << endl;

    basics::matricesStack* temp = circle3D.localToGlobalZ(pnp1,rank,size);
    basics::matricesStack* temp2 = circle3D.localToGlobalZ(*F,rank,size,false,true);
    eval.solve(*temp,*temp2);
    circle3D.globalToLocalZ(pnp1,temp,rank,size);
    utilities::g_manager.unlock(temp2);

    //        cout << "iterations velocity " << utilities::CGSolver::solve(pnp1,velocity3D,velocity3Dpre,circle3D.getDotter(),1.e-10) << endl;
    utilities::g_profiler.pause(eliptic);
    pnp1 += pn;
  }

  eval.m_nu = 1;
  string errx = errorReport(phi.get(n),circle3D,sourceX->phi,grid,t0+n*Dt,&velocity3D);
  if( rank == 0 ) {
    cout << errx;
    cout << utilities::g_profiler.report();
  }

  return 0; // AY-OH-KAY
}

