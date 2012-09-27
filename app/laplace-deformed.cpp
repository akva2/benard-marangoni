#include "mnl/vector.h"
#include "mnl/gll.h"
#include "mnl/matrix.h"
#include "mnl/gordonhall.h"
#include "mnl/function.h"
#include "mnl/cgsolver.h"
#include "mnl/field.h"
#include "mnl/timer.h"

#include "bm/poissonsolver-lgw.h"
#include "bm/sem.h"

#ifdef HAS_MPI
#include <mpi.h>
#endif

using namespace mnl;
using namespace std;
using namespace legendreLegendreW;

#include "bm/functions3d.h"
#include "bm/bigcircle.h"

#include "bm/laplacian.h"

class source3 : public basics::function3D {
  public:
    source3(basics::function3D& n_X,Real n_nu, SEMLaplacianEvaluator::BC bc) : 
      X(n_X), nu(n_nu), m_bc(bc)
  {
  }

    Real val(Real x, Real y, Real z, Real t) const
    {
      return 1;

      if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
        if( m_bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
          return 0;
        else
          return -3*M_PI*M_PI*cos(M_PI*z)+2*M_PI*M_PI;

      return( -(X.diff2x(x,y,z,t)+X.diff2y(x,y,z,t)+X.diff2z(x,y,z,t)) + nu*X.val(x,y,z,t) );
    }
  private:
    basics::function3D& X;
    Real nu;
    SEMLaplacianEvaluator::BC m_bc;
};

#undef RADIUS
#define RADIUS 2

class testTapered : public basics::function3D {
  public:
    Real val(Real x, Real y, Real z, Real t) const
    {
      Real Rz = RADIUS*(z/2+Real(3)/4);
      Real r2 = x*x+y*y;

      return (Rz*Rz-r2)*sin(M_PI*z);
    }

    Real diff2x(Real x, Real y, Real z, Real t) const
    {
      return -2*sin(M_PI*z);
    }

    Real diff2y(Real x, Real y, Real z, Real t) const
    {
      return -2*sin(M_PI*z);
    }

    Real diff2z(Real x, Real y, Real z, Real t) const
    {
      Real Rz = RADIUS*(z/2+Real(3)/4);
      Real r2 = x*x+y*y;

      return -M_PI*M_PI*(Rz*Rz-r2)*sin(M_PI*z);
    }
};

class articleTest1XV : public basics::function3D {
  public:
    Real val(Real x, Real y, Real z, Real t) const
    {
      SETUPLOCALS

        Real z2 = 1-cos(M_PI*z);

      return( cospr*z2 );
    }

    Real diff2x(Real x, Real y, Real z, Real t) const
    {
      SETUPLOCALS

        Real z2 = 1-cos(M_PI*z);

      return z2*(-cospr*pow(M_PI*x,2)/r2+sinpr*M_PI*pow(x,2)/pow(r,3)-sinpr*M_PI/r);
    }

    Real diff2y(Real x, Real y, Real z, Real t) const
    {
      SETUPLOCALS

        Real z2 = 1-cos(M_PI*z);

      return z2*(-cospr*pow(M_PI*y,2)/r2+sinpr*M_PI*pow(y,2)/pow(r,3)-sinpr*M_PI/r);
    }

    Real diff2z(Real x, Real y, Real z, Real t) const
    {
      SETUPLOCALS

        cospz = cos(M_PI*z);
      return cospr*cospz*M_PI*M_PI;
    }
};

int main(int argc, char** argv)
{
  int rank=0, size=1;

#ifdef HAS_MPI
  int aquired;
  MPI_Init_thread(&argc, &argv,MPI_THREAD_SERIALIZED,&aquired);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  cout << "inited mpi in single threaded mode" << endl;
#endif

  int N;
  if( argc < 2 )
    N = 13;
  else
    N=atoi(argv[1]);
  bool preconditioned=true;
  if( argc > 2 )
    preconditioned=(atoi(argv[2])==1?false:true);
  int d3D=0;
  if( argc > 3 )
    d3D=(atoi(argv[3]));
  bool neumann=false;
  if( argc > 4 )
    neumann=(atoi(argv[4])==1?true:false);

  basics::Vector grid("GLL grid",N+1);
  basics::Vector weight("GLL weight",N+1);
  utilities::GLL::GaussLobattoLegendreGrid(grid);
  utilities::GLL::GaussLobattoLegendreWeights(weight,grid);

  Real Lz = 1;

  bigCircleGeometry circle(N,N,grid,weight);
  bigCircleGeometry3D circle3D(weight,circle);
  circle.setLz(Lz);
  Real nu = 0;

  SEMLaplacianEvaluator::BC bc2d = SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET;
  SEMLaplacianEvaluator::BC bc3d = SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET;
  poissonSolver::BC bcp = poissonSolver::Homogenous;
  articleTest1 article;
  articleTest1XV articleV;

  basics::function3D* fun = &article;
  if( neumann ) {
    cout << "using neumann bc's" << endl;
    bc2d = SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN;
    bc3d = SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM;
    bcp = poissonSolver::HomogenousLeft;
    fun = &articleV;
  }
  poissonSolver  SP(N,4,0,bcp,true,Lz,-1,-1,circle3D.m_level);
  poissonSolver SP2(N,4,0,poissonSolver::Nonhomogenous,true,Lz,-1,-1,circle3D.m_level);

  //    testTapered article;
  source3 source(*fun,nu,bc3d);

  //    testTrivial mixed;
  //    solTrivial msol;

  SEMLaplacianEvaluator eval(circle,SP.Dx(),SP.weightX(),NULL,0,1,bc2d);
  extrudedLaplacianEvaluator eval25(eval,SP,nu,&SP2.Ax(),true,rank,size,bc3d);
  extrudedLaplacianPreconditioner eval3(eval25,circle3D);

  SEM3DLaplacianEvaluator eval2(circle3D,SP.weightX(),SP.Dx(),nu,rank,size,bc3d);

  SEM3DFEMLaplacianPreconditioner pre3D(circle3D,SP.weightX(),grid,nu,rank,size);
  pre3D.m_bc = eval2.m_bc;

  basics::coarseGrid desc = circle3D.getDivisionInfo(size);
  basics::matricesStack u("velocity"  ,N+1,N+1,(N+1),desc[rank].elements.size());
  basics::matricesStack u2("velocity 2",N+1,N+1,(N+1),desc[rank].elements.size());
  basics::matricesStack exact("velocity 2",N+1,N+1,(N+1),desc[rank].elements.size());

  circle3D.evaluate(u,source,M_PI/2,desc[rank].elements);
  circle3D.mass(u,desc[rank].elements);
  mask(u,circle3D,rank,size,eval2.m_bc);
  exact = u;
  int whole  = utilities::g_profiler.add("whole");
  const int solves=10;
  for( int i=0;i<solves;++i ) {
    u = exact;
    u2 = exact;
    dssum(u,circle3D,rank,size);
    utilities::g_profiler.start(whole);
    eval2.m_nosum = eval25.m_nosum = pre3D.m_nosum = &u2;
    MPIGeometryDotter<basics::geometryStack3D> dotter(circle3D,desc[rank].elements,rank,size);
    if( preconditioned )
      if( d3D == 1)  {
        cout << "iters " << utilities::CGSolver::solve(u,eval2,pre3D,dotter,1.e-6) << endl;
      } else if( d3D == 2 ) {
        basics::matricesStack* temp = circle3D.localToGlobalZ(u,rank,size);
        basics::matricesStack* temp2 = circle3D.localToGlobalZ(u2,rank,size,false,true);
        eval25.solve(*temp,*temp2);
        circle3D.globalToLocalZ(u,temp,rank,size);
        utilities::g_manager.unlock(temp2);
      }
      else {
        cout << "iters " << utilities::CGSolver::solve(u,eval2,eval3,dotter,1.e-6) << endl;
      }
    else {
      if( d3D == 2 ) {
        exact = u;
        basics::matricesStack* temp  = circle3D.localToGlobalZ(u,rank,size);
        basics::matricesStack* temp2 = utilities::g_manager.clone(*temp);
        eval25.solve(*temp,*temp2);
        circle3D.globalToLocalZ(u,temp,rank,size);
        utilities::g_manager.unlock(temp2);
      } else
        cout << "iters " << utilities::CGSolver::solve(u,eval2,dotter,1.e-6) << endl;
    }

    utilities::g_profiler.pause(whole);
  }

  string report = errorReport(u,circle3D,*fun,grid,M_PI/2,&eval2);
  if( rank == 0 ) {
    cout << report << endl;
    cout << utilities::g_profiler.report();
    cout << "we ran " << solves << " times" << endl;
  }

#ifdef HAS_MPI
  MPI_Finalize();
#endif

  return 0; // AY-OH-KAY
}

