#include "mnl/vector.h"
#include "mnl/gll.h"
#include "mnl/matrix.h"
#include "mnl/gordonhall.h"
#include "mnl/function.h"
#include "mnl/cgsolver.h"
#include "mnl/field.h"

using namespace mnl;
using namespace std;

#include "bm/bigcircle.h"
#include "bm/twoelement.h"
#include "bm/functions2d.h"
#include "bm/poissonsolver-lgw.h"
#include "bm/laplacian.h"
#include "bm/sem.h"

#ifdef HAS_MPI
#include <mpi.h>
#endif

class simplesin : public basics::function2D {
  public:
    Real val(Real x, Real y, Real t) const
    {
      return sin(M_PI*x)*sin(M_PI*y);
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return -sin(M_PI*x)*M_PI*M_PI*sin(M_PI*y);
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return -sin(M_PI*x)*M_PI*M_PI*sin(M_PI*y);
    }
};

class source2 : public basics::function2D {
  public:
    source2(basics::function2D& n_X,Real n_nu) : X(n_X), nu(n_nu) 
  {
  }

    Real val(Real x, Real y, Real t) const
    {
      if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
        return 0;

      return( -(X.diff2x(x,y,t)+X.diff2y(x,y,t))+nu*X.val(x,y,t) );
    }

    static int id() 
    {
      return( 2 );
    }

  private:
    basics::function2D& X;
    Real nu;
};

int main(int argc, char** argv)
{
  int rank=0, size=1;

#ifdef HAS_MPI
  int aquired;
  MPI_Init_thread(&argc, &argv,MPI_THREAD_MULTIPLE,&aquired);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if (aquired != MPI_THREAD_MULTIPLE )
    cout << "threading unavailable! mpi/openmp combo will not work." << endl;
#endif

  int N;
  if( argc < 2)
    N=12;
  else
    N = atoi(argv[1]);
  bool preconditioned = true;
  if( argc > 2 )
    preconditioned=(atoi(argv[2])==1?false:true);

  legendreLegendreW::poissonSolver SP(N,4,0,legendreLegendreW::poissonSolver::Homogenous,true,1);
  bigCircleGeometry circle(N,N,SP.gridX(),SP.weightX());

  SEMLaplacianEvaluator eval(circle,SP.Dx(),SP.weightX(),NULL,rank,size,SEMLaplacianEvaluator::PERIODIC);
  SEMFEMLaplacianPreconditioner pre(circle,SP.gridX(),SP.weightX(),0,rank,size);
  eval.m_nu = 0;

  basics::matrixStack& u     = *utilities::g_manager.aquireMatrixStack("laplacian",N+1,N+1,circle.size()/size);
  basics::matrixStack& u2    = *utilities::g_manager.aquireMatrixStack("buffer",N+1,N+1,circle.size()/size);
  basics::matrixStack& exact = *utilities::g_manager.aquireMatrixStack("buffer",N+1,N+1,circle.size()/size);

  simplesin article;
  source2 source(article,eval.m_nu);
  basics::coarseGrid desc = circle.getDivisionInfo(size);
  circle.evaluate(u,SP.gridX(),source,M_PI/2,desc[rank].elements);
  circle.mass(u,desc[rank].elements);
  if( eval.m_bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
    circle.mask(u,rank,size);
  u2 = u;
  if( eval.m_bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
    circle.dssum(u,rank,size);
  else
    circle.periodicDssum(u);
  //    MPIGeometryDotter dotter(circle,desc[rank].elements,rank,size,13+32+2);
  int whole = utilities::g_profiler.add("whole");
  utilities::g_profiler.start(whole);
  if( preconditioned ) {
    eval.m_nosum = &u2;
    pre.m_nosum = &u2;
    std::cout << "iterations " << utilities::CGSolver::solve(u,eval,pre,circle.getDotter(eval.m_bc==SEMLaplacianEvaluator::PERIODIC?true:false),1.e-10) << std::endl;
  }
  else
    std::cout << "iterations " << utilities::CGSolver::solve(u,eval,circle.getDotter(eval.m_bc==SEMLaplacianEvaluator::PERIODIC?true:false),1.e-10) << std::endl;

  utilities::g_profiler.pause(whole);

  eval.m_nu = 1;
  string report = errorReport(u,circle,article,SP.gridX(),M_PI/2,&eval);
  if( rank == 0 ) {
    cout << report;
    cout << utilities::g_profiler.report();
  }

#ifdef HAS_MPI
  MPI_Finalize();
#endif
  return( 0 ); // AY-OH-KAY
}

