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
#include "bm/functions2d.h"
#include "bm/uzawa.h"
#include "bm/mass.h"
#include "bm/sem.h"

#ifdef HAS_MPI
#include <mpi.h>
#endif

class source2 : public basics::function2D {
  public:
    source2(basics::function2D& n_X, basics::function2D& n_P, int n_dir) :
      X(n_X), P(n_P), dir(n_dir)
  {
  }

    Real val(Real x, Real y, Real t) const
    {
      if( x == 0 && y == 0 && dir == 0)
        return 0;
      if( x == 0 && y == 0 && dir == 1)
        return -Real(8)/5*M_PI*M_PI;

      Real result = -(X.diff2x(x,y,t)+X.diff2y(x,y,t));
      if( dir == 0 )
        result += P.diffx(x,y,t);
      if( dir == 1 )
        result += P.diffy(x,y,t);

      return( result );
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
      return( 0.f ); 
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return( 0 );
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return( 0 );
    }
  private:
    basics::function2D& X;
    basics::function2D& P;
    int dir;
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
    preconditioned = (atoi(argv[2])==1?false:true);


  legendreLegendreW::poissonSolver SP(N,N,0,legendreLegendreW::poissonSolver::Homogenous,true,-1);

  basics::Vector gridGL("GL grid",N-1);
  basics::Vector weightGL("GL weight",N-1);
  utilities::GLL::GaussLegendreGrid(gridGL);
  utilities::GLL::GaussLegendreWeights(weightGL,gridGL);
  basics::Matrix GLL2G = utilities::GLL::interpolationMatrix(SP.gridX(),gridGL);

  bigCircleGeometry circle(N,N,SP.gridX(),SP.weightX());

  SEMUzawaEvaluator 		 eval(circle,SP.Dx(),GLL2G,SP.gridX(),SP.weightX(),weightGL,preconditioned,rank,size);
  SEMInverseMassEvaluator  eval2(circle,weightGL,GLL2G,rank,size);

  basics::coarseGrid desc = circle.getDivisionInfo(size);
  basics::Field2<basics::matrixStack>& u  = *utilities::g_manager.aquireMatrixStackField("velocity 1",N+1,N+1,desc[rank].elements.size());
  basics::Field2<basics::matrixStack>& u2 = *utilities::g_manager.aquireMatrixStackField("velocity 2",N+1,N+1,desc[rank].elements.size());

  basics::Field2<basics::matrixStack>& F  = *utilities::g_manager.aquireMatrixStackField("source",N+1,N+1,desc[rank].elements.size());

  basics::matrixStack& p  = *utilities::g_manager.aquireMatrixStack("pressure 1",N-1,N-1,desc[rank].elements.size());

  articleTest1X articleX;
  articleTest1Y articleY;
  articleTest1P articleP;
  source2 sourceX(articleX,articleP,0);
  source2 sourceY(articleY,articleP,1);
  circle.evaluate(F,SP.gridX(),sourceX,sourceY,M_PI/2,desc[rank].elements);
  int whole = utilities::g_profiler.add("whole");
  utilities::g_profiler.start(whole);
  /* B */
  circle.mass(F,desc[rank].elements);
  u = F;
  circle.maskField(u,rank,size);
  u2 = u;
  circle.dssum(u,rank,size);
  /* A^-1 */
  eval.m_laplacian.m_nu = 0;
  MPIGeometryDotter dotter(circle,desc[rank].elements,rank,size);
  MPIDotter pdotter(rank,size);
  if( preconditioned ) {
    *eval.m_laplacian.m_nosum = u2.X();
    utilities::CGSolver::solve(u.X(),eval.m_laplacian,*eval.m_pre,dotter,1.e-8);
    *eval.m_laplacian.m_nosum = u2.Y();
    utilities::CGSolver::solve(u.Y(),eval.m_laplacian,*eval.m_pre,dotter,1.e-8);
  } else {
    eval.m_laplacian.m_nosum = NULL;
    utilities::CGSolver::solve(u.X(),eval.m_laplacian,dotter,1.e-8);
    utilities::CGSolver::solve(u.Y(),eval.m_laplacian,dotter,1.e-8);
  }
  /* -D */
  eval.m_divergence.evaluate(p,u);
  p *= -1;

  std::cout << "iterations pressure " << utilities::CGSolver::solve(p,eval,eval2,pdotter,1.e-8) << std::endl;

  eval.m_gradient.evaluate(u,p,false);
  u += F;
  circle.maskField(u,rank,size);
  u2 = u;
  circle.dssum(u,rank,size);

  if( preconditioned ) {
    eval.m_laplacian.m_nosum = eval.m_pre->m_nosum = &u2.X();
    std::cout << "iterations velocity x " << utilities::CGSolver::solve(u.X(),eval.m_laplacian,*eval.m_pre,dotter,1.e-10) << std::endl;
    eval.m_laplacian.m_nosum = eval.m_pre->m_nosum = &u2.Y();
    std::cout << "iterations velocity y " << utilities::CGSolver::solve(u.Y(),eval.m_laplacian,*eval.m_pre,dotter,1.e-10) << std::endl;
  } else {
    std::cout << "iterations velocity x " << utilities::CGSolver::solve(u.X(),eval.m_laplacian,dotter,1.e-10) << std::endl;
    std::cout << "iterations velocity y " << utilities::CGSolver::solve(u.Y(),eval.m_laplacian,dotter,1.e-10) << std::endl;
  }
  utilities::g_profiler.pause(whole);

  eval.m_laplacian.m_nu = 1;
  string report = errorReport(u,p,circle,articleX,articleY,articleP,SP.gridX(),gridGL,M_PI/2,eval2,&eval.m_laplacian);
  if( rank == 0 ) {
    cout << report;
    cout << utilities::g_profiler.report();
  }

#ifdef HAS_MPI
  MPI_Finalize();
#endif

  return 0; // AY-OH-KAY
}
