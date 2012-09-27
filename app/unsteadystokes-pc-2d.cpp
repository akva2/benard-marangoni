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
#include "bm/mass.h"
#include "bm/consistent.h"
#include "bm/sem.h"

#ifdef HAS_MPI
#include <mpi.h>
#endif

class shenTest2DXS : public basics::function2D {
  public:
    shenTest2DXS()
    {
      m_id = 2;
    }

    Real val(Real x, Real y, Real t) const
    {
      return M_PI*sin(2*M_PI*y)*sin(M_PI*x)*sin(M_PI*x)*sin(t);
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return 2*M_PI*M_PI*sin(2*M_PI*y)*sin(M_PI*x)*cos(M_PI*x)*sin(t);
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return 2*M_PI*M_PI*cos(2*M_PI*y)*sin(M_PI*x)*sin(M_PI*x)*sin(t);
    }

    Real difft(Real x, Real y, Real t) const
    {
      return M_PI*sin(2*M_PI*y)*sin(M_PI*x)*sin(M_PI*x)*cos(t);
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return 2*pow(M_PI,3.f)*sin(2*M_PI*y)*cos(M_PI*x)*cos(M_PI*x)*sin(t)-2*pow(M_PI,3.f)*sin(2*M_PI*y)*sin(M_PI*x)*sin(M_PI*x)*sin(t);
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return -4*pow(M_PI,3.f)*sin(2*M_PI*y)*sin(M_PI*x)*sin(M_PI*x)*sin(t);
    }
};

class shenTest2DYS : public basics::function2D {
  public:
    shenTest2DYS()
    {
      m_id = 2;
    }

    Real val(Real x, Real y, Real t) const
    {
      return( -M_PI*sin(2*M_PI*x)*sin(M_PI*y)*sin(M_PI*y)*sin(t) );
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return -2*M_PI*M_PI*cos(2*M_PI*x)*sin(M_PI*y)*sin(M_PI*y)*sin(t);
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return -2*M_PI*M_PI*sin(2*M_PI*x)*sin(M_PI*y)*cos(M_PI*y)*sin(t);
    }

    Real difft(Real x, Real y, Real t) const
    {
      return( -M_PI*sin(2*M_PI*x)*sin(M_PI*y)*sin(M_PI*y)*cos(t) );
    }

    Real diff2x(Real x, Real y, Real t) const
    {
      return 4*pow(M_PI,3.f)*sin(2*M_PI*x)*sin(M_PI*y)*sin(M_PI*y)*sin(t);
    }

    Real diff2y(Real x, Real y, Real t) const
    {
      return -2*pow(M_PI,3.f)*sin(2*M_PI*x)*cos(M_PI*y)*cos(M_PI*y)*sin(t) + 2*pow(M_PI,3.f)*sin(2*M_PI*x)*sin(M_PI*y)*sin(M_PI*y)*sin(t);
    }
};

class shenTest2DPS : public basics::function2D {
  public:
    shenTest2DPS()
    {
      m_id = 2;
    }

    Real val(Real x, Real y, Real t) const
    {
      return cos(M_PI*x)*sin(M_PI*y)*sin(t);
    }

    Real diffx(Real x, Real y, Real t) const
    {
      return -M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(t);
    }

    Real diffy(Real x, Real y, Real t) const
    {
      return M_PI*cos(M_PI*x)*cos(M_PI*y)*sin(t);
    }
};

class source2 : public basics::function2D {
  public:
    source2(basics::function2D& n_X, basics::function2D& n_P, int n_dir) :
      X(n_X), P(n_P), dir(n_dir)
  {
  }

    Real val(Real x, Real y, Real t) const
    {
      //            if( abs(x) < 1.e-14 && abs(y) < 1.e-14 && dir == 0)
      //                return 0;
      //            if( abs(x) < 1.e-14 && abs(y) < 1.e-14 && dir == 1)
      //                return -Real(8)/5*M_PI*M_PI*sin(t);

      Real result = -(X.diff2x(x,y,t)+X.diff2y(x,y,t))+X.difft(x,y,t);
      if( dir == 0 )
        result += P.diffx(x,y,t);
      if( dir == 1 )
        result += P.diffy(x,y,t);

      return( result );
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
  Real Dt = 1.e-2f;
  Real T = 1.f;
  if( argc > 3 )
    Dt = Real(atoi(argv[2]))/atoi(argv[3]);
  if( argc > 5 )
    T = Real(atoi(argv[4]))/atoi(argv[5]);
  bool preconditioned = true;
  if( argc > 6 )
    preconditioned = (atoi(argv[6])==1?false:true);
  bool secondorder = true;
  if( argc > 7 )
    preconditioned = (atoi(argv[7])==1?true:false);

  legendreLegendreW::poissonSolver SP(N,N,0,legendreLegendreW::poissonSolver::Homogenous,true,-1);

  basics::Vector gridGL("GL grid",N-1);
  basics::Vector weightGL("GL weight",N-1);
  utilities::GLL::GaussLegendreGrid(gridGL);
  utilities::GLL::GaussLegendreWeights(weightGL,gridGL);
  basics::Matrix GLL2G = utilities::GLL::interpolationMatrix(SP.gridX(),gridGL);

  bigCircleGeometry circle(N,N,SP.gridX(),SP.weightX());

  SEMConsistentPressureEvaluator eval(circle,SP.Dx(),GLL2G,weightGL,SP.weightX(),NULL,0,SEMLaplacianEvaluator::PERIODIC);
  SEMFEMConsistentPressurePreconditioner pre3(circle,SP.weightX(),weightGL,SP.gridX(),gridGL,0,rank,size);
  SEMLaplacianEvaluator eval2(circle,SP.Dx(),SP.weightX(),NULL,rank,size,SEMLaplacianEvaluator::PERIODIC);
  if( secondorder )
    eval2.m_nu = Real(3)/(Real(2)*Dt);
  else
    eval2.m_nu = Real(1)/(Dt);
  SEMFEMLaplacianPreconditioner pre(circle,SP.gridX(),SP.weightX(),eval2.m_nu,rank,size);

  utilities::ringBuffer< basics::Field2<basics::matrixStack> > u;
  basics::coarseGrid desc = circle.getDivisionInfo(size);
  u.add(*utilities::g_manager.aquireMatrixStackField("velocity 1",N+1,N+1,desc[rank].elements.size()));
  u.add(*utilities::g_manager.aquireMatrixStackField("velocity 2",N+1,N+1,desc[rank].elements.size()));
  u.add(*utilities::g_manager.aquireMatrixStackField("velocity 3",N+1,N+1,desc[rank].elements.size()));

  basics::Field2<basics::matrixStack>& F  = *utilities::g_manager.aquireMatrixStackField("source",N+1,N+1,desc[rank].elements.size());

  utilities::ringBuffer<basics::matrixStack> p;
  p.add(*utilities::g_manager.aquireMatrixStack("pressure 1",N-1,N-1,desc[rank].elements.size()));
  p.add(*utilities::g_manager.aquireMatrixStack("pressure 2",N-1,N-1,desc[rank].elements.size()));
  p.add(*utilities::g_manager.aquireMatrixStack("pressure 3",N-1,N-1,desc[rank].elements.size()));

  //    articleTest1X articleX;
  //    articleTest1Y articleY;
  //    articleTest1P articleP;
  shenTest2DXS articleX;
  shenTest2DYS articleY;
  shenTest2DPS articleP;
  source2 sourceX(articleX,articleP,0);
  source2 sourceY(articleY,articleP,1);

  /* set up initial conditions */
  Real t0=0;
  circle.evaluate(u.get(0,-1),SP.gridX(),articleX,articleY,t0-Dt,desc[rank].elements);
  circle.evaluate(u.get(0, 0),SP.gridX(),articleX,articleY,t0,desc[rank].elements);
  circle.evaluate(p.get(0,-1),gridGL,articleP,t0-Dt,desc[rank].elements);
  circle.evaluate(p.get(0, 0),gridGL,articleP,t0,desc[rank].elements);

  basics::geometryStack::geometryDotter dotter = circle.getDotter(eval2.m_bc==SEMLaplacianEvaluator::PERIODIC?true:false);
  SEMInverseMassEvaluator invmass(circle,weightGL,GLL2G,rank,size);

  cout << "time step:\t\t"  	 << Dt << endl
    <<	"final time:\t\t" 	 <<  T << endl
    <<	"preconditioned:\t\t"<< (preconditioned?"yes":"no")  << endl
    <<	"order :\t\t\t"	 	 << (secondorder?2:1) << endl;

  int n;
  int whole = utilities::g_profiler.add("whole");
  for( n=0;n<floor((T-t0)/Dt);++n ) {
    std::cout << "Time: " << t0+(n+1)*Dt << std::endl;

    /* set up references */
    basics::Field2<basics::matrixStack>& unm1 	= u.get(n,-1);
    basics::Field2<basics::matrixStack>& un  	= u.get(n, 0);
    basics::Field2<basics::matrixStack>& unp1 	= u.get(n, 1);
    basics::matrixStack& pnm1 					= p.get(n,-1);
    basics::matrixStack& pn  					= p.get(n, 0);
    basics::matrixStack& pnp1 					= p.get(n, 1);

    circle.evaluate(F,SP.gridX(),sourceX,sourceY,t0+(n+1)*Dt,desc[rank].elements);
    utilities::g_profiler.start(whole);

    /* BDF2 */
    if( secondorder ) {
      F.axpy(Real(2)/Dt,un);
      F.axpy(Real(-1)/(2*Dt),unm1);
    } else
      F.axpy(Real(1)/Dt,un);
    circle.mass(F,desc[rank].elements);

    /* extrapolate pressure */
    if( secondorder )
      utilities::extrapolate(pnp1,p,n,utilities::FIRST_ORDER);
    else
      utilities::extrapolate(pnp1,p,n,utilities::ZEROTH_ORDER);

    eval.m_gradient.evaluate(unp1,pnp1);
    unp1 += F;

    /* delta */
    eval2.evaluate(F,un,false,false);
    unp1 -= F;

    if( eval2.m_bc != SEMLaplacianEvaluator::PERIODIC )
      circle.maskField(unp1,rank,size);
    unm1 = unp1;
    if( eval2.m_bc == SEMLaplacianEvaluator::PERIODIC ) {
      circle.periodicDssum(unp1.X());
      circle.periodicDssum(unp1.Y());
    } else
      circle.dssum(unp1,rank,size);
    if( preconditioned ) {
      eval2.m_nosum = pre.m_nosum = &unm1.X();
      cout << "iterations velocity x " << utilities::CGSolver::solve(unp1.X(),eval2,pre,dotter,1.e-8) << endl;
      eval2.m_nosum = pre.m_nosum = &unm1.Y();
      cout << "iterations velocity y " << utilities::CGSolver::solve(unp1.Y(),eval2,pre,dotter,1.e-8) << endl;
    } else {
      eval2.m_nosum = NULL;
      cout << "iterations velocity x " << utilities::CGSolver::solve(unp1.X(),eval2,dotter,1.e-8) << endl;
      cout << "iterations velocity y " << utilities::CGSolver::solve(unp1.Y(),eval2,dotter,1.e-8) << endl;
    }
    unp1 += un;
    circle.maskField(unp1);
    eval.m_divergence.evaluate(pnm1,unp1);
    pnm1 *= -eval2.m_nu;

    if( preconditioned ) {
      cout << "iterations pressure " 
        << utilities::CGSolver::solve(pnm1,eval,pre3,pnm1.getDotter(),1.e-8) 
        << endl;
      //            eval.m_deflated = true;
      //            eval.solve(pnm1,&pre3);
    } else {
      cout << "iterations pressure " 
        << utilities::CGSolver::solve(pnm1,eval,pnm1.getDotter(),1.e-8) 
        << endl;
      //            cout << "iterations pressure " << eval.solve(pnm1) << endl;
    }

    pnp1 += pnm1;

    /* PCR */
    eval.m_gradient.evaluate(unm1,pnm1);
    if( eval2.m_bc == SEMLaplacianEvaluator::PERIODIC ) {
      circle.periodicDssum(unm1.X());
      circle.periodicDssum(unm1.Y());
    } else
      circle.dssum(unm1,rank,size);

    /* update velocity */
    if( eval2.m_bc == SEMLaplacianEvaluator::PERIODIC ) {
      circle.invMassP(unm1.X());
      circle.invMassP(unm1.Y());
    } else {
      circle.invMass(unm1,desc[rank].elements);
      circle.maskField(unm1,rank,size);
    }

    unp1.axpy(Real(1)/eval2.m_nu,unm1);
    utilities::g_profiler.pause(whole);
  }

  eval2.m_nu = 1;

  string report = errorReport(u.get(n),p.get(n),circle,articleX,articleY,articleP,SP.gridX(),
      gridGL,T,invmass,&eval2);
  if( rank == 0 ) {
    cout << report;
    cout << utilities::g_profiler.report();
  }

#ifdef HAS_MPI
  MPI_Finalize();
#endif

  return 0; // AY-OH-KAY
}

