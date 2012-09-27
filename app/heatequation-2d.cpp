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
#include "bm//functions2d.h"
#include "bm//poissonsolver-lgw.h"
#include "bm/laplacian.h"
#include "bm/sem.h"

class source2 : public basics::function2D {
  public:
    source2(basics::function2D& n_X,Real n_nu) : X(n_X), nu(n_nu) 
  {
  }

    Real val(Real x, Real y, Real t) const
    {
      if( x == 0 && y == 0 )
        return 0;

      return( X.difft(x,y,t) - nu*(X.diff2x(x,y,t)+X.diff2y(x,y,t)) );
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
    preconditioned=(atoi(argv[6])==1?false:true);

  legendreLegendreW::poissonSolver SP(N,4,0,legendreLegendreW::poissonSolver::Homogenous,true,-1);

  bigCircleGeometry circle(N,N,SP.gridX(),SP.weightX());

  SEMLaplacianEvaluator eval(circle,SP.Dx(),SP.weightX());
  SEMFEMLaplacianPreconditioner pre(circle,SP.gridX(),SP.weightX(),Real(3)/(Real(2)*Dt));
  eval.m_nu = Real(3)/(Real(2)*Dt);

  utilities::ringBuffer<basics::matrixStack> u;
  u.add(*utilities::g_manager.aquireMatrixStack("velocity",circle));
  u.add(*utilities::g_manager.aquireMatrixStack("velocity",circle));
  u.add(*utilities::g_manager.aquireMatrixStack("velocity",circle));
  basics::matrixStack& buf = *utilities::g_manager.aquireMatrixStack(  "buffer",circle);

  articleTest1X article;
  source2 source(article,1);

  /* set up initial conditions */
  Real t0=0;
  circle.evaluate(u.get(0,-1),SP.gridX(),article,t0-Dt);
  circle.evaluate(u.get(0, 0),SP.gridX(),article,t0);
  circle.mask(u.get(0, 0));
  circle.mask(u.get(0,-1));

  int n;
  for( n=0;n<floor((T-t0)/Dt);++n ) {
    std::cout << "Time: " << t0+(n+1)*Dt << std::endl;

    /* set up references */
    basics::matrixStack& unm1 	= u.get(n,-1);
    basics::matrixStack& un  	= u.get(n, 0);
    basics::matrixStack& unp1 	= u.get(n, 1);

    /* source */
    circle.evaluate(unp1,SP.gridX(),source,t0+(n+1)*Dt);

    /* BDF2 */
    unp1.axpy(Real(2)/Dt,un);
    unp1.axpy(Real(-1)/(Real(2)*Dt),unm1);

    circle.mass(unp1,false);

    /* delta value */
    eval.m_nosum = &buf;	
    eval.evaluate(unm1,un,false,true);
    unp1 -= buf;

    circle.mask(unp1);
    buf = unp1;
    circle.dssum(unp1);

    if( preconditioned )
      std::cout << "iterations " << utilities::CGSolver::solve(unp1,eval,pre,circle.getDotter(),1.e-10) << std::endl;
    else
      std::cout << "iterations " << utilities::CGSolver::solve(unp1,eval,circle.getDotter(),1.e-10) << std::endl;
    unp1 += un;
  }

  eval.m_nu = 1;
  std::cout << errorReport(u.get(n),circle,article,SP.gridX(),T,&eval);

  return 0; // AY-OH-KAY
}

