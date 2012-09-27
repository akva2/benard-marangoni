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
#include "bm/poissonsolver-lgw.h"
#include "bm/laplacian.h"
#include "bm/sem.h"
#include "bm/convection.h"

class source : public basics::function2D {
  public:
    source(basics::function2D& n_X, basics::function2D& n_Y,
           basics::function2D& n_phi, Real n_nu) : 
      X(n_X), Y(n_Y), phi(n_phi), nu(n_nu) 
  {
  }

    Real val(Real x, Real y, Real t) const
    {
      if( x == 0 && y == 0 )
        return 0;

      return( phi.difft(x,y,t) - nu*(phi.diff2x(x,y,t)+phi.diff2y(x,y,t)) + X.val(x,y,t)*phi.diffx(x,y,t)+Y.val(x,y,t)*phi.diffy(x,y,t) );
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
    basics::function2D& Y;
    basics::function2D& phi;
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

  convectionEvaluator conv(circle,SP.Dx(),SP.weightX());

  utilities::ringBuffer< basics::Field2<basics::matrixStack> > u;
  u.add(*utilities::g_manager.aquireMatrixStackField("velocity",circle));
  u.add(*utilities::g_manager.aquireMatrixStackField("velocity",circle));
  u.add(*utilities::g_manager.aquireMatrixStackField("velocity",circle));
  basics::Field2<basics::matrixStack>* F = utilities::g_manager.aquireMatrixStackField("buffer",circle);
  basics::Field2<basics::matrixStack>* buf = utilities::g_manager.aquireMatrixStackField("buffer",circle);

  articleTest1X articleX;
  articleTest1Y articleY;
  source sourceX(articleX,articleY,articleX,1);
  source sourceY(articleX,articleY,articleY,1);

  /* set up initial conditions */
  Real t0=0;
  circle.evaluate(u.get(0,-1),SP.gridX(),articleX,articleY,t0-Dt);
  circle.evaluate(u.get(0, 0),SP.gridX(),articleX,articleY,t0);
  circle.maskField(u.get(0, 0));
  circle.maskField(u.get(0,-1));

  int n;
  for( n=0;n<floor((T-t0)/Dt);++n ) {
    std::cout << "Time: " << t0+(n+1)*Dt << std::endl;

    /* set up references */
    basics::Field2<basics::matrixStack>& unm1 	= u.get(n,-1);
    basics::Field2<basics::matrixStack>& un  	= u.get(n, 0);
    basics::Field2<basics::matrixStack>& unp1 	= u.get(n, 1);

    /* source */
    circle.evaluate(*F,SP.gridX(),sourceX,sourceY,t0+(n+1)*Dt);

    /* solve convection problems */
    conv.solve(unp1.X(),un.X(),u,n,Dt,Dt,0);
    conv.solve(unp1.Y(),un.Y(),u,n,Dt,Dt,0);

    conv.solve(buf->X(),unm1.X(),u,n,Dt,Dt,-Dt);
    conv.solve(buf->Y(),unm1.Y(),u,n,Dt,Dt,-Dt);

    /* BDF2 */
    unp1 *= Real(2)/Dt;
    unp1.axpy(Real(-1)/(Real(2)*Dt),*buf);
    unp1 += *F;

    circle.mass(unp1,false);

    /* delta value */
    eval.m_nosum = &buf->X();
    eval.evaluate(unm1.X(),un.X(),false,true);
    eval.m_nosum = &buf->Y();
    eval.evaluate(unm1.Y(),un.Y(),false,true);
    unp1 -= *buf;

    circle.maskField(unp1);
    *buf = unp1;
    circle.dssum(unp1);

    eval.m_nosum = pre.m_nosum = &buf->X();
    if( preconditioned )
      std::cout << "iterations " << utilities::CGSolver::solve(unp1.X(),eval,pre,circle.getDotter(),1.e-10) << std::endl;
    else
      std::cout << "iterations " << utilities::CGSolver::solve(unp1.X(),eval,circle.getDotter(),1.e-10) << std::endl;
    eval.m_nosum = pre.m_nosum = &buf->Y();
    if( preconditioned )
      std::cout << "iterations " << utilities::CGSolver::solve(unp1.Y(),eval,pre,circle.getDotter(),1.e-10) << std::endl;
    else
      std::cout << "iterations " << utilities::CGSolver::solve(unp1.Y(),eval,circle.getDotter(),1.e-10) << std::endl;
    unp1 += un;
  }

  eval.m_nu = 1;
  cout << errorReport(u.get(n),circle,articleX,articleY,SP.gridX(),T,&eval);

  return 0; // AY-OH-KAY
}

