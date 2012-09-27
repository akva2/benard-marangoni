#include "mnl/vector.h"
#include "mnl/gll.h"
#include "mnl/matrix.h"
#include "mnl/gordonhall.h"
#include "mnl/function.h"
#include "mnl/cgsolver.h"
#include "mnl/field.h"
#include "mnl/hdf5.h"
#include "mnl/timer.h"
#include "mnl/keyboard.h"

using namespace mnl;
using namespace std;

#include "bm/functions3d.h"
#include "bm/legendrelegendrew.h"
#include "bm/bigcircle.h"
#include "bm/quadratic.h"
#include "bm/consistent.h"
#include "bm/mass.h"
#include "bm/sem.h"
#include "bm/convection.h"

void saveState(const utilities::ringBuffer<basics::Field3<basics::matricesStack> >& u,
    const utilities::ringBuffer<basics::matricesStack>& p, int n, Real Dt)
{
  stringstream str;
  str << "state-ns-pc-" << n << ".hdf5";

  cout << "saving state to " << str.str() << endl;

  HDF5::HDF5Writer writer(str.str());
  str.clear();
  writer.add(u.get(n, 0).X(),"velocity X 1");
  writer.add(u.get(n,-1).X(),"velocity X 2");
  writer.add(u.get(n, 0).Y(),"velocity Y 1");
  writer.add(u.get(n, 0).Y(),"velocity Y 2");
  writer.add(u.get(n, 0).Z(),"velocity Z 1");
  writer.add(u.get(n,-1).Z(),"velocity Z 2");
  writer.add(p.get(n, 0)    ,"pressure 1");
  writer.add(p.get(n,-1)    ,"pressure 2");
  basics::Vector time("time info",2);
  time[0] = n;
  time[1] = Dt;
  writer.add(time);
}

class source3 : public basics::function3D {
  public:
    source3(basics::function3D& n_X, basics::function3D& n_Y,
        basics::function3D& n_Z, basics::function3D& n_phi,
        basics::function3D& n_P, int n_dim) :
      X(n_X), Y(n_Y), Z(n_Z), phi(n_phi), P(n_P), dim(n_dim)
  {
  }

    Real val(Real x, Real y, Real z, Real t) const
    {
      //            if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
      //                return 0;

      Real result = -(phi.diff2x(x,y,z,t)+phi.diff2y(x,y,z,t)+phi.diff2z(x,y,z,t))+phi.difft(x,y,z,t)
        +X.val(x,y,z,t)*phi.diffx(x,y,z,t)+Y.val(x,y,z,t)*phi.diffy(x,y,z,t)
        +Z.val(x,y,z,t)*phi.diffz(x,y,z,t);
      if( dim == 0 )
        result += P.diffx(x,y,z,t);
      if( dim == 1 )
        result += P.diffy(x,y,z,t);
      if( dim == 2 )
        result += P.diffz(x,y,z,t);

      return( result );
    }
  private:
    basics::function3D& X;
    basics::function3D& Y;
    basics::function3D& Z;
    basics::function3D& phi;
    basics::function3D& P;
    int dim;
};

class simplyX : public basics::function3D {
  public:
    simplyX(int index) :
      m_index(index)
  {
  }

    Real val(Real x, Real y, Real z, Real t) const
    {
      switch( m_index ) {
        case  0: return x;
        case  1: return y;
        case  2: return z;
        case  3: return t;
        default: return 0;
      }
    }

    int m_index;
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
  bool preconditioned=true;
  if( argc > 6 )
    preconditioned = atoi(argv[6])==1?false:true;
  bool secondorder=true;
  if( argc > 7 )
    secondorder = atoi(argv[7])==1?false:true;

  Real Lz = 1;

  legendreLegendreW::poissonSolver SP(N,N,0,legendreLegendreW::poissonSolver::Homogenous,true,Lz);
  legendreLegendreW::poissonSolver SP2(N,4,0,legendreLegendreW::poissonSolver::Nonhomogenous,true,Lz);

  basics::Vector gridGL("GL grid",N-1);
  basics::Vector grid("GLL grid",N+1);
  basics::Vector weightGL("GL weight",N-1);
  utilities::GLL::GaussLegendreGrid(gridGL);
  utilities::GLL::GaussLegendreWeights(weightGL,gridGL);
  utilities::GLL::GaussLobattoLegendreGrid(grid);
  basics::Matrix GLL2G = utilities::GLL::interpolationMatrix(grid,gridGL);

  //    bigCircleGeometry circle(N,N,SP.gridX(),SP.weightX());
  quadraticGeometry circle(N,N,grid,SP.weightX());
  circle.setLz(Lz);

  Real nu;
  if( secondorder )
    nu = Real(3)/(2*Dt);
  else
    nu = Real(1)/Dt;

  extrudedConsistentPressureEvaluator eval(circle,SP.Dx(),GLL2G,weightGL,SP.weightX(),grid,gridGL,preconditioned);
  SEMLaplacianEvaluator eval3(circle,SP.Dx(),SP.weightX());
  extrudedLaplacianEvaluator eval2(eval3,SP,nu,&SP2.Ax(),preconditioned);

  convectionEvaluator conv(circle,SP.Dx(),SP.weightX());

  SEMInverseMassEvaluator invmass2D(circle,weightGL,GLL2G);
  extrudedInverseMassEvaluator invmass3D(invmass2D);

  utilities::ringBuffer< basics::Field3<basics::matricesStack> > u;
  u.add(*utilities::g_manager.aquireMatricesStackField("velocity 1",circle));
  u.add(*utilities::g_manager.aquireMatricesStackField("velocity 2",circle));
  u.add(*utilities::g_manager.aquireMatricesStackField("velocity 3",circle));
  basics::Field3<basics::matricesStack>* ut = utilities::g_manager.aquireMatricesStackField("buffer",circle);

  basics::Field3<basics::matricesStack>& F  = *utilities::g_manager.aquireMatricesStackField("source",circle);

  utilities::ringBuffer<basics::matricesStack> p;
  p.add(*utilities::g_manager.aquireMatricesStack("pressure 1",N-1,N-1,N-1,circle.size()));
  p.add(*utilities::g_manager.aquireMatricesStack("pressure 2",N-1,N-1,N-1,circle.size()));
  p.add(*utilities::g_manager.aquireMatricesStack("pressure 3",N-1,N-1,N-1,circle.size()));
  basics::matricesStack& pt  = *utilities::g_manager.aquireMatricesStack("pressure temp",N-1,N-1,N-1,circle.size());
  basics::matricesStack& pt2 = *utilities::g_manager.aquireMatricesStack("pressure temp",N-1,N-1,N-1,circle.size());

  //    articleTest1X articleX;
  //    articleTest1Y articleY;
  //    articleTest1Z articleZ;
  //    articleTest1P articleP;
  boxTest3DX articleX;
  boxTest3DY articleY;
  boxTest3DZ articleZ;
  boxTest3DP articleP;
  source3 sourceX(articleX,articleY,articleZ,articleX,articleP,0);
  source3 sourceY(articleX,articleY,articleZ,articleY,articleP,1);
  source3 sourceZ(articleX,articleY,articleZ,articleZ,articleP,2);

  /* set up initial conditions */
  Real t0=0;
  circle.evaluate(u.get(0,-1),grid,articleX,articleY,articleZ,t0-Dt);
  circle.evaluate(u.get(0, 0),grid,articleX,articleY,articleZ,t0);
  circle.evaluate(p.get(0,-1),gridGL,articleP,t0-Dt);
  circle.evaluate(p.get(0, 0),gridGL,articleP,t0);

  /* setup keyboard */
  utilities::set_raw_tty();

  int whole = utilities::g_profiler.add("whole");
  int n;
  for( n=0;n<(T-t0)/Dt;++n ) {
    std::cout << "Time: " << t0+(n+1)*Dt << std::endl;

    /* check keyboard */
    if( utilities::keypress() == 's' )
      saveState(u,p,n,Dt);
    if( utilities::keypress() == 'q' )
      break;

    /* set up references */
    basics::Field3<basics::matricesStack>& unm1	= u.get(n,-1);
    basics::Field3<basics::matricesStack>& un  	= u.get(n, 0);
    basics::Field3<basics::matricesStack>& unp1	= u.get(n, 1);
    basics::matricesStack& pnm1 				= p.get(n,-1);
    basics::matricesStack& pn  					= p.get(n, 0);
    basics::matricesStack& pnp1 				= p.get(n, 1);

    circle.evaluate(F,grid,sourceX,sourceY,sourceZ,t0+(n+1)*Dt);
    utilities::g_profiler.start(whole);

    /* solve convection problems */
    conv.solve(unp1,un,u,n,Dt,Dt,0,Dt);
    if( secondorder )
      conv.solve(*ut,unm1,u,n,Dt,Dt,-Dt,Dt);

    /* BDF2 */
    if( secondorder ) {
      F.axpy(Real(2)/Dt,unp1);
      F.axpy(Real(-1)/(2*Dt),*ut);
    } else
      F.axpy(Real(1)/Dt,unp1);

    circle.mass(F,false);

    /* extrapolate pressure */
    if( secondorder )
      utilities::extrapolate(pnp1,p,n,utilities::FIRST_ORDER);
    else
      utilities::extrapolate(pnp1,p,n,utilities::ZEROTH_ORDER);

    eval.m_gradient.evaluate(unp1,pnp1);
    unp1 += F;

    /* delta */
    eval2.evaluate(unm1,un,false,false);
    unp1 -= unm1;

    circle.maskField(unp1);
    F = unp1;
    circle.dssum(unp1);
    eval2.solve(unp1,F);
    unp1 += un;

    eval.m_divergence.evaluate(pt,unp1);
    invmass3D.evaluate(pt2,pt);
    pt *= -nu;

    //        if( preconditioned )
    //            cout << "iterations pressure " << utilities::CGSolver::solve(pt2,eval,eval4,pt2.getDotter(),1.e-8) << endl;
    //        else
    //            cout << "iterations pressure " << utilities::CGSolver::solve(pt2,eval,pt2.getDotter(),1.e-8) << endl;
    eval.solve(pt);

    pnp1 += pt;

    /* update velocity */
    eval.m_gradient.evaluate(unm1,pt);
    dssum(unm1,circle,0,1);
    mask(unm1,circle,0,1,eval2.m_bc);
    circle.invMass(unm1);
    unp1.axpy(Real(1)/nu,unm1);
    utilities::g_profiler.pause(whole);

    /* rotational */
    //        pnp1 -= pt2;
  }

  /* restore keyboard */
  utilities::set_normal_tty();

  saveState(u,p,n,Dt);

  eval2.m_nu = 0;
  cout << errorReport(u.get(n),p.get(n),circle,articleX,articleY,articleZ,
      articleP,grid,gridGL,n*Dt,invmass3D,&eval2);

  cout << utilities::g_profiler.report();
  return 0; // AY-OH-KAY
}

