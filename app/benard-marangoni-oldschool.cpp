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

class initialTemperature : public basics::function2D {
  public:
    initialTemperature(const Real theta0) :
      m_theta0(theta0)
  {
    srand(time(NULL));
  }

    Real val(Real x, Real y, Real t) const
    {
      return (m_theta0*rand())/(RAND_MAX);
    }
  protected:
    Real m_theta0;
};

void loadState(utilities::ringBuffer<basics::Field3<basics::matricesStack> >& u,
    utilities::ringBuffer<basics::matricesStack>& p,
    utilities::ringBuffer<basics::matricesStack>& theta,
    int& n, Real& Dt, Real& Ra, Real& Ma, Real& Pr, Real& Lz,
    int rank, const string& file)
{
  stringstream str;
  str << file << "-" << rank << ".hdf5";

  cout << "loading state from " << str.str() << endl;

  HDF5::HDF5Reader reader(str.str());
  basics::Vector time("time info",6);
  reader.read(time,"info");
  n = time[0];
  Dt = time[1];
  Ra = time[2];
  Ma = time[3];
  Pr = time[4];
  Lz = time[5];

  reader.read(u.get(n, 0).X(),"velocity X 1");
  reader.read(u.get(n,-1).X(),"velocity X 2");
  reader.read(u.get(n, 0).Y(),"velocity Y 1");
  reader.read(u.get(n,-1).Y(),"velocity Y 2");
  reader.read(u.get(n, 0).Z(),"velocity Z 1");
  reader.read(u.get(n,-1).Z(),"velocity Z 2");
  reader.read(p.get(n, 0)    ,"pressure 1");
  reader.read(p.get(n,-1)    ,"pressure 2");
  reader.read(theta.get(n, 0),"temperature 1");
  reader.read(theta.get(n,-1),"temperature 2");
}

int addn;

void saveState(const utilities::ringBuffer<basics::Field3<basics::matricesStack> >& u,
    const utilities::ringBuffer<basics::matricesStack>& p,
    const utilities::ringBuffer<basics::matricesStack>& theta, int n,
    Real Dt, Real Ra, Real Ma, Real Pr, Real Lz, int rank)
{
  stringstream str;
  str << "state-rbN-" << n+addn << "-" << rank << ".hdf5";

  cout << "saving state to " << str.str() << endl;

  HDF5::HDF5Writer writer(str.str());
  str.clear();
  writer.add(u.get(n, 0).X(),"velocity X 1");
  writer.add(u.get(n,-1).X(),"velocity X 2");
  writer.add(u.get(n, 0).Y(),"velocity Y 1");
  writer.add(u.get(n,-1).Y(),"velocity Y 2");
  writer.add(u.get(n, 0).Z(),"velocity Z 1");
  writer.add(u.get(n,-1).Z(),"velocity Z 2");
  writer.add(p.get(n, 0)    ,"pressure 1");
  writer.add(p.get(n,-1)    ,"pressure 2");
  writer.add(theta.get(n, 0),"temperature 1");
  writer.add(theta.get(n,-1),"temperature 2");
  basics::Vector time("info",6);
  time[0] = n;
  time[1] = Dt;
  time[2] = Ra;
  time[3] = Ma;
  time[4] = Pr;
  time[5] = Lz;
  writer.add(time);
}

int main(int argc, char** argv)
{
  int rank=0, size=1;

#ifdef HAS_MPI
  int aquired;
  MPI_Init_thread(&argc, &argv,MPI_THREAD_SERIALIZED,&aquired);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
#endif
  int N;
  if( argc < 2)
    N=12;
  else
    N = atoi(argv[1]);

  Real Dt = 1.e-2f;
  Real T = 1.f;
  bool load=false;
  if( argc > 2 && strstr(argv[2],"state"))
    load = true;
  if( argc > 3 && !load)
    Dt = Real(atoi(argv[2]))/atoi(argv[3]);
  if( argc > 5 || (load && argc > 3))
    T = Real(atoi(argv[4-load]))/atoi(argv[5-load]);
  Real Ra = 27.f;
  Real Ma = 101.f;
  Real Pr = 890;
  if (argc > 6 )
    Ra = atof(argv[6]);
  if (argc > 7 )
    Ma = atof(argv[7]);
  if (argc > 8 )
    Pr = atof(argv[8]);
  Real Lz=1;
  if( argc > 9 )
    Lz=atof(argv[9]);;
  bool preconditioned=true;
  if( argc > 10 )
    preconditioned = atoi(argv[10])==1?false:true;
  bool lumped=true;
  if( argc > 10 )
    lumped = atoi(argv[11])==1?false:true;
  int symmetry=0;
  if( argc > 11 )
    symmetry = atoi(argv[12]);

  legendreLegendreW::poissonSolver SP(N,N,0,legendreLegendreW::poissonSolver::Homogenous,true,Lz);
  legendreLegendreW::poissonSolver SP2(N,4,0,legendreLegendreW::poissonSolver::Nonhomogenous,true,Lz);
  legendreLegendreW::poissonSolver SP3(N,4,0,legendreLegendreW::poissonSolver::HomogenousLeft,true,Lz);

  basics::Vector gridGL("GL grid",N-1);
  basics::Vector grid("GLL grid",N+1);
  basics::Vector weightGL("GL weight",N-1);
  utilities::GLL::GaussLegendreGrid(gridGL);
  utilities::GLL::GaussLobattoLegendreGrid(grid);
  utilities::GLL::GaussLegendreWeights(weightGL,gridGL);
  basics::Matrix GLL2G = utilities::GLL::interpolationMatrix(grid,gridGL);

  bigCircleGeometry circle(N,N,grid,SP.weightX());

  basics::coarseGrid desc = circle.getDivisionInfo(size);
  utilities::ringBuffer< basics::Field3<basics::matricesStack> > u;
  u.add(*utilities::g_manager.aquireMatricesStackField("velocity 1",N+1,N+1,N+1,
        desc[rank].elements.size()));
  u.add(*utilities::g_manager.aquireMatricesStackField("velocity 2",N+1,N+1,N+1,
        desc[rank].elements.size()));
  u.add(*utilities::g_manager.aquireMatricesStackField("velocity 3",N+1,N+1,N+1,
        desc[rank].elements.size()));
  basics::Field3<basics::matricesStack>& buf =
    *utilities::g_manager.aquireMatricesStackField("buffer",N+1,N+1,N+1,
        desc[rank].elements.size());

  utilities::ringBuffer<basics::matricesStack> p;
  p.add(*utilities::g_manager.aquireMatricesStack("pressure 1",N-1,N-1,N-1,
        desc[rank].elements.size()));
  p.add(*utilities::g_manager.aquireMatricesStack("pressure 2",N-1,N-1,N-1,
        desc[rank].elements.size()));
  p.add(*utilities::g_manager.aquireMatricesStack("pressure 3",N-1,N-1,N-1,
        desc[rank].elements.size()));
  basics::matricesStack& pt =
    *utilities::g_manager.aquireMatricesStack("pressure temp",N-1,N-1,N-1,
        desc[rank].elements.size());

  utilities::ringBuffer<basics::matricesStack> theta;
  theta.add(*utilities::g_manager.aquireMatricesStack("temperature 1",N+1,N+1,N+1,
        desc[rank].elements.size()));
  theta.add(*utilities::g_manager.aquireMatricesStack("temperature 2",N+1,N+1,N+1,
        desc[rank].elements.size()));
  theta.add(*utilities::g_manager.aquireMatricesStack("temperature 3",N+1,N+1,N+1,
        desc[rank].elements.size()));
  int n=0;
  Real t0=0;
  if( load )
    loadState(u,p,theta,n,Dt,Ra,Ma,Pr,Lz,rank,argv[2]);
  else {

    /* set up initial conditions */
    Real theta0=1.e-3*100;
    initialTemperature temp(theta0);

    circle.evaluate(theta.get(0).at(0),grid,temp,t0,desc[rank].elements);

    dssum(theta.get(0),circle,rank,size);
    for( int n=0;n<desc[rank].elements.size();++n )
      basics::multPointwise(buf.X().at(0)[n],theta.get(0).at(0)[n],(*circle.m_mult)[desc[rank].elements[n]]);

    for( int l=0;l<N+1;++l ) {
      Real z = (grid[l]+1)/2;
      theta.get(0).at(l) = buf.X().at(0);
      theta.get(0).at(l) *= z*(2-z);
    }
  }

  circle.setLz(Lz);

  bool secondorder=true;

  /* pressure solver */
  extrudedConsistentPressureEvaluator consistent3D(circle,SP.Dx(),GLL2G,weightGL,SP.weightX(),
      grid,gridGL,preconditioned,rank,size,
      SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM);
  extrudedFEMConsistentPressurePreconditioner ppre(circle,SP.weightX(),weightGL,grid,
      gridGL,rank,size);

  /* velocity solver */
  Real nu;
  if( secondorder )
    nu = Real(3)/(2*Dt*Pr);
  else
    nu = Real(1)/(Dt*Pr);

  SEMLaplacianEvaluator velocity2D(circle,SP.Dx(),SP.weightX(),NULL,0,1,
      SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
  extrudedLaplacianEvaluator velocity3D1(velocity2D,SP3,nu,&SP2.Ax(),preconditioned,
      rank,size,SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM);
  extrudedLaplacianEvaluator velocity3D2(velocity2D,SP3,nu,&SP2.Ax(),preconditioned,
      rank,size,SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM);
  extrudedLaplacianEvaluator velocity3D3(velocity2D,SP,nu,&SP2.Ax(),preconditioned,
      rank,size,SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
  extrudedFEMLaplacianPreconditioner upre1(circle,SP.weightX(),grid,nu,
      rank,size,SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM);
  extrudedFEMLaplacianPreconditioner upre3(circle,SP.weightX(),grid,nu,
      rank,size,SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);

  /* temperature solver */
  SEMLaplacianEvaluator temperature2D(circle,SP.Dx(),SP.weightX(),NULL,0,1,
      SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN);
  extrudedLaplacianEvaluator temperature3D(temperature2D,SP3,nu*Pr,&SP2.Ax(),preconditioned,
      rank,size,SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM);
  extrudedFEMLaplacianPreconditioner tpre(circle,SP.weightX(),grid,nu*Pr,
      rank,size,SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM);
  temperature3D.m_name = "temperature";

  if( n == 0 ) {
    velocity3D1.m_nu = velocity3D2.m_nu = velocity3D3.m_nu = Real(1)/(Pr*Dt);
    temperature3D.m_nu = Real(1)/Dt;
  }

  /* convection solver */
  convectionEvaluator conv(circle,SP.Dx(),SP.weightX(),rank,size);
  if( n == 0 || !secondorder )
    conv.m_order = 0;
  else
    conv.m_order = 1;

  MPIDotter pdotter(rank,size);
  MPIGeometryDotter<basics::geometryStack> udotter(circle,desc[rank].elements,rank,size);

  if( rank == 0 )
    cout << "time step:\t\t"  	 << Dt << endl
      <<	"rayleigh number:\t" << Ra << endl
      <<	"maragoni number:\t" << Ma << endl
      <<	"prandtl number:\t\t"<< Pr << endl
      <<	"final time:\t\t" 	 <<  T << endl
      <<	"Lz:\t\t\t"  		 << Lz << endl
      <<	"lumped:\t\t\t"		 << (lumped?"yes":"no") << endl
      <<	"symmetry:\t\t"	 	 << (symmetry>0?(symmetry==2?"octant":"quadrant"):"none") << endl
      <<	"preconditioned:\t\t"<< (preconditioned?"yes":"no") << endl;

  /* setup keyboard */
  utilities::set_raw_tty();
  int vel = utilities::g_profiler.add("velocity solve");
  int pre = utilities::g_profiler.add("pressure solve");
  int tmp = utilities::g_profiler.add("temperature solve");
  int cnv = utilities::g_profiler.add("convection problems");
  for( n;n<(T-t0)/Dt;++n ) {
    if( rank == 0 )
      std::cout << "Time: " << t0+(n+1)*Dt << std::endl;

    /* check keyboard */
    if( utilities::keypress() == 's' )
      saveState(u,p,theta,n,Dt,Ra,Ma,Pr,Lz,rank);
    if( utilities::keypress() == 'q' )
      break;

    if( n % 100 == 0 ) // save every 100th timestep
      saveState(u,p,theta,n,Dt,Ra,Ma,Pr,Lz,rank);

    /* set up references */
    basics::Field3<basics::matricesStack>& unm1	= u.get(n,-1);
    basics::Field3<basics::matricesStack>& un  	= u.get(n, 0);
    basics::Field3<basics::matricesStack>& unp1	= u.get(n, 1);
    basics::matricesStack& pnm1  				= p.get(n,-1);
    basics::matricesStack& pn  					= p.get(n, 0);
    basics::matricesStack& pnp1 				= p.get(n, 1);
    basics::matricesStack& tnm1					= theta.get(n,-1);
    basics::matricesStack& tn  					= theta.get(n, 0);
    basics::matricesStack& tnp1 				= theta.get(n, 1);

    /* velocity convection problems */
    if( n > 10 )
      utilities::g_profiler.start(cnv);
    conv.m_bc = velocity3D1.m_bc;
    conv.solve(unp1.X(),un.X(),u,n,Dt,Dt,0,Dt);
    conv.solve(unp1.Y(),un.Y(),u,n,Dt,Dt,0,Dt);
    conv.m_bc = velocity3D3.m_bc;
    conv.solve(unp1.Z(),un.Z(),u,n,Dt,Dt,0,Dt);
    if( conv.m_order > 0 ) {
      conv.m_bc = velocity3D1.m_bc;
      conv.solve(buf.X(),unm1.X(),u,n,Dt,Dt,-Dt,Dt);
      conv.solve(buf.Y(),unm1.Y(),u,n,Dt,Dt,-Dt,Dt);
      conv.m_bc = velocity3D3.m_bc;
      conv.solve(buf.Z(),unm1.Z(),u,n,Dt,Dt,-Dt,Dt);
    }
    utilities::g_profiler.pause(cnv);

    /* BDF */
    if( conv.m_order > 0 ) {
      unp1 *= Real(2)/(Dt*Pr);
      unp1.axpy(Real(-1)/(2*Dt*Pr),buf);
    }
    else
      unp1 *= Real(1)/(Dt*Pr);

    /* temperature */
    unp1.Z().axpy(Ra,tn);
    circle.mass(unp1,desc[rank].elements);

    /* upper tangential boundary conditions */
    int max=desc[rank].elements.size();
#pragma omp parallel for schedule(static)
    for( int i=0;i<max;++i ) {
      int elem=desc[rank].elements[i];
      basics::multTranspose(buf.X().at(0)[i],SP.Dx(),tn.at(N)[i],'N','N');
      basics::multTranspose(buf.Y().at(0)[i],tn.at(N)[i],SP.Dx(),'N','T');

      /* construct dtheta/dx */
      basics::multPointwise(buf.X().at(1)[i],buf.X().at(0)[i],
          circle[elem].getGH().getGeometryDerivatives()[3]);
      basics::multPointwise(buf.Y().at(1)[i],buf.Y().at(0)[i],
          circle[elem].getGH().getGeometryDerivatives()[2]);
      buf.X().at(1)[i] -= buf.Y().at(1)[i];

      /* construct dtheta/dy */
      basics::multPointwise(buf.Y().at(1)[i],buf.Y().at(0)[i],
          circle[elem].getGH().getGeometryDerivatives()[0]);
      basics::multPointwise(buf.Y().at(0)[i],buf.X().at(0)[i],
          circle[elem].getGH().getGeometryDerivatives()[1]);
      buf.Y().at(1)[i] -= buf.Y().at(0)[i];

      massReference(buf.X().at(1)[i],SP.weightX());
      massReference(buf.Y().at(1)[i],SP.weightX());
    }
    unp1.X().at(N).axpy(-Ma,buf.X().at(1));
    unp1.Y().at(N).axpy(-Ma,buf.Y().at(1));

    /* delta */
    velocity3D1.evaluate(buf,un,false,false);
    unp1 -= buf;
    /* extrapolate pressure */
    if( conv.m_order == 0 )
      utilities::extrapolate(pnp1,p,n,utilities::ZEROTH_ORDER);
    else
      utilities::extrapolate(pnp1,p,n,utilities::FIRST_ORDER);
    consistent3D.m_gradient.evaluate(buf,pnp1);
    unp1 += buf;
    /* solve helmholtz problems */
    mask(unp1,circle,rank,size,velocity3D1.m_bc);

    buf = unp1;
    dssum(unp1,circle,rank,size);
    utilities::g_profiler.start(vel);

    velocity3D1.m_name = "velocity X";
    velocity3D1.m_nosum = upre1.m_nosum = &buf.X();
    cout << "iters velocity X " << utilities::CGSolver::solve(unp1.X(),velocity3D1,
        upre1,udotter,1.e-10) << endl;
    velocity3D2.m_nosum = upre1.m_nosum = &buf.Y();
    cout << "iters velocity Y " << utilities::CGSolver::solve(unp1.Y(),velocity3D2,
        upre1,udotter,1.e-10) << endl;
    velocity3D3.m_nosum = upre3.m_nosum = &buf.Z();
    cout << "iters velocity Z " << utilities::CGSolver::solve(unp1.Z(),velocity3D3,
        upre3,udotter,1.e-10) << endl;
    unp1 += un;
    utilities::g_profiler.pause(vel);

    /* solve pressure problem */
    consistent3D.m_divergence.evaluate(pt,unp1);
    pt *= -velocity3D1.m_nu;
    utilities::g_profiler.start(pre);
    cout << "iterations pressure " << utilities::CGSolver::solve(pt,consistent3D,ppre,
        pdotter,1.e-10) << endl;
    utilities::g_profiler.pause(pre);
    pnp1 += pt;

    /* update velocity */
    consistent3D.m_gradient.evaluate(buf,pt);
    dssum(buf,circle,rank,size);
    circle.invMass(buf,desc[rank].elements);
    mask(buf,circle,rank,size,velocity3D1.m_bc);
    unp1.axpy(Real(1)/velocity3D1.m_nu,buf);

    /* solve temperature convection problem */
    utilities::g_profiler.start(cnv);
    conv.m_bc = temperature3D.m_bc;
    conv.solve(tnp1,tn,u,n,Dt,Dt,0,Dt);
    if( conv.m_order > 0 )
      conv.solve(buf.X(),tnm1,u,n,Dt,Dt,-Dt,Dt);
    utilities::g_profiler.pause(cnv);
    if( conv.m_order > 0) {
      tnp1 *= Real(2)/Dt;
      tnp1.axpy(Real(-1)/(2*Dt),buf.X());
    }
    else
      tnp1 *= Real(1)/Dt;

    tnp1 += unp1.Z();
    circle.mass(tnp1,desc[rank].elements);

    /* delta */
    temperature3D.evaluate(buf.X(),tn,false,false);
    tnp1 -= buf.X();
    mask(tnp1,circle,rank,size,temperature3D.m_bc);
    buf.X() = tnp1;
    dssum(tnp1,circle,rank,size);

    /* solve the elliptic temperature problem */
    utilities::g_profiler.start(tmp);
    temperature3D.m_nosum = tpre.m_nosum = &buf.X();
    cout << "iters temperature " << utilities::CGSolver::solve(tnp1,temperature3D,
        tpre,udotter,1.e-10) << endl;
    utilities::g_profiler.pause(tmp);
    tnp1 += tn;

    if( secondorder ) {
      velocity3D1.m_nu = velocity3D2.m_nu = velocity3D3.m_nu = Real(3)/(2*Dt*Pr);
      temperature3D.m_nu = Real(3)/(2*Dt);
      conv.m_order = 1;
    }
  }

  cout << utilities::g_profiler.report();

  /* restore keyboard */
  utilities::set_normal_tty();

  saveState(u,p,theta,n,Dt,Ra,Ma,Pr,Lz,rank);
#ifdef HAS_MPI
  MPI_Finalize();
#endif

  return 0; // AY-OH-KAY
}

