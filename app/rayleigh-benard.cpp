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
#include "bm/consistent.h"
#include "bm/mass.h"
#include "bm/sem.h"
#include "bm/convection.h"

class initialTemperature : public basics::function3D {
  public:
    initialTemperature(const Real theta0) :
      m_theta0(theta0)
  {
    srand(time(NULL));
  }

    Real val(Real x, Real y, Real z, Real t) const
    {
      //            return -1;
      return 1-(z+1)/2;
      //            return (m_theta0*rand())*(1-(z+1)/2)/(RAND_MAX);
    }
  protected:
    Real m_theta0;
};

void loadState(utilities::ringBuffer<basics::Field3<basics::matricesStack> >& u,
    utilities::ringBuffer<basics::matricesStack>& p,
    utilities::ringBuffer<basics::matricesStack>& theta,
    int& n, Real& Dt, Real& Ra, const string& file)
{
  HDF5::HDF5Reader reader(file);
  basics::Vector time("time info",3);
  reader.read(time,"time info");
  n = time[0];
  Dt = time[1];
  Ra = time[2];
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

void saveState(const utilities::ringBuffer<basics::Field3<basics::matricesStack> >& u,
    const utilities::ringBuffer<basics::matricesStack>& p,
    const utilities::ringBuffer<basics::matricesStack>& theta, int n, Real Dt, Real Ra)
{
  return;
  stringstream str;
  str << "state-rb-" << n << ".hdf5";

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
  basics::Vector time("time info",3);
  time[0] = n;
  time[1] = Dt;
  time[2] = Ra;
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
  if( argc > 2 && strstr(argv[2],"hdf5"))
    load = true;
  if( argc > 3 )
    Dt = Real(atoi(argv[2]))/atoi(argv[3]);
  if( argc > 5 )
    T = Real(atoi(argv[4]))/atoi(argv[5]);
  Real Ra = 2000;
  if (argc > 7 )
    Ra = Real(atoi(argv[6]))/atoi(argv[7]);
  bool preconditioned=true;
  if( argc > 8 )
    preconditioned = atoi(argv[8])==1?false:true;

  legendreLegendreW::poissonSolver SP(N,N,0,legendreLegendreW::poissonSolver::Homogenous,true,-1);
  legendreLegendreW::poissonSolver SP2(N,4,0,legendreLegendreW::poissonSolver::Nonhomogenous,true,-1);

  basics::Vector gridGL("GL grid",N-1);
  basics::Vector weightGL("GL weight",N-1);
  utilities::GLL::GaussLegendreGrid(gridGL);
  utilities::GLL::GaussLegendreWeights(weightGL,gridGL);
  basics::Matrix GLL2G = utilities::GLL::interpolationMatrix(SP.gridX(),gridGL);

  bigCircleGeometry circle(N,N,SP.gridX(),SP.weightX());

  /* pressure solver */
  extrudedConsistentPressureEvaluator consistent3D(circle,SP.Dx(),GLL2G,weightGL,SP.weightX(),
      SP.gridX(),gridGL,preconditioned,rank,size);

  //    extrudedFEMConsistentPressurePreconditioner eval4(circle,SP.weightX(),weightGL,SP.gridX(),gridGL);

  /* velocity solver */
  SEMLaplacianEvaluator velocity2D(circle,SP.Dx(),SP.weightX(),NULL,0,1,
      SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
  extrudedLaplacianEvaluator velocity3D(velocity2D,SP,Real(3)/(2*Dt),&SP2.Ax(),preconditioned,
      rank,size,SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
  velocity3D.m_name = "velocity";

  /* temperature solver */
  SEMLaplacianEvaluator temperature2D(circle,SP.Dx(),SP.weightX(),NULL,0,1,
      SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN);
  extrudedLaplacianEvaluator temperature3D(temperature2D,SP,Real(3)/(2*Dt),&SP2.Ax(),preconditioned,
      rank,size,SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
  temperature3D.m_name = "temperature";

  velocity3D.m_nu = temperature3D.m_nu = Real(1)/Dt;

  /* convection solver */
  convectionEvaluator conv(circle,SP.Dx(),SP.weightX(),rank,size);

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

  basics::Matrix iters("iters",N-1,T/Dt+1);
  basics::Matrix counts("counts",size,T/Dt+1);
  basics::Matrix planes("planes",N-1,T/Dt+1);

  int n=0;
  Real t0=0;
  if( load )
    loadState(u,p,theta,n,Dt,Ra,argv[2]);
  else {
    /* set up initial conditions */
    Real theta0=1.e-3;

    initialTemperature temp(theta0);
    circle.evaluate(theta.get(0),SP.gridX(),temp,t0,desc[rank].elements);
  }

  if( rank == 0 ) {
    cout << "time step:\t\t"  	 << Dt << endl
      <<	"rayleigh number:\t" << Ra << endl
      <<	"final time:\t\t" 	 <<  T << endl
      <<	"preconditioned:\t\t"<< (preconditioned?"yes":"no") << endl;
  }

  int elipticu    = utilities::g_profiler.add("velocity solve");
  int convectiveu = utilities::g_profiler.add("velocity convection");
  int elipticp    = utilities::g_profiler.add("pressure solve");
  int eliptict    = utilities::g_profiler.add("temperature solve");
  int convectivet = utilities::g_profiler.add("temperature convection");

  /* setup keyboard */
  utilities::set_raw_tty();
  for( n;n<(T-t0)/Dt;++n ) {
    if( rank == 0 )
      std::cout << "Time: " << t0+(n+1)*Dt << std::endl;

    /* check keyboard */
    if( utilities::keypress() == 's' )
      saveState(u,p,theta,n,Dt,Ra);
    if( utilities::keypress() == 'q' )
      break;

    if( n % 100 == 0 && n ) // save every 100th timestep
      saveState(u,p,theta,n,Dt,Ra);

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
      utilities::g_profiler.start(convectiveu);
    conv.solve(unp1,un,u,n,Dt,Dt,0,Dt);
    if( n > 0 )
      conv.solve(buf,unm1,u,n,Dt,Dt,-Dt,Dt);
    if( n > 10 )
      utilities::g_profiler.pause(convectiveu);

    /* BDF */
    if( n > 0 ) {
      unp1 *= Real(2)/Dt;
      unp1.axpy(Real(-1)/(2*Dt),buf);
    }
    else
      unp1 *= Real(1)/Dt;

    unp1.Z().axpy(Ra,tn);
    circle.mass(unp1,desc[rank].elements);

    /* delta */
    velocity3D.evaluate(buf,un,false,false);
    unp1 -= buf;

    /* extrapolate pressure */
    utilities::extrapolate(pnp1,p,n,utilities::FIRST_ORDER);
    consistent3D.m_gradient.evaluate(buf,pnp1);
    unp1 += buf;

    /* solve helmholtz problems */
    mask(unp1,circle,rank,size,velocity3D.m_bc);
    buf = unp1;
    dssum(unp1,circle,rank,size);
    if( n > 10 )
      utilities::g_profiler.start(elipticu);
    //                for( int i=0;i<N-1;++i )
    //                        iters[n][i] = velocity3D.m_stat[0]->get()[i];
    //                for( int i=0;i<velocity3D.m_counts.first.size();++i )
    //                        counts[n][i] = velocity3D.m_counts.first[i];
    //                int l=0;
    //                for( int k=0;k<velocity3D.m_counts.second.size();++k ) {
    //                        for( int i=0;i<velocity3D.m_counts.second[k].size();++i )
    //                                planes[n][l++] = velocity3D.m_counts.second[k][i];
    //                }
    velocity3D.solve(unp1,buf);
    if( n > 10 )
      utilities::g_profiler.pause(elipticu);
    unp1 += un;

    /* solve pressure problem */
    consistent3D.m_divergence.evaluate(pt,unp1);
    pt *= -velocity3D.m_nu;
    if( n > 10 )
      utilities::g_profiler.start(elipticp);
    //        for( int i=0;i<N-1;++i )
    //            iters[n][i] = consistent3D.m_stat.get()[i];
    //        vector<int> ncounts = consistent3D.m_stat.getCounts(size);
    //        for( int i=0;i<size;++i )
    //            counts[n][i] = ncounts[i];
    //        for( int i=0;i<size;++i )
    //            counts[n][i] = consistent3D.m_counts.first[i];
    //        int l=0;
    //        for( int k=0;k<consistent3D.m_counts.second.size();++k ) {
    //            for( int i=0;i<consistent3D.m_counts.second[k].size();++i )
    //                planes[n][l++] = consistent3D.m_counts.second[k][i];
    //        }
    consistent3D.solve(pt);
    if( n > 10 )
      utilities::g_profiler.pause(elipticp);
    //        if( preconditioned )
    //            cout << "iterations pressure " << utilities::CGSolver::solve(pt,consistent3D,eval4,pt.getDotter(),1.e-10) << endl;
    //        else
    //            cout << "iterations pressure " << utilities::CGSolver::solve(pt,consistent3D,pt.getDotter(),1.e-10) << endl;
    pnp1 += pt;
    //        pnp1 -= pnp1.sum()/pnp1.length();
    //        cout << "mean " << pnp1.sum()/pnp1.length() << endl;

    /* update velocity */
    consistent3D.m_gradient.evaluate(buf,pt);
    dssum(buf,circle,rank,size);
    circle.invMass(buf,desc[rank].elements);
    mask(buf,circle,rank,size,velocity3D.m_bc);
    unp1.axpy(Real(1)/velocity3D.m_nu,buf);

    /* solve temperature convection problem */
    if( n > 10 )
      utilities::g_profiler.start(convectivet);
    conv.solve(tnp1,tn,u,n,Dt,Dt,0,Dt);
    if( n > 0 )
      conv.solve(buf.X(),tnm1,u,n,Dt,Dt,-Dt,Dt);
    if( n > 10 )
      utilities::g_profiler.pause(convectivet);

    if( n > 0 ) {
      tnp1 *= Real(2)/Dt;
      tnp1.axpy(Real(-1)/(2*Dt),buf.X());
    }
    else
      tnp1 *= Real(1)/Dt;

    circle.mass(tnp1,desc[rank].elements);

    /* delta */
    buf.Y() = tn;
    buf.Y().at(0) = 1;
    temperature3D.evaluate(buf.X(),buf.Y(),false,false);
    tnp1 -= buf.X();
    tnp1.at(0).clear();
    tnp1.at(N).clear();
    buf.X() = tnp1;
    dssum(tnp1,circle,rank,size);

    /* solve the elliptic temperature problem */
    if( n > 10 )
      utilities::g_profiler.start(eliptict);
    temperature3D.solve(tnp1,buf.X());
    if( n > 10 )
      utilities::g_profiler.pause(eliptict);
    tnp1 += tn;
    tnp1.at(0) = 1;

    velocity3D.m_nu = temperature3D.m_nu = Real(3)/(2*Dt);
  }
  printf("store counts\n");
  n = iters.cols()-1;
  for( int i=0;i<N-1;++i )
    iters[n][i] = velocity3D.m_stat[0]->get()[i];
  for( int i=0;i<size;++i )
    counts[n][i] = velocity3D.m_counts.first[i];
  int l=0;
  for( int k=0;k<velocity3D.m_counts.second.size();++k ) {
    for( int i=0;i<velocity3D.m_counts.second[k].size();++i )
      planes[n][l++] = velocity3D.m_counts.second[k][i];
  }
  //    for( int i=0;i<size;++i )
  //        counts[n][i] = consistent3D.m_counts.first[i];
  //    int l=0;
  //    for( int k=0;k<consistent3D.m_counts.second.size();++k ) {
  //        for( int i=0;i<consistent3D.m_counts.second[k].size();++i )
  //            planes[n][l++] = consistent3D.m_counts.second[k][i];
  //    }

  printf("done\n");
  if( rank == 0 )
    cout << utilities::g_profiler.report();

  if( rank == 0 ) {
    char temp[100];
    sprintf(temp,"iters%i",size);
    iters.save(temp);
    sprintf(temp,"counts%i",size);
    counts.save(temp);
    sprintf(temp,"planes%i",size);
    planes.save(temp);
  }

  /* restore keyboard */
  utilities::set_normal_tty();

  //        saveState(u,p,theta,n,Dt,Ra);

#ifdef HAS_MPI
  MPI_Finalize();
#endif
  return 0; // AY-OH-KAY
}

