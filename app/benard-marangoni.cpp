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

class testTemperature : public basics::function2D {
  public:
    Real val(Real x, Real y, Real t) const
    {
#undef RADIUS
      Real RADIUS = 2;
      Real r = sqrt(x*x+y*y);

      return .5*y*y;

      if( abs(r) < 1.e-14 )
        return 0;

      //return Real(-1)/6*(-3*RADIUS+2*r)*y*(3*x*x-y*y)/r;
    }
};

class simplyX : public basics::function2D {
  public:
    Real val(Real x, Real y, Real t) const
    {
      return x;
    }
};

class simplyY : public basics::function2D {
  public:
    Real val(Real x, Real y, Real t) const
    {
      return y;
    }
};

class simplyZ : public basics::function3D {
  public:
    Real val(Real x, Real y, Real z, Real t) const
    {
      return z;
    }
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

void flipx(basics::Matrix& A, const basics::Matrix& B)
{
  int N=A.rows()-1;

  for( int j=0;j<A.cols();++j )
    for( int k=0;k<A.rows();++k )
      A[j][k] = B[j][N-k];
}

void flipy(basics::Matrix& A, const basics::Matrix& B)
{
  int N=A.rows()-1;

  for( int j=0;j<A.cols();++j )
    for( int k=0;k<A.rows();++k )
      A[j][k] = B[N-j][k];
}

//void totalflip(basics::Matrix& A, const basics::Matrix& B)
//{
//    int N=A.rows()-1;

//    for( int j=0;j<A.cols();++j )
//        for( int k=0;k<A.rows();++k )
//            A[j][k] = B[k][N-j];
//}

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
  //    quadraticGeometry circle(N,N,grid,SP.weightX());

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

    /* bigcircle connectivity */
    basics::matrixStack& buff = theta.get(0).at(0);

    if( symmetry == 2 ) {
      if( circle.size() == 16 ) {
        /* quadratic connectivity - fully symmetrized */
        buff[0] += theta.get(0).at(0)[0].transposed();
        buff[0] *= Real(1)/2;
        buff[5] += buff[5].transposed();
        buff[5] *= Real(1)/2;
        buff[1] = buff[4].transposed();
      }
      if( circle.size() == 48 ) {
        /* symmetrices middle */
        buff[2] += buff[2].transposed();
        buff[2] *= Real(1)/2;
        buff[1] += buff[1].transposed();
        buff[1] *= Real(1)/2;
        buff[3] = buff[0].transposed();

        /* symmetrices others */
        buff[46] = buff[18].transposed();
        buff[47] = buff[16].transposed();
        buff[44] = buff[19].transposed();
        buff[45] = buff[17].transposed();
      }
    }

    if( symmetry > 0 ) {
      if( circle.size() == 16 ) {
        /* quadratic connectivity */
        flipx(buff[2],buff[1]);
        flipx(buff[3],buff[0]);
        flipx(buff[6],buff[5]);
        flipx(buff[7],buff[4]);
        flipy(buff[8],buff[4]);
        flipy(buff[12],buff[0]);
        flipy(buff[9],buff[5]);
        flipy(buff[13],buff[1]);
        flipy(buff[10],buff[6]);
        flipy(buff[14],buff[2]);
        flipy(buff[11],buff[7]);
        flipy(buff[15],buff[3]);
      }
      if( circle.size() == 48 ) {

        /* quadrant flippling - bigcircle */
        flipy(buff[14],buff[ 0]);
        flipy(buff[15],buff[ 1]);
        flipy(buff[43],buff[45]);
        flipy(buff[42],buff[44]);
        flipy(buff[12],buff[ 2]);
        flipy(buff[13],buff[ 3]);
        flipy(buff[40],buff[46]);
        flipy(buff[41],buff[47]);
        flipy(buff[38],buff[16]);
        flipy(buff[39],buff[17]);
        flipy(buff[36],buff[18]);
        flipy(buff[37],buff[19]);

        flipx(buff[22],buff[19]);
        flipx(buff[20],buff[17]);
        flipx(buff[ 6],buff[ 3]);
        flipx(buff[ 4],buff[ 1]);
        flipx(buff[10],buff[15]);
        flipx(buff[ 8],buff[13]);
        flipx(buff[34],buff[39]);
        flipx(buff[32],buff[37]);

        flipx(buff[23],buff[18]);
        flipx(buff[21],buff[16]);
        flipx(buff[ 7],buff[ 2]);
        flipx(buff[ 5],buff[ 0]);
        flipx(buff[11],buff[14]);
        flipx(buff[ 9],buff[12]);
        flipx(buff[35],buff[38]);
        flipx(buff[33],buff[36]);

        //                flipx(buff[27],buff[46]);
        //                flipx(buff[26],buff[47]);
        //                flipx(buff[24],buff[45]);
        //                flipx(buff[25],buff[44]);
        flipx(buff[30],buff[43]);
        flipx(buff[31],buff[42]);
        flipx(buff[28],buff[41]);
        flipx(buff[29],buff[40]);
      }
      if( circle.size() == 192 ) {
        /* 1<->2 */
        for( int i=0;i<4;++i )
          flipy(buff[60+i],buff[i]);
        for( int i=0;i<4;++i )
          flipy(buff[56+i],buff[4+i]);
        for( int i=0;i<4;++i )
          flipy(buff[52+i],buff[8+i]);
        for( int i=0;i<4;++i )
          flipy(buff[48+i],buff[12+i]);
        for( int i=0;i<4;++i )
          flipy(buff[156+i],buff[64+i]);
        for( int i=0;i<4;++i )
          flipy(buff[152+i],buff[68+i]);
        for( int i=0;i<4;++i )
          flipy(buff[148+i],buff[72+i]);
        for( int i=0;i<4;++i )
          flipy(buff[144+i],buff[76+i]);
        for( int i=0;i<4;++i )
          flipy(buff[172+i],buff[176+i]);
        for( int i=0;i<4;++i )
          flipy(buff[168+i],buff[180+i]);
        for( int i=0;i<4;++i )
          flipy(buff[164+i],buff[184+i]);
        for( int i=0;i<4;++i )
          flipy(buff[160+i],buff[188+i]);

        /* 1 <-> 3 */
        for( int i=0;i<4;++i )
          flipx(buff[16+4*i],buff[3+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[17+4*i],buff[2+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[18+4*i],buff[1+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[19+4*i],buff[0+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[80+4*i],buff[67+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[81+4*i],buff[66+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[82+4*i],buff[65+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[83+4*i],buff[64+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[96+4*i],buff[179+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[97+4*i],buff[178+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[98+4*i],buff[177+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[99+4*i],buff[176+4*i]);

        /* 2 <-> 4 */
        for( int i=0;i<4;++i )
          flipx(buff[32+4*i],buff[51+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[33+4*i],buff[50+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[34+4*i],buff[49+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[35+4*i],buff[48+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[112+4*i],buff[163+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[113+4*i],buff[162+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[114+4*i],buff[161+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[115+4*i],buff[160+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[128+4*i],buff[147+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[129+4*i],buff[146+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[130+4*i],buff[145+4*i]);
        for( int i=0;i<4;++i )
          flipx(buff[131+4*i],buff[144+4*i]);
      }
    }

    dssum(theta.get(0),circle,rank,size);
    for( int n=0;n<desc[rank].elements.size();++n )
      basics::multPointwise(buf.X().at(0)[n],theta.get(0).at(0)[n],(*circle.m_mult)[desc[rank].elements[n]]);

    for( int l=0;l<N+1;++l ) {
      Real z = (grid[l]+1)/2;
      theta.get(0).at(l) = buf.X().at(0);
      theta.get(0).at(l) *= z*(2-z);
    }
  }

  //    basics::matricesStack* stack = utilities::g_manager.aquireMatricesStack("init",circle);
  //    HDF5::HDF5Reader reader("init.hdf5");
  //    reader.read("temperature 1",*stack);
  //    for( int i=0;i<desc[rank].elements.size();++i )
  //        theta.get(0)[i] = (*stack)[desc[rank].elements[i]];

  //    t0 = n*Dt;
  //    addn = n;
  //    n = 0;
  //    Dt = Real(1)/300;
  //    Ma = 85;
  //    Ra = 19;
  //    Dt = Real(1)/500;
  //    saveState(u,p,theta,n,Dt,Ra,Ma,Pr,Lz,rank);

  circle.setLz(Lz);

  bool secondorder=true;

  /* pressure solver */
  extrudedConsistentPressureEvaluator consistent3D(circle,SP.Dx(),GLL2G,weightGL,SP.weightX(),
      grid,gridGL,preconditioned,rank,size,
      SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM);

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

  /* temperature solver */
  SEMLaplacianEvaluator temperature2D(circle,SP.Dx(),SP.weightX(),NULL,0,1,
      SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN);
  extrudedLaplacianEvaluator temperature3D(temperature2D,SP3,nu*Pr,&SP2.Ax(),preconditioned,
      rank,size,SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM);
  temperature3D.m_name = "temperature";

  if( n == 0 ) {
    velocity3D1.m_nu = velocity3D2.m_nu = velocity3D3.m_nu = Real(1)/(Pr*Dt);
    temperature3D.m_nu = Real(1)/Dt;
  }

  //    basics::Matrix iters("iters",3*N-1,T/Dt+1);
  //    basics::Matrix counts("counts",size,T/Dt+1);
  //    basics::Matrix planes("planes",3*N-1,T/Dt+1);

  //    basics::Matrix iters1("iters",N,T/Dt+1);
  //    basics::Matrix counts1("counts",size,T/Dt+1);
  //    basics::Matrix planes1("planes",N,T/Dt+1);
  //    basics::Matrix iters2("iters",N,T/Dt+1);
  //    basics::Matrix counts2("counts",size,T/Dt+1);
  //    basics::Matrix planes2("planes",N,T/Dt+1);
  //    basics::Matrix iters3("iters",N-1,T/Dt+1);
  //    basics::Matrix counts3("counts",size,T/Dt+1);
  //    basics::Matrix planes3("planes",N-1,T/Dt+1);
  //    basics::Matrix itersc("iters",N-1,T/Dt+1);
  //    basics::Matrix countsc("counts",size,T/Dt+1);
  //    basics::Matrix planesc("planes",N-1,T/Dt+1);
  //    basics::Matrix iterst("iters",N,T/Dt+1);
  //    basics::Matrix countst("counts",size,T/Dt+1);
  //    basics::Matrix planest("planes",N,T/Dt+1);

  /* convection solver */
  convectionEvaluator conv(circle,SP.Dx(),SP.weightX(),rank,size);
  if( n == 0 || !secondorder )
    conv.m_order = 0;
  else
    conv.m_order = 1;

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

    if( lumped )
      velocity3D1.solve(unp1,buf,velocity3D1,velocity3D2,velocity3D3);
    else {
      velocity3D1.m_name = "velocity X";
      velocity3D1.solve(unp1.X(),buf.X());
      velocity3D2.m_name = "velocity Y";
      velocity3D2.solve(unp1.Y(),buf.Y());
      velocity3D3.m_name = "velocity Z";
      velocity3D3.solve(unp1.Z(),buf.Z());
    }
    unp1 += un;
    utilities::g_profiler.pause(vel);

    /* solve pressure problem */
    consistent3D.m_divergence.evaluate(pt,unp1);
    pt *= -velocity3D1.m_nu;
    utilities::g_profiler.start(pre);
    consistent3D.solve(pt);
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
    temperature3D.solve(tnp1,buf.X());
    utilities::g_profiler.pause(tmp);
    tnp1 += tn;

    if( secondorder ) {
      velocity3D1.m_nu = velocity3D2.m_nu = velocity3D3.m_nu = Real(3)/(2*Dt*Pr);
      temperature3D.m_nu = Real(3)/(2*Dt);
      conv.m_order = 1;
    }
    //        if( lumped ) {
    //            pair< vector<int>,vector< vector<int> > > count = velocity3D1.m_stat3->getCounts2(size);
    //            for( int i=0;i<velocity3D1.m_stat3->m_entries;++i )
    //                iters[n][i] = velocity3D1.m_stat3->get()[i];
    //            for( int i=0;i<size;++i )
    //                counts[n][i] = velocity3D1.m_counts3.first[i];
    //            int l=0;
    //            for( int k=0;k<velocity3D1.m_counts3.second.size();++k ) {
    //                for( int i=0;i<velocity3D1.m_counts3.second[k].size();++i )
    //                    planes[n][l++] = velocity3D1.m_counts3.second[k][i];
    //            }
    //        } else {
    //            pair< vector<int>,vector< vector<int> > > count 
    //                = velocity3D1.m_stat[0]->getCounts2(size);
    //            for( int i=0;i<N;++i )
    //                iters1[n][i] = velocity3D1.m_stat[0]->get()[i];
    //            for( int i=0;i<size;++i )
    //                counts1[n][i] = velocity3D1.m_counts.first[i];
    //            int l=0;
    //            for( int k=0;k<velocity3D1.m_counts.second.size();++k ) {
    //                for( int i=0;i<velocity3D1.m_counts.second[k].size();++i )
    //                    planes1[n][l++] = velocity3D1.m_counts.second[k][i];
    //            }
    //            count = velocity3D2.m_stat[0]->getCounts2(size);
    //            for( int i=0;i<N;++i )
    //                iters2[n][i] = velocity3D2.m_stat[0]->get()[i];
    //            for( int i=0;i<size;++i )
    //                counts2[n][i] = velocity3D2.m_counts.first[i];
    //            l=0;
    //            for( int k=0;k<velocity3D2.m_counts.second.size();++k ) {
    //                for( int i=0;i<velocity3D2.m_counts.second[k].size();++i )
    //                    planes2[n][l++] = velocity3D2.m_counts.second[k][i];
    //            }
    //            count = velocity3D3.m_stat[0]->getCounts2(size);
    //            for( int i=0;i<N-1;++i )
    //                iters3[n][i] = velocity3D3.m_stat[0]->get()[i];
    //            for( int i=0;i<size;++i )
    //                counts3[n][i] = velocity3D3.m_counts.first[i];
    //            l=0;
    //            for( int k=0;k<velocity3D3.m_counts.second.size();++k ) {
    //                for( int i=0;i<velocity3D3.m_counts.second[k].size();++i )
    //                    planes3[n][l++] = velocity3D3.m_counts.second[k][i];
    //            }
    //        }
    //        pair< vector<int>,vector< vector<int> > > count 
    //            = consistent3D.m_stat.getCounts2(size);
    //        for( int i=0;i<N-1;++i )
    //            itersc[n][i] = consistent3D.m_stat.get()[i];
    //        for( int i=0;i<size;++i )
    //            countsc[n][i] = consistent3D.m_counts.first[i];
    //        int l=0;
    //        for( int k=0;k<consistent3D.m_counts.second.size();++k ) {
    //            for( int i=0;i<consistent3D.m_counts.second[k].size();++i )
    //                planesc[n][l++] = consistent3D.m_counts.second[k][i];
    //        }
    //        count = temperature3D.m_stat[0]->getCounts2(size);
    //        for( int i=0;i<N;++i )
    //            iterst[n][i] = temperature3D.m_stat[0]->get()[i];
    //        for( int i=0;i<size;++i )
    //            countst[n][i] = temperature3D.m_counts.first[i];
    //        l=0;
    //        for( int k=0;k<temperature3D.m_counts.second.size();++k ) {
    //            for( int i=0;i<temperature3D.m_counts.second[k].size();++i )
    //                planest[n][l++] = temperature3D.m_counts.second[k][i];
    //        }
  }

  cout << utilities::g_profiler.report();

  //    if( rank == 0 ) {
  //        char temp[100];
  //        if( lumped ) {
  //            sprintf(temp,"iters-lumped-%i",size);
  //            iters.save(temp);
  //            sprintf(temp,"counts-lumped-%i",size);
  //            counts.save(temp);
  //            sprintf(temp,"planes-lumped-%i",size);
  //            planes.save(temp);
  //        } else {
  //            sprintf(temp,"itersx%i",size);
  //            iters1.save(temp);
  //            sprintf(temp,"countsx%i",size);
  //            counts1.save(temp);
  //            sprintf(temp,"planesx%i",size);
  //            planes1.save(temp);
  //            sprintf(temp,"itersy%i",size);
  //            iters2.save(temp);
  //            sprintf(temp,"countsy%i",size);
  //            counts2.save(temp);
  //            sprintf(temp,"planesy%i",size);
  //            planes2.save(temp);
  //            sprintf(temp,"itersz%i",size);
  //            iters3.save(temp);
  //            sprintf(temp,"countsz%i",size);
  //            counts3.save(temp);
  //            sprintf(temp,"planesz%i",size);
  //            planes3.save(temp);
  //        }
  //        sprintf(temp,"itersc%i",size);
  //        itersc.save(temp);
  //        sprintf(temp,"countsc%i",size);
  //        countsc.save(temp);
  //        sprintf(temp,"planesc%i",size);
  //        planesc.save(temp);
  //        sprintf(temp,"iterst%i",size);
  //        iterst.save(temp);
  //        sprintf(temp,"countst%i",size);
  //        countst.save(temp);
  //        sprintf(temp,"planest%i",size);
  //        planest.save(temp);
  //    }

  /* restore keyboard */
  utilities::set_normal_tty();

  saveState(u,p,theta,n,Dt,Ra,Ma,Pr,Lz,rank);
#ifdef HAS_MPI
  MPI_Finalize();
#endif

  return 0; // AY-OH-KAY
}

