#include "mnl/vector.h"
#include "mnl/gll.h"
#include "mnl/matrix.h"
#include "mnl/gordonhall.h"
#include "mnl/function.h"
#include "mnl/cgsolver.h"
#include "mnl/field.h"
#include "mnl/HDF5/hdf5.h"
#include "mnl/timer.h"
#include "mnl/keyboard.h"

using namespace mnl;
using namespace std;

#include "thesis05/functions3d.h"
#include "thesis05/legendrelegendrew.h"
#include "../bigcircle.h"
#include "../consistent.h"
#include "../mass.h"
#include "../sem.h"

void saveState(const utilities::ringBuffer<basics::Field3<basics::matricesStack> >& u,
			   const utilities::ringBuffer<basics::matricesStack>& p, int n, Real Dt, int rank)
{
	stringstream str;
	str << "state-unsteadystokes-pc-" << n << "-" << rank << ".hdf5";
	
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
		source3(basics::function3D& n_X, basics::function3D& n_P, int n_dim) :
			X(n_X), P(n_P), dim(n_dim)
		{
		}

		Real val(Real x, Real y, Real z, Real t) const
		{
			if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
				return 0;

			Real result = -(X.diff2x(x,y,z,t)+X.diff2y(x,y,z,t)+X.diff2z(x,y,z,t))+X.difft(x,y,z,t);
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
		basics::function3D& P;
		int dim;
};

class bugger : public basics::function3D {
	public:
		Real val(Real x, Real y, Real z, Real t) const
		{
			return z;
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
	if (aquired != MPI_THREAD_SERIALIZED ) {
		cout << "threading unavailable (got " << aquired << "!) mpi/openmp combo will not work." << endl;
//        if( omp_get_num_threads() > 1 )
//            return 1;
	} else
		cout << "init threaded mpi successfull " << MPI_THREAD_MULTIPLE << " " << aquired << endl;
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
	Real Lz=1;
	if( argc > 6 )
		Lz = atof(argv[6]);
	bool preconditioned=true;
	if( argc > 7 )
		preconditioned = atoi(argv[7])==1?false:true;
	int d3D=2;
	if( argc > 8 )
		d3D = atoi(argv[8]);

//    legendreLegendreW::poissonSolver SP3(N,N,0,legendreLegendreW::poissonSolver::HomogenousLeft,true,Lz);

	basics::Vector gridGL("GL grid",N-1);
	basics::Vector grid("GLL grid",N+1);
	basics::Vector weight("GLL weight",N+1);
	basics::Vector weightGL("GL weight",N-1);
	utilities::GLL::GaussLobattoLegendreGrid(grid);
	utilities::GLL::GaussLobattoLegendreWeights(weight,grid);
	utilities::GLL::GaussLegendreGrid(gridGL);
	utilities::GLL::GaussLegendreWeights(weightGL,gridGL);
	basics::Matrix GLL2G = utilities::GLL::interpolationMatrix(grid,gridGL);

	bigCircleGeometry circle(N,N,grid,weight);
	bigCircleGeometry3D circle3D(weight,circle);
	circle.setLz(Lz);

	legendreLegendreW::poissonSolver SP(N,N,0,legendreLegendreW::poissonSolver::Homogenous,true,Lz,-1,-1,circle3D.m_level);
	legendreLegendreW::poissonSolver SP2(N,N,0,legendreLegendreW::poissonSolver::Nonhomogenous,true,Lz,-1,-1,circle3D.m_level);

	Real nu = Real(3)/(2*Dt);
	extrudedConsistentPressureEvaluator pressure3DT(circle,SP.Dx(),GLL2G,weightGL,SP.weightX(),
											 	   grid,gridGL,true,rank,size,
											 	   SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET,
												   circle3D.m_level);
	extrudedConsistentPressurePreconditioner pressure3Dpre(pressure3DT,circle3D);
//    extrudedFEMConsistentPressurePreconditioner P3P(circle,weight,weightGL,grid,gridGL,rank,size);
	SEM3DConsistentPressureEvaluator pressure3D(circle3D,SP.Dx(),GLL2G,weightGL,SP.weightX(),
												rank,size,SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
	SEM3DFEMConsistentPressurePreconditioner pressure3Dpre3(circle3D,SP.weightX(),weightGL,grid,gridGL);
	SEMLaplacianEvaluator eval3(circle,SP.Dx(),SP.weightX(),NULL,0,1,
								SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
	extrudedLaplacianEvaluator velocity3DT(eval3,SP,nu,&SP2.Ax(),preconditioned,rank,size,
										   SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
	SEM3DLaplacianEvaluator velocity3D(circle3D,SP.weightX(),SP.Dx(),nu,rank,size,
								  	   SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
	extrudedLaplacianPreconditioner velocity3Dpre(velocity3DT,circle3D);

	basics::coarseGrid desc = circle3D.getDivisionInfo(size);
	utilities::ringBuffer< basics::Field3<basics::matricesStack> > u;
	u.add(*utilities::g_manager.aquireMatricesStackField("velocity 1",N+1,N+1,N+1,desc[rank].elements.size()));
	u.add(*utilities::g_manager.aquireMatricesStackField("velocity 2",N+1,N+1,N+1,desc[rank].elements.size()));
	u.add(*utilities::g_manager.aquireMatricesStackField("velocity 3",N+1,N+1,N+1,desc[rank].elements.size()));

	basics::Field3<basics::matricesStack>& F  = *utilities::g_manager.aquireMatricesStackField("source",N+1,N+1,N+1,desc[rank].elements.size());

	utilities::ringBuffer<basics::matricesStack> p;
	p.add(*utilities::g_manager.aquireMatricesStack("pressure 1",N-1,N-1,N-1,desc[rank].elements.size()));
	p.add(*utilities::g_manager.aquireMatricesStack("pressure 2",N-1,N-1,N-1,desc[rank].elements.size()));
	p.add(*utilities::g_manager.aquireMatricesStack("pressure 3",N-1,N-1,N-1,desc[rank].elements.size()));
	basics::matricesStack& pt = *utilities::g_manager.aquireMatricesStack("pressure temp",N-1,N-1,N-1,desc[rank].elements.size());

	articleTest1X2 articleX;
	articleTest1Y2 articleY;
	articleTest1Z2 articleZ;
    articleTest1P articleP;
//    solX3 articleX;
//    solY3 articleY;
//    solZ3 articleZ;
//    solP3 articleP;
//    sourceX3 sourceX;
//    sourceY3 sourceY;
//    sourceX3 sourceX;
//    sourceY3 sourceY;
//    sourceZ3 sourceZ;
	source3 sourceX(articleX,articleP,0);
	source3 sourceY(articleY,articleP,1);
	source3 sourceZ(articleZ,articleP,2);

//    bugger bug;
//    circle.evaluate(p.get(0),gridGL,bug,0,desc[rank].elements);

	/* set up initial conditions */
	Real t0=0;
//    circle3D.evaluate(u.get(0,-1).X(),articleP,t0-Dt,desc[rank].elements);
//    circle3D.evaluate(u.get(0,-1).Y(),articleP,t0,desc[rank].elements);
//    for( int i=0;i<u.get(0,-1).X().size();++i ) {
//        circle3D.interpolate(p.get(0,-1)[i],u.get(0,-1).X()[i],GLL2G);
//        circle3D.interpolate(p.get(0, 0)[i],u.get(0,-1).Y()[i],GLL2G);
//    }
	
//    circle3D.evaluate(u.get(0,-1),articleX,articleY,articleZ,t0-Dt,desc[rank].elements);
//    circle3D.evaluate(u.get(0, 0),articleX,articleY,articleZ,t0,desc[rank].elements);
	
	cout << "time step:\t\t"  	 << Dt << endl
		 <<	"final time:\t\t" 	 <<  T << endl
		 <<	"Lz:\t\t\t"  		 << Lz << endl
		 <<	"preconditioned:\t\t"<< (preconditioned?"yes":"no")  << endl
		 << "3D evaluator:\t\t"  << (d3D==2?"no":"yes") << endl
		 << "3D preconditioner:\t"  << (d3D==1?"yes":"no") << endl;

	/* setup keyboard */
	utilities::set_raw_tty();

	int elipticu = utilities::g_profiler.add("eliptic velocity");
	int elipticp = utilities::g_profiler.add("eliptic pressure");
	MPIGeometryDotter<basics::geometryStack3D> udotter(circle3D,desc[rank].elements,rank,size);
	MPIDotter pdotter(rank,size);
	int n;
	for( n=0;n<(T-t0)/Dt;++n ) {
		if( rank == 0 )
			cout << "Time: " << t0+(n+1)*Dt << endl;

		/* set up references */
		basics::Field3<basics::matricesStack>& unm1	= u.get(n,-1);
		basics::Field3<basics::matricesStack>& un  	= u.get(n, 0);
		basics::Field3<basics::matricesStack>& unp1	= u.get(n, 1);
		basics::matricesStack& pnm1 				= p.get(n,-1);
		basics::matricesStack& pn  					= p.get(n, 0);
		basics::matricesStack& pnp1 				= p.get(n, 1);

		/* check keyboard */
		if( utilities::keypress() == 's' )
			saveState(u,p,n,Dt,rank);
		if( utilities::keypress() == 'q' )
			break;

//        circle3D.evaluate(F,sourceX,sourceY,sourceZ,t0+(n+1)*Dt,desc[rank].elements);
		F.clear();
		F.Z() = -1;

//        /* BDF2 */
		F.axpy(Real(2)/Dt,un);
		F.axpy(Real(-1)/(2*Dt),unm1);
		circle3D.mass(F,desc[rank].elements);

		/* extrapolate pressure */
		utilities::extrapolate(pnp1,p,n,utilities::FIRST_ORDER);

		pressure3D.m_gradient.evaluate(unp1,pnp1);
		unp1 += F;

		/* delta */
		velocity3D.evaluate(unm1,un,false,false);
		unp1 -= unm1;

		mask(unp1,circle3D,rank,size,velocity3D.m_bc);
		F = unp1;
		dssum(unp1,circle3D,rank,size);
		utilities::g_profiler.start(elipticu);
//        basics::matricesStack* temp = circle3D.localToGlobalZ(unp1.X(),rank,size);
//        basics::matricesStack* temp2 = circle3D.localToGlobalZ(F.X(),rank,size,false,true);
//        velocity3DT.solve(*temp,*temp2);
//        circle3D.globalToLocalZ(unp1.X(),temp,rank,size);
//        utilities::g_manager.unlock(temp2);
//        temp = circle3D.localToGlobalZ(unp1.Y(),rank,size);
//        temp2 = circle3D.localToGlobalZ(F.Y(),rank,size,false,true);
//        velocity3DT.solve(*temp,*temp2);
//        circle3D.globalToLocalZ(unp1.Y(),temp,rank,size);
//        utilities::g_manager.unlock(temp2);
//        temp = circle3D.localToGlobalZ(unp1.Z(),rank,size);
//        temp2 = circle3D.localToGlobalZ(F.Z(),rank,size,false,true);
//        velocity3DT.solve(*temp,*temp2);
//        circle3D.globalToLocalZ(unp1.Z(),temp,rank,size);
//        utilities::g_manager.unlock(temp2);
		velocity3DT.m_name = "velocity X";
		velocity3D.m_nosum = velocity3DT.m_nosum = &F.X();
		cout << "iters velocity X " << utilities::CGSolver::solve(unp1.X(),velocity3D,velocity3Dpre,udotter,1.e-10) << endl;
		velocity3D.m_nosum = velocity3DT.m_nosum = &F.Y();
		velocity3DT.m_name = "velocity Y";
		cout << "iters velocity Y " << utilities::CGSolver::solve(unp1.Y(),velocity3D,velocity3Dpre,udotter,1.e-10) << endl;
		velocity3D.m_nosum = velocity3DT.m_nosum = &F.Z();
		velocity3DT.m_name = "velocity Z";
		cout << "iters velocity Z " << utilities::CGSolver::solve(unp1.Z(),velocity3D,velocity3Dpre,udotter,1.e-10) << endl;
		utilities::g_profiler.pause(elipticu);
		unp1 += un;

		pressure3D.m_divergence.evaluate(pt,unp1);
		pt *= -Real(3)/(2*Dt);
	
//        circle3D.evaluate(u.get(0,-1).Y(),articleY,t0,desc[rank].elements);
//        for( int i=0;i<u.get(0,-1).X().size();++i )
//            circle3D.interpolate(pt[i],u.get(0,-1).X()[i],GLL2G);

//        pt = 1;

		utilities::g_profiler.start(elipticp);
		if( preconditioned ) {
			if( d3D == 1 ) {
				cout << "iterations pressure " << utilities::CGSolver::solve(pt,pressure3D,pressure3Dpre3,pdotter,1.e-6) << endl;
//                cout << "iterations pressure " << utilities::CGSolver::solve(pt,pressure3D,P3P,pdotter,1.e-6) << endl;
			} else if( d3D == 2 ) {
				basics::matricesStack* temp = circle3D.localToGlobalZ(pt,rank,size,true);
				pressure3DT.solve(*temp);
				circle3D.globalToLocalZ(pt,temp,rank,size,true);
			}
			else
				cout << "iterations pressure " << utilities::CGSolver::solve(pt,pressure3D,pressure3Dpre,pdotter,1.e-6) << endl;
		} else {
			if( d3D == 2 ) {
				basics::matricesStack* temp = circle3D.localToGlobalZ(pt,rank,size,true);
				pressure3DT.solve(*temp);
				circle3D.globalToLocalZ(pt,temp,rank,size,true);
			} else
				cout << "iterations pressure " << utilities::CGSolver::solve(pt,pressure3D,pdotter,1.e-6) << endl;
		}

		utilities::g_profiler.pause(elipticp);
		pnp1 += pt;

		/* update velocity */
		pressure3D.m_gradient.evaluate(unm1,pt);
		dssum(unm1,circle3D,rank,size);
		circle3D.invMass(unm1,desc[rank].elements);
		mask(unm1,circle3D,rank,size,velocity3D.m_bc);
		unp1.axpy(2*Dt/3,unm1);
	}

	/* restore keyboard */
	utilities::set_normal_tty();

//    saveState(u,p,n,Dt,rank);

	velocity3DT.m_nu = 1;
	velocity3DT.m_nosum = &F.X();
	SEM3DInverseMassEvaluator invmass3D(circle3D,weightGL,GLL2G,rank,size);
	string err = errorReport(u.get(n),p.get(n),circle3D,articleX,articleY,articleZ,
							 articleP,grid,gridGL,t0+n*Dt,invmass3D,&velocity3D);

	if( rank == 0 ) {
		cout << err;
		cout << utilities::g_profiler.report();
	}
#ifdef HAS_MPI
	MPI_Finalize();
#endif

	return 0; // AY-OH-KAY
}

