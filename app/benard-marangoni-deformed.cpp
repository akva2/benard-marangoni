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
			   int& n, Real& Dt, Real& Ra, Real& Ma, Real& Pr, Real& Lz, const string& file,
			   int rank)
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
	reader.read(  p.get(n, 0)  ,"pressure 1");
	reader.read(  p.get(n,-1)  ,"pressure 2");
	reader.read(theta.get(n, 0),"temperature 1");
	reader.read(theta.get(n,-1),"temperature 2");
}

void loadState2(utilities::ringBuffer<basics::Field3<basics::matricesStack> >& u,
			    utilities::ringBuffer<basics::matricesStack>& p,
			    utilities::ringBuffer<basics::matricesStack>& theta,
			    int& n, Real& Dt, Real& Ra, Real& Ma, Real& Pr, Real& Lz, const string& file,
			    int rank, const basics::coarseGrid& elements, int elems)
{
	stringstream str;
//    str << file << "-" << rank << ".hdf5";
	str << file;

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

	n = 0;

	int N = theta.get(0)[0].rows()-1;
	basics::matricesStack* temp
		= utilities::g_manager.aquireMatricesStack("temp",N+1,N+1,N+1,elems);
	basics::matricesStack* temp2
		= utilities::g_manager.aquireMatricesStack("temp",N-1,N-1,N-1,elems);
	reader.read(*temp,"velocity X 1");
	for( int i=0;i<elements[rank].elements.size();++i )
		u.get(0).X()[i] = (*temp)[elements[rank].elements[i]];
	reader.read(*temp,"velocity X 2");
	for( int i=0;i<elements[rank].elements.size();++i )
		u.get(0,-1).X()[i] = (*temp)[elements[rank].elements[i]];
	reader.read(*temp,"velocity Y 1");
	for( int i=0;i<elements[rank].elements.size();++i )
		u.get(0).Y()[i] = (*temp)[elements[rank].elements[i]];
	reader.read(*temp,"velocity Y 2");
	for( int i=0;i<elements[rank].elements.size();++i )
		u.get(0,-1).Y()[i] = (*temp)[elements[rank].elements[i]];
	reader.read(*temp,"velocity Z 1");
	for( int i=0;i<elements[rank].elements.size();++i )
		u.get(0).Z()[i] = (*temp)[elements[rank].elements[i]];
	reader.read(*temp,"velocity Z 2");
	for( int i=0;i<elements[rank].elements.size();++i )
		u.get(0,-1).Z()[i] = (*temp)[elements[rank].elements[i]];
	reader.read(  *temp2,"pressure 1");
	for( int i=0;i<elements[rank].elements.size();++i )
		p.get(0,0)[i] = (*temp2)[elements[rank].elements[i]];
	reader.read(  *temp2,"pressure 2");
	for( int i=0;i<elements[rank].elements.size();++i )
		p.get(0,-1)[i] = (*temp2)[elements[rank].elements[i]];
	reader.read(*temp,"temperature 1");
	for( int i=0;i<elements[rank].elements.size();++i )
		theta.get(0,0)[i] = (*temp)[elements[rank].elements[i]];
	reader.read(*temp,"temperature 2");
	for( int i=0;i<elements[rank].elements.size();++i )
		theta.get(0,-1)[i] = (*temp)[elements[rank].elements[i]];

	utilities::g_manager.unlock(temp);
	utilities::g_manager.unlock(temp2);
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
	if( argc > 2 && (strstr(argv[2],"hdf5") || strstr(argv[2],"state")) )
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

	basics::Vector gridGL("GL grid",N-1);
	basics::Vector grid("GLL grid",N+1);
	basics::Vector weight("GLL weight",N+1);
	basics::Vector weightGL("GL weight",N-1);
	utilities::GLL::GaussLegendreGrid(gridGL);
	utilities::GLL::GaussLobattoLegendreGrid(grid);
	utilities::GLL::GaussLobattoLegendreWeights(weight,grid);
	utilities::GLL::GaussLegendreWeights(weightGL,gridGL);
	basics::Matrix GLL2G = utilities::GLL::interpolationMatrix(grid,gridGL);

	bigCircleGeometry circle(N,N,grid,weight);
	bigCircleGeometry3D circle3D(weight,circle);

	legendreLegendreW::poissonSolver SP(N,N,0,legendreLegendreW::poissonSolver::Homogenous,true,Lz,-1,-1,circle3D.m_level);
	legendreLegendreW::poissonSolver SP2(N,4,0,legendreLegendreW::poissonSolver::Nonhomogenous,true,Lz,-1,-1,circle3D.m_level);
	legendreLegendreW::poissonSolver SP3(N,4,0,legendreLegendreW::poissonSolver::HomogenousLeft,true,Lz,-1,-1,circle3D.m_level);

	basics::coarseGrid desc = circle3D.getDivisionInfo(size);
	basics::coarseGrid desc2D = circle.getDivisionInfo(size);
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
	if( load ) {
		if( strstr(argv[2],"init"))
			loadState2(u,p,theta,n,Dt,Ra,Ma,Pr,Lz,argv[2],rank,desc,circle3D.size());
		else
			loadState(u,p,theta,n,Dt,Ra,Ma,Pr,Lz,argv[2],rank);
	} else {
        /* set up initial conditions */
		Real theta0=1.e-3*100;
		initialTemperature temp(theta0);
//        simplyX temp;
		circle.evaluate(theta.get(0).at(0),grid,temp,t0,desc[rank].elements);
		dssum(theta.get(0),circle,rank,size);

		for( int i=0;i<desc[rank].elements.size();++i )
			basics::multPointwise(theta.get(0).at(0)[i],(*circle.m_mult)[desc[rank].elements[i]]);

		basics::matricesStack* fool = circle3D.localToGlobalZ(theta.get(0),rank,size);
		for( int l=0;l<(*fool)[0].matrices();++l ) {
			for( int i=0;i<desc2D[rank].elements.size();++i ) {
				fool->at(l)[i] = theta.get(0).at(0)[i];
				basics::Matrix& Z0 = circle3D[desc2D[rank].elements[i]].getGH().getMapping().Z()[0];
				basics::Matrix& Z  = circle3D[desc2D[rank].elements[i]].getGH().getMapping().Z()[l];
				for( int j=0;j<(*fool)[0].cols();++j )
					for( int k=0;k<(*fool)[0].rows();++k )
						fool->at(l)[i][j][k] *= (Z[j][k]-Z0[j][k])*(2-Z[j][k]);
			}
		}
		circle3D.globalToLocalZ(theta.get(0),fool,rank,size);
	}

//    basics::matricesStack* stack =
//        utilities::g_manager.aquireMatricesStack("lu",N+1,N+1,N+1,circle3D.size());
//    HDF5::HDF5Reader reader("init.hdf5");
//    reader.read("temperature 1",*stack);
//    for( int i=0;i<desc[rank].elements.size();++i )
//        theta.get(0)[i] = (*stack)[desc[rank].elements[i]];

	bool secondorder=true;

	circle.setLz(Lz);
	/* pressure solver */
	extrudedConsistentPressureEvaluator consistent3DT(circle,SP.Dx(),GLL2G,weightGL,SP.weightX(),
													 grid,gridGL,preconditioned,rank,size,
													 SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM,
													 circle3D.m_level);
    extrudedConsistentPressurePreconditioner consistent3Dpre(consistent3DT,circle3D);
	SEM3DConsistentPressureEvaluator consistent3D(circle3D,SP.Dx(),GLL2G,weightGL,SP.weightX(),
												  rank,size,SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM);

//    SEM3DFEMConsistentPressurePreconditioner consistent3Dpre(circle3D,SP.weightX(),weightGL,grid,gridGL);

	/* velocity solver */
	Real nu;
	if( secondorder )
		nu = Real(3)/(2*Dt*Pr);
	else
		nu = Real(1)/(Dt*Pr);

	SEMLaplacianEvaluator velocity2D(circle,SP.Dx(),SP.weightX(),NULL,0,1,
									 SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
	extrudedLaplacianEvaluator velocity3DT(velocity2D,SP,nu,&SP2.Ax(),preconditioned,
										  rank,size,SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
	extrudedLaplacianEvaluator velocity3DT2(velocity2D,SP,nu,&SP2.Ax(),preconditioned,
										  rank,size,SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
	extrudedLaplacianEvaluator velocity3DNT(velocity2D,SP3,nu,&SP2.Ax(),preconditioned,
										    rank,size,SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM);
	extrudedLaplacianPreconditioner velocity3Dpre(velocity3DT,circle3D);
	extrudedLaplacianPreconditioner velocity3DNpre(velocity3DNT,circle3D);
	SEM3DLaplacianEvaluator velocity3D(circle3D,SP.weightX(),SP.Dx(),nu,rank,size,
									   SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM);
	velocity3DT.m_name = "velocity";

	/* temperature solver */
	SEMLaplacianEvaluator temperature2D(circle,SP.Dx(),SP.weightX(),NULL,0,1,
									  	SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN);
	extrudedLaplacianEvaluator temperature3DT(temperature2D,SP3,nu*Pr,&SP2.Ax(),preconditioned,
										   	 rank,size,SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM);
	extrudedLaplacianPreconditioner temperature3Dpre(temperature3DT,circle3D);
	temperature3DT.m_name = "temperature";

	if( n == 0 ) {
		velocity3DNT.m_nu = velocity3DT.m_nu = Real(1)/(Pr*Dt);
		temperature3DT.m_nu = Real(1)/Dt;
	}

	/* convection solver */
	convectionEvaluator conv(circle,SP.Dx(),SP.weightX(),rank,size,&circle3D);
	if( n == 0 )
		conv.m_order = 0;

	if( rank == 0 )
		cout << "time step:\t\t"  	 << Dt << endl
			 <<	"rayleigh number:\t" << Ra << endl
			 <<	"maragoni number:\t" << Ma << endl
			 <<	"prandtl number:\t\t"<< Pr << endl
			 <<	"final time:\t\t" 	 <<  T << endl
			 <<	"Lz:\t\t\t"  		 << Lz << endl
			 << "level:\t\t\t" 		 << circle3D.m_level << endl
			 <<	"preconditioned:\t\t"<< (preconditioned?"yes":"no") << endl;

	MPIGeometryDotter<basics::geometryStack3D> udotter(circle3D,desc[rank].elements,rank,size);
	MPIDotter pdotter(rank,size);

	/* setup keyboard */
	utilities::set_raw_tty();

	int pressure = utilities::g_profiler.add("pressure solve");
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
		conv.m_bc = velocity3DNT.m_bc;
		conv.solve(unp1.X(),un.X(),u,n,Dt,Dt,0,Dt);
		conv.solve(unp1.Y(),un.Y(),u,n,Dt,Dt,0,Dt);
		conv.m_bc = velocity3DT.m_bc;
		conv.solve(unp1.Z(),un.Z(),u,n,Dt,Dt,0,Dt);
		if( conv.m_order > 0 ) {
			conv.m_bc = velocity3DNT.m_bc;
			conv.solve(buf.X(),unm1.X(),u,n,Dt,Dt,-Dt,Dt);
			conv.solve(buf.Y(),unm1.Y(),u,n,Dt,Dt,-Dt,Dt);
			conv.m_bc = velocity3DT.m_bc;
			conv.solve(buf.Z(),unm1.Z(),u,n,Dt,Dt,-Dt,Dt);
		}

		/* BDF */
		if( conv.m_order > 0 ) {
			unp1 *= Real(2)/(Dt*Pr);
			unp1.axpy(Real(-1)/(2*Dt*Pr),buf);
		}
		else
			unp1 *= Real(1)/(Dt*Pr);

		/* temperature */
		unp1.Z().axpy(Ra,tn);

		circle3D.mass(unp1,desc[rank].elements);
		/* upper tangential boundary conditions */
		vector<basics::matricesStack*> t2D = circle3D.setupFakeStacks(tn,circle3D.m_level,desc2D[rank].elements.size());
		vector<basics::matricesStack*> ux2D = circle3D.setupFakeStacks(unp1.X(),circle3D.m_level,desc2D[rank].elements.size());
		vector<basics::matricesStack*> uy2D = circle3D.setupFakeStacks(unp1.Y(),circle3D.m_level,desc2D[rank].elements.size());
		basics::matrixStack& tcurr = t2D[t2D.size()-1]->at(N);
		basics::matrixStack& uxcurr = ux2D[ux2D.size()-1]->at(N);
		basics::matrixStack& uycurr = uy2D[uy2D.size()-1]->at(N);
		basics::matrixStack& buf1 = buf.X().at(0);
		basics::matrixStack& buf2 = buf.Y().at(0);
		basics::matrixStack& buf3 = buf.X().at(1);
		basics::matrixStack& buf4 = buf.Y().at(1);
		int max = desc2D[rank].elements.size();
#pragma omp parallel for schedule(static)
		for( int i=0;i<max;++i ) {
			int eval = desc2D[rank].elements[i];
			basics::multTranspose(buf1[i],SP.Dx(),tcurr[i],'N','N');
			basics::multTranspose(buf2[i],tcurr[i],SP.Dx(),'N','T');

			/* construct dtheta/dx */
			basics::multPointwise(buf3[i],buf1[i],
								  circle[eval].getGH().getGeometryDerivatives()[3]);
			basics::multPointwise(buf4[i],buf2[i],
								  circle[eval].getGH().getGeometryDerivatives()[2]);
			buf3[i] -= buf4[i];

			/* construct dtheta/dy */
			basics::multPointwise(buf4[i],buf2[i],
								  circle[eval].getGH().getGeometryDerivatives()[0]);
			basics::multPointwise(buf2[i],buf1[i],
								  circle[eval].getGH().getGeometryDerivatives()[1]);
			buf4[i] -= buf2[i];

			massReference(buf3[i],SP.weightX());
			massReference(buf4[i],SP.weightX());

			uxcurr[i].axpy(-Ma,buf3[i]);
			uycurr[i].axpy(-Ma,buf4[i]);
		}
		circle3D.killFakeStacks(t2D);
		circle3D.killFakeStacks(ux2D);
		circle3D.killFakeStacks(uy2D);

		/* delta */
		velocity3D.m_nu = velocity3DNT.m_nu;
		velocity3D.evaluate(buf,un,false,false);
		unp1 -= buf;
		/* extrapolate pressure */
		if( conv.m_order == 0 )
			utilities::extrapolate(pnp1,p,n,utilities::ZEROTH_ORDER);
		else
			utilities::extrapolate(pnp1,p,n,utilities::FIRST_ORDER);

		consistent3D.m_gradient.evaluate(buf,pnp1);
		unp1 += buf;

		/* solve helmholtz problems */
		mask(unp1,circle3D,rank,size,velocity3DNT.m_bc);

		buf = unp1;
		dssum(unp1,circle3D,rank,size);
		if( lumped ) {
			velocity3DNT.m_nosum = velocity3D.m_nosum = &buf.X();
			velocity3D.m_bc = velocity3DNT.m_bc;
			int iter = utilities::CGSolver::solve(unp1.X(),velocity3D,velocity3DNpre,udotter,1.e-8);
			if( rank == 0 )
				cout << "iterations velocity X " << iter << endl;
			velocity3DNT.m_nosum = velocity3D.m_nosum = &buf.Y();
			iter = utilities::CGSolver::solve(unp1.Y(),velocity3D,velocity3DNpre,udotter,1.e-8);
			if( rank == 0 )
				cout << "iterations velocity Y " << iter << endl;
			velocity3DT.m_nosum = velocity3D.m_nosum = &buf.Z();
			velocity3D.m_bc = velocity3DT.m_bc;
			iter = utilities::CGSolver::solve(unp1.Z(),velocity3D,velocity3Dpre,udotter,1.e-8);
			if( rank == 0 )
				cout << "iterations velocity Z " << iter << endl;
		} else {
			velocity3DNT.m_name = "velocity X";
			basics::matricesStack* temp  = circle3D.localToGlobalZ(unp1.X(),rank,size);
			basics::matricesStack* temp2 = circle3D.localToGlobalZ(buf.X(),rank,size,false,true);
			velocity3DNT.solve(*temp,*temp2,0);
			circle3D.globalToLocalZ(unp1.X(),temp,rank,size);
			utilities::g_manager.unlock(temp2);
			velocity3DNT.m_name = "velocity Y";
			temp  = circle3D.localToGlobalZ(unp1.Y(),rank,size);
			temp2 = circle3D.localToGlobalZ(buf.Y(),rank,size,false,true);
			velocity3DNT.solve(*temp,*temp2,1);
			circle3D.globalToLocalZ(unp1.Y(),temp,rank,size);
			utilities::g_manager.unlock(temp2);
			velocity3DT.m_name = "velocity Z";
			temp  = circle3D.localToGlobalZ(unp1.Z(),rank,size);
			temp2 = circle3D.localToGlobalZ(buf.Z(),rank,size,false,true);
			velocity3DT.solve(*temp,*temp2,0);
			circle3D.globalToLocalZ(unp1.Z(),temp,rank,size);
			utilities::g_manager.unlock(temp2);
		}
		unp1 += un;
		mask(unp1,circle3D,rank,size,velocity3DNT.m_bc);

		/* solve pressure problem */
		consistent3D.m_divergence.evaluate(pt,unp1);
		pt *= -velocity3DT.m_nu;
		utilities::g_profiler.start(pressure);
		int iter = utilities::CGSolver::solve(pt,consistent3D,consistent3Dpre,pdotter,1.e-8);
		if( rank == 0 )
			cout << "iterations pressure " << iter << endl;

		/* update velocity */
		consistent3D.m_gradient.evaluate(buf,pt);
		dssum(buf,circle3D,rank,size);
		circle3D.invMass(buf,desc[rank].elements);
		mask(buf,circle3D,rank,size,velocity3DNT.m_bc);
		unp1.axpy(Real(1)/velocity3DT.m_nu,buf);
		pnp1 += pt;

		/* solve temperature convection problem */
		conv.m_bc = temperature3DT.m_bc;
		conv.solve(tnp1,tn,u,n,Dt,Dt,0,Dt);
		if( conv.m_order > 0 )
			conv.solve(buf.X(),tnm1,u,n,Dt,Dt,-Dt,Dt);

		if( conv.m_order > 0 ) {
			tnp1 *= Real(2)/Dt;
			tnp1.axpy(Real(-1)/(2*Dt),buf.X());
		}
		else
			tnp1 *= Real(1)/Dt;
		tnp1 += unp1.Z();
		circle3D.mass(tnp1,desc[rank].elements);

		/* delta */
		velocity3D.m_nu = temperature3DT.m_nu;
		velocity3D.evaluate(buf.X(),tn,false,false);
		tnp1 -= buf.X();
		mask(tnp1,circle3D,rank,size,temperature3DT.m_bc);
		buf.X() = tnp1;
		dssum(tnp1,circle3D,rank,size);

		/* solve the elliptic temperature problem */
		velocity3D.m_bc = temperature3DT.m_bc;
		velocity3D.m_nosum = temperature3DT.m_nosum = &buf.X();
		iter = utilities::CGSolver::solve(tnp1,velocity3D,temperature3Dpre,udotter,1.e-8);
		tnp1 += tn;
		if( rank == 0 )
			cout << "iterations temperature " << iter << endl;

		if( secondorder ) {
			velocity3DT.m_nu = velocity3DNT.m_nu = Real(3)/(2*Dt*Pr);
			temperature3DT.m_nu = Real(3)/(2*Dt);
			conv.m_order = 1;
		}
	}

	if( rank == 0 )
		cout << utilities::g_profiler.report();

	/* restore keyboard */
	utilities::set_normal_tty();

    saveState(u,p,theta,n,Dt,Ra,Ma,Pr,Lz,rank);

#ifdef HAS_MPI
	MPI_Finalize();
#endif

	return 0; // AY-OH-KAY
}

