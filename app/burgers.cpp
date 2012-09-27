#include "mnl/vector.h"
#include "mnl/gll.h"
#include "mnl/matrix.h"
#include "mnl/gordonhall.h"
#include "mnl/function.h"
#include "mnl/cgsolver.h"
#include "mnl/field.h"
#include "mnl/timer.h"

#include "thesis05/poissonsolver-lgw.h"

using namespace mnl;
using namespace std;
using namespace legendreLegendreW;

#include "thesis05/functions3d.h"
#include "../bigcircle.h"
#include "../laplacian.h"
#include "../sem.h"
#include "../convection.h"

#ifdef HAS_MPI
#include <mpi.h>
#endif

class source3 : public basics::function3D {
	public:
		source3(basics::function3D& n_X, basics::function3D& n_Y,
				basics::function3D& n_Z, basics::function3D& n_phi, Real n_nu) : 
			X(n_X), Y(n_Y), Z(n_Z), phi(n_phi), nu(n_nu) 
		{
		}
		
		Real val(Real x, Real y, Real z, Real t) const
		{
			if( x == 0 && y == 0 )
				return 0;

			return( phi.difft(x,y,z,t) - nu*(phi.diff2x(x,y,z,t)+phi.diff2y(x,y,z,t)+phi.diff2z(x,y,z,t)) + X.val(x,y,z,t)*phi.diffx(x,y,z,t) + Y.val(x,y,z,t)*phi.diffy(x,y,z,t) + Z.val(x,y,z,t)*phi.diffz(x,y,z,t) );
		}
	private:
		basics::function3D& X;
		basics::function3D& Y;
		basics::function3D& Z;
		basics::function3D& phi;
		Real nu;
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
	bool preconditioned=true;
	if( argc > 6 )
		preconditioned=(atoi(argv[6])==1?false:true);

	legendreLegendreW::poissonSolver  SP(N,N,0,legendreLegendreW::poissonSolver::Homogenous,true,-1);
	legendreLegendreW::poissonSolver SP2(N,N,0,legendreLegendreW::poissonSolver::Nonhomogenous,true,-1);

	bigCircleGeometry circle(N,N,SP.gridX(),SP.weightX());

	SEMLaplacianEvaluator eval2(circle,SP.Dx(),SP.weightX());
	extrudedLaplacianEvaluator eval(eval2,SP,Real(3)/(Real(2)*Dt),&SP2.Ax(),preconditioned,rank,size);
	convectionEvaluator conv(circle,SP.Dx(),SP.weightX(),rank,size);

    basics::coarseGrid desc = circle.getDivisionInfo(size);
	utilities::ringBuffer<basics::Field3<basics::matricesStack> > u;
	u.add(*utilities::g_manager.aquireMatricesStackField("velocity 1",N+1,N+1,N+1,desc[rank].elements.size()));
	u.add(*utilities::g_manager.aquireMatricesStackField("velocity 2",N+1,N+1,N+1,desc[rank].elements.size()));
	u.add(*utilities::g_manager.aquireMatricesStackField("velocity 3",N+1,N+1,N+1,desc[rank].elements.size()));
	basics::Field3<basics::matricesStack>* F   = utilities::g_manager.aquireMatricesStackField("source",N+1,N+1,N+1,desc[rank].elements.size());
	basics::Field3<basics::matricesStack>* buf = utilities::g_manager.aquireMatricesStackField("source",N+1,N+1,N+1,desc[rank].elements.size());
	
	articleTest1X articleX;
	articleTest1Y articleY;
	articleTest1Z articleZ;
	source3 sourceX(articleX,articleY,articleZ,articleX,1);
	source3 sourceY(articleX,articleY,articleZ,articleY,1);
	source3 sourceZ(articleX,articleY,articleZ,articleZ,1);
	
	/* set up initial conditions */
	Real t0=0;
	circle.evaluate(u.get(0,-1),SP.gridX(),articleX,articleY,articleZ,t0-Dt,desc[rank].elements);
	circle.evaluate(u.get(0, 0),SP.gridX(),articleX,articleY,articleZ,t0,desc[rank].elements);
	circle.maskField(u.get(0,-1),rank,size);
	circle.maskField(u.get(0, 0),rank,size);

	int convection = utilities::g_profiler.add("convection");
	int eliptic = utilities::g_profiler.add("eliptic");

	int n;
	for( n=0;n<(T-t0)/Dt;++n ) {
		cout << "Time: " << t0+(n+1)*Dt << endl;

		/* set up references */
		basics::Field3<basics::matricesStack>& unm1 	= u.get(n,-1);
		basics::Field3<basics::matricesStack>& un  		= u.get(n, 0);
		basics::Field3<basics::matricesStack>& unp1 	= u.get(n, 1);

		/* source */
		circle.evaluate(*F,SP.gridX(),sourceX,sourceY,sourceZ,t0+(n+1)*Dt,desc[rank].elements);

#ifdef HAS_MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif

		/* solve convection problems */
		utilities::g_profiler.start(convection);
		conv.solve(unp1,un,u,n,Dt,Dt,0);
		conv.solve(*buf,unm1,u,n,Dt,Dt,-Dt);
		utilities::g_profiler.pause(convection);
		
		/* BDF2 */
		unp1 *= Real(2)/Dt;
		unp1.axpy(Real(-1)/(Real(2)*Dt),*buf);
		unp1 += *F;

		circle.mass(unp1,desc[rank].elements);

		/* delta value */
		eval.evaluate(unm1,un,false,false);
		unp1 -= unm1;

		circle.maskField(unp1,rank,size);

		*buf = unp1;
		dssum(unp1,circle,size);
		utilities::g_profiler.start(eliptic);
		eval.solve(unp1,*buf);
		utilities::g_profiler.pause(eliptic);
		unp1 += un;
	}

	eval.m_nu = 1;
	string errx = errorReport(u.get(n),circle,articleX,articleY,articleZ,SP.gridX(),n*Dt,&eval);
	if( rank == 0 ) {
		cout << errx;
		cout << utilities::g_profiler.report();
	}

	return 0; // AY-OH-KAY
}

