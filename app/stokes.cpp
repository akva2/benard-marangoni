#include "mnl/vector.h"
#include "mnl/gll.h"
#include "mnl/matrix.h"
#include "mnl/gordonhall.h"
#include "mnl/function.h"
#include "mnl/cgsolver.h"
#include "mnl/field.h"
#include "mnl/hdf5.h"
#include "mnl/timer.h"

using namespace mnl;
using namespace std;

#include "bm/bigcircle.h"
#include "bm/functions3d.h"
#include "bm/legendrelegendrew.h"
#include "bm/uzawa.h"
#include "bm/mass.h"
#include "bm/sem.h"

#ifdef HAS_MPI
#include <mpi.h>
#endif

using namespace legendreLegendreW;
		
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

			Real result = -(X.diff2x(x,y,z,t)+X.diff2y(x,y,z,t)+X.diff2z(x,y,z,t));
			if( dim == 0 )
				result += P.diffx(x,y,z,t);
			if( dim == 1 )
				result += P.diffy(x,y,z,t);
			if( dim == 2 )
				result += P.diffz(x,y,z,t);

			return( result );
		}

		Real diffx(Real x, Real y, Real z, Real t) const
		{
			return( 0 );
		}

		Real diffy(Real x, Real y, Real z, Real t) const
		{
			return( 0 );
		}
		
		Real diffz(Real x, Real y, Real z, Real t) const
		{
			return( 0 );
		}

		Real difft(Real x, Real y, Real z, Real t) const
		{
			return( 0.f ); 
		}

		Real diff2x(Real x, Real y, Real z, Real t) const
		{
			return( 0 );
		}

		Real diff2y(Real x, Real y, Real z, Real t) const
		{
			return( 0 );
		}
		
		Real diff2z(Real x, Real y, Real z, Real t) const
		{
			return( 0 );
		}
	private:
		basics::function3D& X;
		basics::function3D& P;
		int dim;
};


int main(int argc, char** argv)
{
	int rank=0, size=1;

#ifdef HAS_MPI
	int aquired;
	MPI_Init_thread(&argc, &argv,MPI_THREAD_MULTIPLE,&aquired);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	if (aquired != MPI_THREAD_MULTIPLE ) {
		cout << "threading unavailable (got " << aquired << "!) mpi/openmp combo will not work." << endl;
		if( omp_get_num_threads() > 1 && size > 1)
			return 1;
	} else
		cout << "init threaded mpi successfull " << MPI_THREAD_MULTIPLE << " " << aquired << endl;
#endif
	int N;
	if( argc < 2)
		N=12;
	else
		N = atoi(argv[1]);
	bool preconditioned=true;
	if( argc > 2 )
		preconditioned = (atoi(argv[2])==1?false:true);

	basics::Vector gridGL("GL grid",N-1);
	basics::Vector weightGL("GL weight",N-1);
	utilities::GLL::GaussLegendreGrid(gridGL);
	utilities::GLL::GaussLegendreWeights(weightGL,gridGL);
	poissonSolver  SP(N,4,0,poissonSolver::Homogenous,true,-1);
	poissonSolver SP2(N,4,0,poissonSolver::Nonhomogenous,true,-1);
	basics::Matrix GLL2G = utilities::GLL::interpolationMatrix(SP.gridX(),gridGL);

	bigCircleGeometry circle(N,N,SP.gridX(),SP.weightX());

	extrudedUzawaEvaluator    	  eval (circle,GLL2G,weightGL,SP,&SP2.Ax(),preconditioned,0,rank,size);
	SEMInverseMassEvaluator 	  eval3(circle,weightGL,GLL2G,rank,size);
	extrudedInverseMassEvaluator  eval2(eval3);
	vector<basics::geometryStack::coarseDescriptor> desc = circle.getDivisionInfo(size);

	int total = utilities::g_profiler.add("total time");

	basics::Field3<basics::matricesStack>& u   = *utilities::g_manager.aquireMatricesStackField("velocity 1",N+1,N+1,N+1,desc[rank].elements.size());
	basics::matricesStack& u2  = *utilities::g_manager.aquireMatricesStack("velocity 2",N+1,N+1,N+1,desc[rank].elements.size());
	
	basics::Field3<basics::matricesStack>& F = *utilities::g_manager.aquireMatricesStackField("source x",N+1,N+1,N+1,desc[rank].elements.size());

	basics::matricesStack& p  = *utilities::g_manager.aquireMatricesStack("pressure 1",N-1,N-1,N-1,desc[rank].elements.size());

	articleTest1X articleX;
	articleTest1Y articleY;
	articleTest1Z articleZ;
	articleTest1P articleP;
	source3 sourceX(articleX,articleP,0);
	source3 sourceY(articleY,articleP,1);
	source3 sourceZ(articleZ,articleP,2);

	circle.evaluate(F,SP.gridX(),sourceX,sourceY,sourceZ,M_PI/2,desc[rank].elements);
	MPIDotter pdotter(rank,size);

	/* B */
	circle.mass(F,desc[rank].elements);
	circle.maskField(F,rank,size);
	u = F;
	/* A^-1 */
	u2 = F.X();
	circle.dssum(u.X(),rank,size);
	eval.m_laplacian.solve(u.X(),u2);
	u2 = F.Y();
	circle.dssum(u.Y(),rank,size);
	eval.m_laplacian.solve(u.Y(),u2);
	u2 = F.Z();
	circle.dssum(u.Z(),rank,size);
	eval.m_laplacian.solve(u.Z(),u2);
	/* -D */
	eval.m_divergence.evaluate(p,u);
	p *= -1;

	utilities::g_profiler.start(total);
	std::cout << "outer iterations " << utilities::CGSolver::solve(p,eval,eval2,pdotter,1.e-8) << std::endl;

	eval.m_gradient.evaluate(u,p,false);
	u += F;
	circle.maskField(u,rank,size);
	u2 = u.X();
	circle.dssum(u.X(),rank,size);
	eval.m_laplacian.solve(u.X(),u2);
	u2 = u.Y();
	circle.dssum(u.Y(),rank,size);
	eval.m_laplacian.solve(u.Y(),u2);
	u2 = u.Z();
	circle.dssum(u.Z(),rank,size);
	eval.m_laplacian.solve(u.Z(),u2);
	utilities::g_profiler.pause(total);

	eval.m_laplacian.m_nu = 1;
	string report = errorReport(u,p,circle,articleX,articleY,articleZ,articleP,
								SP.gridX(),gridGL,M_PI/2,eval2,&eval.m_laplacian);
	if( rank == 0 ) {
		cout << report;
		cout << utilities::g_profiler.report();
	}

#ifdef HAS_MPI
	MPI_Finalize();
#endif

	return 0; // AY-OH-KAY
}

