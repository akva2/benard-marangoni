#include "mnl/vector.h"
#include "mnl/gll.h"
#include "mnl/matrix.h"
#include "mnl/gordonhall.h"
#include "mnl/function.h"
#include "mnl/cgsolver.h"
#include "mnl/field.h"
#include "mnl/timer.h"

#include "bm/poissonsolver-lgw.h"
#include "bm/sem.h"

#ifdef HAS_MPI
#include <mpi.h>
#endif

#ifdef OPENMP
#include <omp.h>
#endif

using namespace mnl;
using namespace std;
using namespace legendreLegendreW;

#include "bm/functions3d.h"
#include "bm/bigcircle.h"
#include "bm/quadratic.h"

#include "bm/laplacian.h"

#undef RADIUS
#define RADIUS 2

#undef SETUPLOCALS

#define SETUPLOCALS \
		Real r2		= x*x+y*y; \
		Real r		= sqrt(r2); \
		Real sinpr	= sin(M_PI*r); \
		Real cospr	= cos(M_PI*r); \
		Real sinpz	= sin(M_PI*z); \
		Real cospz	= cos(M_PI*z); \
		Real sint	= sin(t); \
		Real cost	= cos(t);

class solTrivial : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return (RADIUS*RADIUS-r2)*sin(M_PI*z);
	}
};

class testMixed : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return 4 + 4*sin(M_PI*z/2)+Real(1)/4*sin(M_PI*z/2)*M_PI*M_PI*(RADIUS*RADIUS-r2);
	}
};

class solMixed : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		return (RADIUS*RADIUS-r2)*(1+sin(M_PI*z));
	}
};

class testMixedN : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		Real sz = sin(M_PI*z/2);

//        if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
			return 0;

//        return 	 4*(1+sz)*cospr*M_PI*M_PI*x*x/(r2*RADIUS*RADIUS)-2*(1+sz)*sinpr*M_PI*x*x/(pow(r2,1.5f)*RADIUS)
//                +4*(1+sz)*sinpr*M_PI/(r*RADIUS)+4*(1+sz)*cospr*M_PI*M_PI*y*y/(r2*RADIUS*RADIUS)
//                -2*(1+sz)*sinpr*M_PI*y*y/(pow(r2,1.5f)*RADIUS)+Real(1)/4*sz*M_PI*M_PI*cospr;

//        if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
//            return -2*RADIUS-2*sz*RADIUS;

//        return -(1+sz)*(RADIUS-x*x/r-r)-(1+sz)*(RADIUS-y*y/r-r)+Real(1)/4*sz*M_PI*M_PI*(RADIUS/2*r2-pow(r,3)/3);
	}
};

class solMixedN : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
	
//        return (1+sin(M_PI*z/2))*cospr;
        return (1+sin(M_PI*z/2))*(RADIUS*r2/2-pow(r,3)/3);
	}
};

class source3 : public basics::function3D {
	public:
		source3(basics::function3D& n_X,Real n_nu) : X(n_X), nu(n_nu) 
		{
		}
		
		Real val(Real x, Real y, Real z, Real t) const
		{
//            return 1;
			if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
				return 0;
//                return 4*M_PI*M_PI-5*M_PI*M_PI*cos(M_PI*z);

//            return 1;
			return( -(X.diff2x(x,y,z,t)+X.diff2y(x,y,z,t)+X.diff2z(x,y,z,t)) + nu*X.val(x,y,z,t) );
		}
	private:
		basics::function3D& X;
		Real nu;
};

class testPeriodicD : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		return sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);	
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		return M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
	}
	
	Real diffy(Real x, Real y, Real z, Real t) const
	{
		return M_PI*sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*z);
	}
	
	Real diffz(Real x, Real y, Real z, Real t) const
	{
		return M_PI*sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*z);
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		return -M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);	
	}

	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		return -M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
	}
	
	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		return -M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);	
	}
};

class testPeriodicV : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		return sin(M_PI*x)*sin(M_PI*y)*(1-cos(M_PI*z));
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		return M_PI*cos(M_PI*x)*sin(M_PI*y)*(1-cos(M_PI*z));
	}
	
	Real diffy(Real x, Real y, Real z, Real t) const
	{
		return M_PI*sin(M_PI*x)*cos(M_PI*y)*(1-cos(M_PI*z));
	}
	
	Real diffz(Real x, Real y, Real z, Real t) const
	{
		return M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		return -M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*(1-cos(M_PI*z));
	}

	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		return -M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*(1-cos(M_PI*z));
	}
	
	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		return M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*z);
	}
};

class articleTest1XV : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		Real z2 = 1-cos(M_PI*z);

		return( cospr*z2 );
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		Real z2 = 1-cos(M_PI*z);

		return z2*(-cospr*pow(M_PI*x,2)/r2+sinpr*M_PI*pow(x,2)/pow(r,3)-sinpr*M_PI/r);
	}

	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		Real z2 = 1-cos(M_PI*z);

		return z2*(-cospr*pow(M_PI*y,2)/r2+sinpr*M_PI*pow(y,2)/pow(r,3)-sinpr*M_PI/r);
	}

	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
	
		cospz = cos(M_PI*z);
		return cospr*cospz*M_PI*M_PI;
	}
};

class articleTest1XV2 : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		Real z2 = 1-cos(M_PI*z);

		return( (1-r*r)*z2 );
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		Real z2 = 1-cos(M_PI*z);

		return -2*z2;
	}

	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		Real z2 = 1-cos(M_PI*z);

		return -2*z2;
	}

	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
	
		cospz = cos(M_PI*z);
		return (1-r*r)*cospz*M_PI*M_PI;
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
	cout << "inited mpi in single threaded mode" << endl;
#endif

	int N;
	if( argc < 2 )
		N = 13;
	else
		N=atoi(argv[1]);
	bool preconditioned=true;
	if( argc > 2 )
		preconditioned=(atoi(argv[2])==1?false:true);
	int d3D=0;
	if( argc > 3 )
		d3D=(atoi(argv[3]));

	basics::Vector grid("GLL grid",N+1);
	utilities::GLL::GaussLobattoLegendreGrid(grid);

	Real Lz = 1;

	poissonSolver::BC bcp = poissonSolver::HomogenousLeft;
	SEMLaplacianEvaluator::BC bc2d = SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN;
	SEMLaplacianEvaluator::BC bc3d = SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM;

	poissonSolver  SP(N,4,0,bcp,true,Lz);
	poissonSolver SP2(N,4,0,poissonSolver::Nonhomogenous,true,Lz);

	bigCircleGeometry circle(N,N,grid,SP.weightX());
//    vector<basics::coarseGrid> group = circle.getCoarseGroups(size);
//    if( rank == 0 ) {
//        for( int n=0;n<size;++n ) {
//            cout << "--------" << n << "------" << endl;
//            for(int k=0;k<group[n].size();++k) {
//                cout << k << " (" << group[n][k].size3 << "): ";
//                for(int l=0;l<group[n][k].elements.size();++l)
//                    cout << group[n][k].elements[l] << " ";
//                cout << endl;
//            }
//        }
//    }
//        exit(1);
//    quadraticGeometry circle(N,N,grid,SP.weightX());
	circle.setLz(Lz);
	Real nu = 100;

//    boxTest3DX article;
//    testPeriodicV article;
	articleTest1XV article;
	source3 source(article,nu);
//    testTrivial mixed;
//    solTrivial msol;

	SEMLaplacianEvaluator eval(circle,SP.Dx(),SP.weightX(),NULL,0,1,bc2d);
	extrudedLaplacianEvaluator eval2(eval,SP,nu,&SP2.Ax(),preconditioned,
									 rank,size,bc3d);
	extrudedFEMLaplacianPreconditioner pre(circle,SP.weightX(),grid,nu,
										   rank,size,bc3d);

	basics::coarseGrid desc = circle.getDivisionInfo(size);
	basics::matricesStack& u     = *utilities::g_manager.aquireMatricesStack("velocity"  ,N+1,N+1,(N+1),desc[rank].elements.size());
	basics::matricesStack& u2    = *utilities::g_manager.aquireMatricesStack("velocity 2",N+1,N+1,(N+1),desc[rank].elements.size());
	basics::matricesStack& exact = *utilities::g_manager.aquireMatricesStack("velocity 2",N+1,N+1,(N+1),desc[rank].elements.size());

	circle.evaluate(u,grid,source,M_PI/2,desc[rank].elements);
//    u = 1;
	circle.mass(u,desc[rank].elements);
	mask(u,circle,rank,size,eval2.m_bc);
	u2 = u;
	dssum(u,circle,rank,size,eval2.m_bc);

	MPIGeometryDotter<basics::geometryStack> udotter(circle,desc[rank].elements,rank,size);

	int whole  = utilities::g_profiler.add("whole");
	int step13 = utilities::g_profiler.add("step13");
	utilities::g_profiler.start(whole);
	if( d3D ) {
		if( preconditioned ) {
			eval2.m_nosum = pre.m_nosum = &u2;
			cout << "iters " << utilities::CGSolver::solve(u,eval2,pre,udotter,1.e-8) << endl;
		} else
			cout << "iters " << utilities::CGSolver::solve(u,eval2,udotter,1.e-8) << endl;
	}
	else
		eval2.solve(u,u2);

	utilities::g_profiler.pause(whole);

	string report = errorReport(u,circle,article,grid,M_PI/2,&eval2);
	if( rank == 0 ) {
		cout << report << endl;
		cout << utilities::g_profiler.report();
	}

#ifdef HAS_MPI
	MPI_Finalize();
#endif

	return 0; // AY-OH-KAY
}

