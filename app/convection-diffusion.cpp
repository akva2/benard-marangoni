#include "mnl/vector.h"
#include "mnl/gll.h"
#include "mnl/matrix.h"
#include "mnl/gordonhall.h"
#include "mnl/function.h"
#include "mnl/cgsolver.h"
#include "mnl/field.h"
#include "mnl/timer.h"

#include "bm/poissonsolver-lgw.h"

using namespace mnl;
using namespace std;
using namespace legendreLegendreW;

#include "bm/functions3d.h"
#include "bm/bigcircle.h"
#include "bm/quadratic.h"
#include "bm/laplacian.h"
#include "bm/sem.h"
#include "bm/convection.h"

#ifdef HAS_MPI
#include <mpi.h>
#endif

class source3 : public basics::function3D {
	public:
		source3(basics::function3D& n_X, basics::function3D& n_Y,
				basics::function3D& n_Z, basics::function3D& n_phi) : 
			X(n_X), Y(n_Y), Z(n_Z), phi(n_phi) 
		{
		}

		Real val(Real x, Real y, Real z, Real t) const
		{
//            if( x == 0 && y == 0 )
//                return 0;

			return(  phi.difft(x,y,z,t) - (phi.diff2x(x,y,z,t)+phi.diff2y(x,y,z,t)+phi.diff2z(x,y,z,t)) 
					+X.val(x,y,z,t)*phi.diffx(x,y,z,t)
					+Y.val(x,y,z,t)*phi.diffy(x,y,z,t)
					+Z.val(x,y,z,t)*phi.diffz(x,y,z,t) );
		}
	private:
		basics::function3D& X;
		basics::function3D& Y;
		basics::function3D& Z;
		basics::function3D& phi;
};

//class articleTest1XV : public basics::function3D {
//public:
//    Real val(Real x, Real y, Real z, Real t) const
//    {
//        SETUPLOCALS

//        cospz = 1-cos(2*M_PI*z);

//        if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
//            return 0;

//        return( (Real(2)/5)*pow(cospr,2)*y*cospz*sint*x/r2 );
//    }

//    Real diffx(Real x, Real y, Real z, Real t) const
//    {
//        SETUPLOCALS
//        if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
//            return 0;

//        cospz = 1-cos(2*M_PI*z);
//        return( -(Real(4)/5)*sinpr*y*cospz*sint*pow(x,2)*cospr*M_PI/pow(r2,(Real(3)/2))
//                +(Real(2)/5)*pow(cospr,2)*y*cospz*sint/r2
//                -(Real(4)/5)*pow(cospr,2)*y*cospz*sint*pow(x,2)/pow(r2,2) );
//    }

//    Real diffy(Real x, Real y, Real z, Real t) const
//    {
//        SETUPLOCALS
//        if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
//            return 0;

//        cospz = 1-cos(2*M_PI*z);
//        return -(Real(4)/5)*sinpr*pow(y,2)*cospz*sint*x*cospr*M_PI/pow(r2,Real(3)/2)
//               +(Real(2)/5)*pow(cospr,2)*cospz*sint*x/r2
//               -(Real(4)/5)*pow(cospr,2)*pow(y,2)*cospz*sint*x/pow(r2,2);
//    }

//    Real diffz(Real x, Real y, Real z, Real t) const
//    {
//        SETUPLOCALS
//        if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
//            return 0;

//        sinpz = sin(2*M_PI*z);
//        return (Real(4)/5)*pow(cospr,2)*y*sinpz*M_PI*sint*x/r2;
//    }

//    Real difft(Real x, Real y, Real z, Real t) const
//    {
//        if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
//            return 0;

//        SETUPLOCALS

//        cospz = 1-cos(2*M_PI*z);
//        return( (Real(2)/5)*pow(cospr,2)*y*cospz*cost*x/r2 );
//    }

//    Real diff2x(Real x, Real y, Real z, Real t) const
//    {
//        SETUPLOCALS
//        if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
//            return 0;

//        cospz = 1-cos(2*M_PI*z);
//        return( -(Real( 4)/5)*pow(cospr,2)*pow(M_PI,2)*pow(x,3)*y*cospz*sint/pow(r2,2)
//                -(Real(12)/5)*sinpr*y*cospz*sint*x*cospr*M_PI/pow(r2,Real(3)/2)
//                +(Real( 4)/1)*sinpr*y*cospz*sint*pow(x,3)*cospr*M_PI/pow(r2,Real(5)/2)
//                +(Real( 4)/5)*pow(sinpr,2)*y*cospz*sint*pow(x,3)*pow(M_PI,2)/pow(r2,2)
//                -(Real(12)/5)*pow(cospr,2)*y*cospz*sint*x/pow(r2,2)
//                +(Real(16)/5)*pow(cospr,2)*y*cospz*sint*pow(x,3)/pow(r2,3) );
//    }

//    Real diff2y(Real x, Real y, Real z, Real t) const
//    {
//        SETUPLOCALS
//        if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
//            return 0;

//        cospz = 1-cos(2*M_PI*z);
//        return( -(Real( 4)/5)*pow(cospr,2)*pow(M_PI,2)*pow(y,3)*cospz*sint*x/pow(r2,2)
//                -(Real(12)/5)*sinpr*y*cospz*sint*x*cospr*M_PI/pow(r2,Real(3)/2)
//                +(Real( 4)/1)*sinpr*pow(y,3)*cospz*sint*x*cospr*M_PI/pow(r2,Real(5)/2)
//                +(Real( 4)/5)*pow(sinpr,2)*pow(y,3)*cospz*sint*x*pow(M_PI,2)/pow(r2,2)
//                -(Real(12)/5)*pow(cospr,2)*y*cospz*sint*x/pow(r2,2)
//                +(Real(16)/5)*pow(cospr,2)*pow(y,3)*cospz*sint*x/pow(r2,3) );
//    }

//    Real diff2z(Real x, Real y, Real z, Real t) const
//    {
//        SETUPLOCALS
//        if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
//            return 0;

//        cospz = cos(2*M_PI*z);
//        return( (Real(8)/5)*pow(cospr,2)*y*cospz*pow(M_PI,2)*sint*x/r2 );
//    }
//};

//class articleTest1XV : public basics::function3D {
//public:
//    Real val(Real x, Real y, Real z, Real t) const
//    {
//        return sin(t)*(z*z/2-pow(z,3)/3);
//    }

//    Real difft(Real x, Real y, Real z, Real t) const
//    {
//        return cos(t)*(z*z/2-pow(z,3)/3);
//    }

//    Real diffz(Real x, Real y, Real z, Real t) const
//    {
//        return sin(t)*(z-z*z);
//    }

//    Real diff2z(Real x, Real y, Real z, Real t) const
//    {
//        return sin(t)*(1-2*z);
//    }
//};

//class articleTest1XV : public basics::function3D {
//public:
//    Real val(Real x, Real y, Real z, Real t) const
//    {
//        return sin(t)*(4-x*x-y*y)*(z*z/2-pow(z,3)/3);
//    }

//    Real difft(Real x, Real y, Real z, Real t) const
//    {
//        return cos(t)*(4-x*x-y*y)*(z*z/2-pow(z,3)/3);
//    }

//    Real diffx(Real x, Real y, Real z, Real t) const
//    {
//        return -2*x*(z*z/2-pow(z,3)/3)*sin(t);
//    }
//    
//    Real diffy(Real x, Real y, Real z, Real t) const
//    {
//        return -2*y*(z*z/2-pow(z,3)/3)*sin(t);
//    }

//    Real diffz(Real x, Real y, Real z, Real t) const
//    {
//        return sin(t)*(4-x*x-y*y)*(z-z*z);
//    }
//    
//    Real diff2x(Real x, Real y, Real z, Real t) const
//    {
//        return (-z*z+Real(2)/3*pow(z,3))*sin(t);
//    }
//    
//    Real diff2y(Real x, Real y, Real z, Real t) const
//    {
//        return (-z*z+Real(2)/3*pow(z,3))*sin(t);
//    }

//    Real diff2z(Real x, Real y, Real z, Real t) const
//    {
//        return sin(t)*(4-x*x-y*y)*(1-2*z);
//    }
//};

class articleTest1XV : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		return sin(t)*(4-x*x-y*y)*z*(1-z);
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		return cos(t)*(4-x*x-y*y)*z*(1-z);
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		return -2*x*z*(1-z)*sin(t);
	}
	
	Real diffy(Real x, Real y, Real z, Real t) const
	{
		return -2*y*z*(1-z)*sin(t);
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		return sin(t)*(4-x*x-y*y)*(1-2*z);
	}
	
	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		return -2*z*(1-z)*sin(t);
	}
	
	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		return -2*z*(1-z)*sin(t);
	}

	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		return (-8+2*x*x+2*y*y)*sin(t);
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

	Real Lz = 1;

	legendreLegendreW::poissonSolver  SP(N,N,0,legendreLegendreW::poissonSolver::Homogenous,true,Lz);
	legendreLegendreW::poissonSolver SP2(N,N,0,legendreLegendreW::poissonSolver::Nonhomogenous,true,Lz);
	
	basics::Vector grid("GLL grid",N+1);
	utilities::GLL::GaussLobattoLegendreGrid(grid);

//    bigCircleGeometry circle(N,N,grid,SP.weightX());
	quadraticGeometry circle(N,N,grid,SP.weightX());
	circle.setLz(Lz);

//    Real nu = Real(1)/(Dt);
	Real nu = Real(3)/(2*Dt);
	SEMLaplacianEvaluator eval2(circle,SP.Dx(),SP.weightX(),NULL,0,1,
								SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
	extrudedLaplacianEvaluator eval(eval2,SP,nu,&SP2.Ax(),preconditioned,rank,size,
									SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
	convectionEvaluator conv(circle,SP.Dx(),SP.weightX(),rank,size);

    basics::coarseGrid desc = circle.getDivisionInfo(size);
	utilities::ringBuffer<basics::Field3<basics::matricesStack> > u;
	u.add(*utilities::g_manager.aquireMatricesStackField("velocity 1",N+1,N+1,N+1,desc[rank].elements.size()));
	u.add(*utilities::g_manager.aquireMatricesStackField("velocity 2",N+1,N+1,N+1,desc[rank].elements.size()));
	utilities::ringBuffer<basics::matricesStack> phi;
	phi.add(*utilities::g_manager.aquireMatricesStack("phi 1",N+1,N+1,N+1,desc[rank].elements.size()));
	phi.add(*utilities::g_manager.aquireMatricesStack("phi 2",N+1,N+1,N+1,desc[rank].elements.size()));
	phi.add(*utilities::g_manager.aquireMatricesStack("phi 3",N+1,N+1,N+1,desc[rank].elements.size()));
	basics::matricesStack* F   = utilities::g_manager.aquireMatricesStack("source",N+1,N+1,N+1,desc[rank].elements.size());

//    articleTest1XV articleX;
//    articleTest1XV articleY;
//    articleTest1XV article;
//    articleTest1XV phifunc;
	boxTest3DX articleX;
	boxTest3DY articleY;
	boxTest3DZ articleZ;
	boxTest3DX phifunc;
	source3 sourceX(articleX,articleY,articleZ,phifunc);

	/* set up initial conditions */
	Real t0=0;
	circle.evaluate(  u.get(0,-1),grid,articleX,articleY,articleZ,t0-Dt,desc[rank].elements);
	circle.evaluate(  u.get(0, 0),grid,articleX,articleY,articleZ,t0,desc[rank].elements);
	circle.evaluate(phi.get(0,-1),grid,phifunc,t0-Dt,desc[rank].elements);
	circle.evaluate(phi.get(0, 0),grid,phifunc,t0,desc[rank].elements);

	int convection = utilities::g_profiler.add("convection");
	int eliptic = utilities::g_profiler.add("eliptic");

	int n;
	for( n=0;n<(T-t0)/Dt;++n ) {
		cout << "Time: " << t0+(n+1)*Dt << endl;

		/* set up references */
		basics::Field3<basics::matricesStack>& unm1 = u.get(n,-1);
		basics::Field3<basics::matricesStack>& un  	= u.get(n, 0);
		basics::matricesStack& pnm1 				= phi.get(n,-1);
		basics::matricesStack& pn  					= phi.get(n, 0);
		basics::matricesStack& pnp1 				= phi.get(n, 1);

		/* source */
		circle.evaluate(  un,grid,articleX,articleY,articleZ,t0+n*Dt,desc[rank].elements);

		/* solve convection problems */
		utilities::g_profiler.start(convection);
		conv.m_bc = eval.m_bc;
		conv.solve(pnp1,pn,u,n,Dt,Dt,0);
		conv.solve(*F,pnm1,u,n,Dt,Dt,-Dt);
		utilities::g_profiler.pause(convection);

//        /* BDF1 */
//        pnp1 *= Real(1)/Dt;

		/* BDF2 */
		pnp1 *= Real(2)/Dt;
		pnp1.axpy(Real(-1)/(Real(2)*Dt),*F);

		circle.evaluate(*F,grid,sourceX,t0+(n+1)*Dt,desc[rank].elements);
		pnp1 += *F;

		circle.mass(pnp1,desc[rank].elements);

		/* delta value */
		eval.evaluate(pnm1,pn,false,false);
		pnp1 -= pnm1;

		mask(pnp1,circle,rank,size,eval.m_bc);
		*F = pnp1;
		dssum(pnp1,circle,rank,size);

		utilities::g_profiler.start(eliptic);
		eval.solve(pnp1,*F);
		utilities::g_profiler.pause(eliptic);
		pnp1 += pn;
	}

	eval.m_nu = 1;
	string errx = errorReport(phi.get(n),circle,phifunc,grid,t0+n*Dt,&eval);
	if( rank == 0 ) {
		cout << errx;
		cout << utilities::g_profiler.report();
	}

	return 0; // AY-OH-KAY
}

