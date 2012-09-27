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
#include "../quadratic.h"
#include "../laplacian.h"
#include "../sem.h"

class source3 : public basics::function3D {
	public:
		source3(basics::function3D& n_X,Real n_nu) : X(n_X), nu(n_nu) 
		{
		}
		
		Real val(Real x, Real y, Real z, Real t) const
		{
			return( X.difft(x,y,z,t) - nu*(X.diff2x(x,y,z,t)+X.diff2y(x,y,z,t)+X.diff2z(x,y,z,t)) );
		}
	private:
		basics::function3D& X;
		Real nu;
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

	Real Lz = 1;

	basics::Vector grid("GLL grid",N+1);
	utilities::GLL::GaussLobattoLegendreGrid(grid);

	legendreLegendreW::poissonSolver  SP(N,N,0,legendreLegendreW::poissonSolver::Homogenous,true,Lz);
	legendreLegendreW::poissonSolver SP2(N,N,0,legendreLegendreW::poissonSolver::Nonhomogenous,true,Lz);

//    bigCircleGeometry circle(N,N,grid,SP.weightX());
	quadraticGeometry circle(N,N,grid,SP.weightX());
	circle.setLz(Lz);

	SEMLaplacianEvaluator eval2(circle,SP.Dx(),SP.weightX());
	extrudedLaplacianEvaluator eval(eval2,SP,Real(3)/(Real(2)*Dt),&SP2.Ax(),false);

	utilities::ringBuffer<basics::matricesStack> u;
	u.add(*utilities::g_manager.aquireMatricesStack("velocity",circle));
	u.add(*utilities::g_manager.aquireMatricesStack("velocity",circle));
	u.add(*utilities::g_manager.aquireMatricesStack("velocity",circle));

//    articleTest1X article;
	boxTest3DX article;
	source3 source(article,1);
	
	/* set up initial conditions */
	Real t0=0;
	circle.evaluate(u.get(0,-1),grid,article,t0-Dt);
	circle.evaluate(u.get(0, 0),grid,article,t0);
	circle.mask(u.get(0, 0));
	circle.mask(u.get(0,-1));

	int n;
	for( n=0;n<floor((T-t0)/Dt);++n ) {
		std::cout << "Time: " << t0+(n+1)*Dt << std::endl;

		/* set up references */
		basics::matricesStack& unm1 	= u.get(n,-1);
		basics::matricesStack& un  		= u.get(n, 0);
		basics::matricesStack& unp1 	= u.get(n, 1);

		/* source */
		circle.evaluate(unp1,grid,source,t0+(n+1)*Dt);

		/* BDF2 */
		unp1.axpy(Real(2)/Dt,un);
		unp1.axpy(Real(-1)/(Real(2)*Dt),unm1);

		circle.mass(unp1,false);

		/* delta value */
		eval.evaluate(unm1,un,false,false);
		unp1 -= unm1;

		circle.mask(unp1);
		unm1 = unp1;
		circle.dssum(unp1);
		eval.solve(unp1,unm1);
		unp1 += un;
	}

	eval.m_nu = 1;
	std::cout << errorReport(u.get(n),circle,article,grid,T,&eval);

	return 0; // AY-OH-KAY
}

