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

#include <signal.h>

using namespace mnl;
using namespace std;

#include "bm/functions3d.h"
#include "bm/legendrelegendrew.h"
#include "bm/bigcircle.h"
#include "bm/quadratic.h"
#include "bm/consistent.h"
#include "bm/mass.h"
#include "bm/sem.h"
	
int n;
utilities::ringBuffer< basics::Field3<basics::matricesStack> > u;
utilities::ringBuffer<basics::matricesStack> p;
Real Dt = 1.e-2f;

void saveState(const utilities::ringBuffer<basics::Field3<basics::matricesStack> >& u,
			   const utilities::ringBuffer<basics::matricesStack>& p, int n, Real Dt, int rank)
{
	stringstream str;
	str << "state-unsteadystokes-pc-" << n << "-" << rank << ".hdf5";
	
	cout << "saving state to " << str.str() << endl;

	HDF5::HDF5Writer writer(str.str());
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

void loadState(utilities::ringBuffer<basics::Field3<basics::matricesStack> >& u,
			   utilities::ringBuffer<basics::matricesStack>& p, int& n, Real& Dt, int rank,
			   const char* file)
{
	cout << "loading state from " << file << endl;

	HDF5::HDF5Reader reader(file);
	reader.read(u.get(n, 0).X(),"velocity X 1");
	reader.read(u.get(n,-1).X(),"velocity X 2");
	reader.read(u.get(n, 0).Y(),"velocity Y 1");
	reader.read(u.get(n, 0).Y(),"velocity Y 2");
	reader.read(u.get(n, 0).Z(),"velocity Z 1");
	reader.read(u.get(n,-1).Z(),"velocity Z 2");
	reader.read(p.get(n, 0)    ,"pressure 1");
	reader.read(p.get(n,-1)    ,"pressure 2");
	basics::Vector time("time info",2);
	reader.read(time,"time info");
	n = time[0];
	Dt = time[1];
}

class source3 : public basics::function3D {
	public:
		source3(basics::function3D& n_X, basics::function3D& n_P, int n_dim) :
			X(n_X), P(n_P), dim(n_dim)
		{
		}

		Real val(Real x, Real y, Real z, Real t) const
		{
			if( abs(x) < 1.e-10 && abs(y) < 1.e-10 )
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

/*
class sourceX3 : public basics::function3D {
	public:

		Real val(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS

			return 	-Real( 4)/5*cospr*cospr*pow(M_PI,2)*pow(x,3)*y*sinpz*sint/pow(r2,2)
					+Real( 4)/1*sinpr*y*sinpz*sint*pow(x,3)*cospr*M_PI/pow(r2,Real(5)/2)
					-Real(24)/5*sinpr*y*sinpz*sint*cospr*M_PI/pow(r2,Real(3)/2)
					+Real( 4)/5*sinpr*sinpr*y*sinpz*sint*pow(x,3)*M_PI*M_PI/pow(r2,2)
					-Real(16)/5*sinpr*sinpr*sinpz*sint*pow(x,3)/pow(r2,3)
					+Real(24)/5*sinpr*sinpr*y*sinpz*sint*x/pow(r2,2)
					-Real( 4)/5*cospr*cospr*M_PI*M_PI*pow(y,3)*sinpz*sint*x/pow(r2,2)
					+Real( 4)/1*sinpr*pow(y,3)*sinpz*sint*x*cospr*M_PI/pow(r2,Real(5)/2)
					+Real( 4)/5*sinpr*sinpr*pow(y,3)*sint*x*M_PI*M_PI/pow(r2,2)
					-Real(16)/5*sinpr*sinpr*pow(y,3)*sinpz*sint*x/pow(r2,3)
					+Real( 8)/5*sinpr*sinpr*y*sinpz*M_PI*M_PI*sint*x/r2
					+Real( 2)/1*sinpr*sinpz*sint*cospr*M_PI*x/r
					+Real( 2)/5*sinpr*sinpr*y*sinpz*cost*x/r2;
		}
};

class sourceY3 : public basics::function3D {
	public:
		Real val(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS

			return 	 Real(12)/5*sinpr*sinpr*y*y*sinpz*sint/pow(r2,2)
					-Real( 1)/5*sinpr*sinpr*x*x*sinpz*cost/r2
					+Real( 2)/1*sinpr*pow(y,4)*sinpz*sint*cospr*M_PI/pow(r2,Real(5)/2)
					+Real( 2)/5*sinpr*sinpr*pow(y,4)*sinpz*sint*M_PI*M_PI/pow(r2,2)
					+Real(12)/5*sinpr*x*x*sinpz*sint*cospr*M_PI/pow(r2,Real(3)/2)
					+Real( 2)/5*cospr*cospr*M_PI*M_PI*pow(x,4)*sinpz*sint/pow(r2,2)
					-Real( 2)/1*sinpr*pow(x,4)*sinpz*sint*cospr*M_PI/pow(r2,Real(5)/2)
					-Real( 8)/5*sinpr*sinpr*pow(y,4)*sinpz*sint/pow(r2,3)
					-Real( 2)/5*sinpr*sinpr*pow(x,4)*sinpz*sint*M_PI*M_PI/pow(r2,2)
					-Real(12)/5*sinpr*sinpr*x*x*sinpz*sint/pow(r2,2)
					+Real( 4)/5*sinpr*sinpr*y*y*sinpz*M_PI*M_PI*sint/r2
					+Real( 8)/5*sinpr*sinpr*pow(x,4)*sinpz*sint/pow(r2,3)
					-Real( 4)/5*sinpr*sinpr*x*x*sinpz*M_PI*M_PI*sint/r2
					+Real( 2)/1*sinpr*sinpz*sint*cospr*M_PI*y/r
					+Real( 1)/5*sinpr*sinpr*y*y*sinpz*cost/r2
					-Real( 2)/5*cospr*cospr*M_PI*M_PI*pow(y,4)*sinpz*sint/pow(r2,2)
					-Real(12)/5*sinpr*y*y*sinpz*sint*cospr*M_PI/pow(r2,Real(3)/2);
		}
};
*/

class sourceX3 : public basics::function3D {
	public:
		Real val(Real x, Real y, Real z, Real t) const
		{
			if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
				return 0;

			SETUPLOCALS

			cospz = cos(M_PI*(z+1)/2);
			
			return 	 Real(12)/5*sinpr*y*cospz*sint*x*cospr*M_PI/pow(r,3)
					+Real( 2)/5*pow(cospr,2)*pow(M_PI,2)*pow(x,3)*y*cospz*sint/pow(r,4)
					-Real(1)/20*pow(sinpr,2)*y*cospz*pow(M_PI,2)*sint*x/r2
					-Real(1)/20*pow(sinpr,2)*x*cospz*pow(M_PI,2)*sint*y/r2
					-Real( 1)/5*pow(sinpr,2)*y*cospz*cost*x/r2
					-Real( 1)/5*pow(sinpr,2)*x*cospz*cost*y/r2
					-Real( 2)/1*sinpr*y*cospz*sint*pow(x,3)*cospr*M_PI/pow(r,5)
					-Real( 2)/5*pow(sinpr,2)*y*cospz*sint*pow(x,3)*pow(M_PI,2)/pow(r,4)
					-Real(12)/5*pow(sinpr,2)*y*cospz*sint*x/pow(r,4)
					+Real( 8)/5*pow(sinpr,2)*y*cospz*sint*pow(x,3)/pow(r,6)
					+Real( 2)/5*pow(cospr,2)*pow(M_PI,2)*pow(y,3)*x*cospz*sint/pow(r,4)
					-Real( 2)/1*sinpr*x*cospz*sint*pow(y,3)*cospr*M_PI/pow(r,5)
					-Real( 2)/5*pow(sinpr,2)*x*cospz*sint*pow(y,3)*pow(M_PI,2)/pow(r,4)
					+Real( 8)/5*pow(sinpr,2)*x*cospz*sint*pow(y,3)/pow(r,6)
					+Real( 2)/5*pow(cospr,2)*pow(M_PI,2)*pow(x,3)*cospz*sint*y/pow(r,4)
					-Real( 2)/1*sinpr*pow(x,3)*cospz*sint*y*cospr*M_PI/pow(r,5)
					-Real( 2)/5*pow(sinpr,2)*x*cospz*sint*pow(y,3)*pow(M_PI,2)/pow(r,4)
					-Real(12)/5*pow(sinpr,2)*cospz*sint*y*x/pow(r,4)
					+Real( 8)/5*pow(sinpr,2)*pow(x,3)*cospz*sint*y/pow(r,6)
					+Real( 2)/5*pow(cospz,2)*pow(M_PI,2)*pow(y,3)*cospz*sint*x/pow(r,4)
					-Real( 2)/1*sinpr*pow(y,3)*cospz*sint*x*cospr*M_PI/pow(r,5)
					-Real( 2)/5*pow(sinpr,2)*pow(y,3)*cospz*sint*x*pow(M_PI,2)/pow(r,4);
		}
};

class sourceY3 : public basics::function3D {
	public:
		Real val(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS
			cospz = cos(M_PI*(z+1)/2);

			if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
				return 	-Real(2)/5*sin(M_PI*z/2)*pow(M_PI,2)*sint
						+Real(2)/5*sint*pow(M_PI,2)*sin(M_PI*z/2);
			

			return 	-Real( 2)/1*sinpr*pow(y,4)*cospz*sint*cospr*M_PI/pow(r,5)
					-Real( 2)/5*pow(sinpr,2)*pow(y,4)*cospz*sint*pow(M_PI,2)/pow(r,4)
					+Real( 2)/5*pow(cospr,2)*pow(M_PI,2)*pow(y,4)*cospz*sint/pow(r,4)
					+Real( 2)/5*pow(sinpr,2)*cospz*sint/r2
					+Real( 8)/5*pow(sinpr,2)*pow(y,4)*cospz*sint/pow(r,6)
					+Real( 2)/5*pow(sinpr,2)*pow(x,2)*cospz*sint*pow(M_PI,2)*pow(y,2)/pow(r,4)
					-Real( 2)/5*pow(cospr,2)*pow(M_PI,2)*pow(y,2)*pow(x,2)*cospz*sint/pow(r,4)
					+Real( 2)/1*sinpr*pow(x,2)*cospz*sint*cospr*M_PI*pow(y,2)/pow(r,5)
					-Real(12)/5*sinpr*pow(x,2)*cospz*sint*cospr*M_PI/pow(r,3)
					-Real( 8)/5*pow(sinpr,2)*pow(x,2)*cospz*sint*pow(y,2)/pow(r,6)
					-Real(1)/20*pow(sinpr,2)*pow(y,2)*cospz*pow(M_PI,2)*sint/r2
					+Real(1)/20*pow(sinpr,2)*pow(x,2)*cospz*pow(M_PI,2)*sint/r2
					-Real( 1)/5*pow(sinpr,2)*pow(y,2)*cospz*cost/r2
					+Real( 1)/5*pow(sinpr,2)*pow(x,2)*cospz*cost/r2
					-Real( 2)/5*pow(cospr,2)*pow(M_PI,2)*pow(x,4)*cospz*sint/pow(r,4)
					+Real( 2)/1*sinpr*pow(x,4)*cospz*sint*cospr*M_PI/pow(r,5)
					+Real( 2)/5*pow(sinpr,2)*pow(x,4)*cospz*sint*pow(M_PI,2)/pow(r,4)
					+Real(12)/5*sinpr*pow(y,2)*cospz*sint*cospr*M_PI/pow(r,3)
					-Real( 2)/1*sinpr*pow(y,2)*cospz*sint*cospr*M_PI*pow(x,2)/pow(r,5)
					+Real( 2)/5*pow(cospr,2)*pow(M_PI,2)*pow(x,2)*pow(y,2)*cospz*sint/pow(r,4)
					-Real( 2)/5*pow(sinpr,2)*cospz*sint/r2
					+Real(12)/5*pow(sinpr,2)*pow(x,2)*cospz*sint/pow(r,4)
					-Real( 8)/5*pow(sinpr,2)*pow(x,4)*cospz*sint/pow(r,6)
					-Real( 2)/5*pow(sinpr,2)*pow(y,2)*cospz*sint*pow(M_PI,2)*pow(x,2)/pow(r,4)
					-Real(12)/5*pow(sinpr,2)*pow(y,2)*cospz*sint/pow(r,4)
					+Real( 8)/5*pow(sinpr,2)*pow(y,2)*cospz*sint*pow(x,2)/pow(r,6);
		}
};

class sourceZ3 : public basics::function3D {
	public:
		Real val(Real x, Real y, Real z, Real t) const
		{
			if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
				return 0;

			SETUPLOCALS
			
			sinpz = sin(M_PI*(z+1)/2);
			cospz = cos(M_PI*(z+1)/2);

			return ( Real(2)/5*sinpr*pow(M_PI,2)*pow(x,2)*(2*M_PI*cospr+2*sinpr/r)*y*sinpz/(M_PI*pow(r,3))
					+Real(6)/5*cospr*M_PI*pow(x,2)*(2*M_PI*cospr+2*sinpr/r)*y*sinpz/(M_PI*pow(r,4))
					-Real(8)/5*cospr*M_PI*(2*M_PI*cospr+2*sinpr/r)*y*sinpz/(M_PI*pow(r,4))
					-Real(4)/5*cospr*M_PI*x*(-2*M_PI*sinpr*M_PI*x/r-2*sinpr*x/pow(r,3)+2*cospr*M_PI*x/r2)*y*sinpz/(M_PI*r2)
					-Real(2)/(5*r*M_PI)*(sinpr*(-2*M_PI*cospr*pow(M_PI,2)*pow(x,2)/r2+2*M_PI*sinpr*M_PI*pow(x,2)/pow(r,3)-2*M_PI*sinpr*M_PI/r+6*sinpr*pow(x,2)/pow(r,5)-6*cospr*M_PI*pow(x,2)/pow(r,4)-2*sinpr/pow(r,3)-2*sinpr*pow(M_PI,2)*pow(x,2)/pow(r,3)+2*cospr*M_PI*x/r2)*y*sinpz)
					+Real(4)/5*sinpr*(-2*M_PI*sinpr*M_PI*x/r-2*sinpr*x/pow(r,3)+2*cospr*M_PI*x/r2)*y*sinpz*x/(M_PI*pow(r,3))
					-Real(6)/5*sinpr*(2*M_PI*cospr+2*sinpr/r)*y*sinpz*pow(x,2)/(M_PI*pow(r,5))
					+Real(8)/5*sinpr*(2*M_PI*cospr+2*sinpr/r)*y*sinpz/(M_PI*pow(r,3))
					+Real(2)/5*sinpr*pow(M_PI,2)*pow(y,3)*(2*M_PI*cospr+2*sinpr/r)*sinpz/(M_PI*pow(r,3))
					+Real(6)/5*cospr*M_PI*pow(y,3)*(2*M_PI*cospr+2*sinpr/r)*sinpz/(M_PI*pow(r,4))
					-Real(4)/5*cospr*M_PI*pow(y,2)*(-2*M_PI*sinpr*M_PI*y/r-2*sinpr*y/pow(r,3)+2*cospr*M_PI*y/r2)*sinpz/(M_PI*r2)
					-Real(2)/(5*M_PI*r)*(sinpr*(-2*M_PI*cospr*pow(M_PI,2)*pow(y,2)/r2+2*M_PI*sinpr*M_PI*pow(y,2)/pow(r,3)-2*M_PI*sinpr*M_PI/r+6*sinpr*pow(y,2)/pow(r,5)-6*cospr*M_PI*pow(y,2)/pow(r,4)-2*sinpr/pow(r,3)-2*sinpr*pow(M_PI,2)*pow(y,2)/pow(r,3)+2*cospr*M_PI/r2)*y*sinpz)
					-Real(4)/5*sinpr*(-2*M_PI*sinpr*M_PI*y/r-2*sinpr*y/pow(r,3)+2*cospr*M_PI*y/r2)*sinpz/(M_PI*r)
					+Real(4)/5*sinpr*(-2*M_PI*sinpr*M_PI*y/r-2*sinpr*y/pow(r,3)+2*cospr*M_PI*y/r2)*pow(y,2)*sinpz/(M_PI*pow(r,3))
					-Real(6)/5*sinpr*(2*M_PI*cospr+2*sinpr/r)*pow(y,3)*sinpz/(M_PI*pow(r,5)))*sint
					+Real(1)/10*sinpr*(2*M_PI*cospr+2*sinpr/r)*y*sinpz*M_PI/r+pow(sinpr,2)*cospz*M_PI*sint;
					+Real(2)/5*sinpr*(2*M_PI*cospr+2*sinpr/r)*y*sinpz*cost/(M_PI*r);
		}
};

class solX3 : public basics::function3D {
	public:
		Real val(Real x, Real y, Real z, Real t) const
		{
			if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
				return 0;
			SETUPLOCALS
			
			cospz = cos(M_PI*(z+1)/2);

			return 	-Real(1)/5*pow(sinpr,2)*y*cospz*sint*x/r2
					+Real(1)/5*pow(sinpr,2)*pow(x,2)*cospz*sint/r2;
		}

		Real diffx(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS
			cospz = cos(M_PI*(z+1)/2);

			return	-Real(4)/5*sinpr*y*cospz*sint*pow(x,2)*cospr*M_PI/pow(r,3)
					+Real(4)/5*pow(sinpr,2)*cospz*sint*pow(x,2)/pow(r,4)
					-Real(2)/5*pow(sinpr,2)*y*cospz*sint/r2;
		}

		Real diffy(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS
			
			cospz = cos(M_PI*(z+1)/2);
			return	-Real(4)/5*sinpr*pow(y,2)*cospz*sint*x*cospr*M_PI/pow(r,3)
					-Real(2)/5*pow(sinpr,2)*cospz*sint*x/r2
					+Real(4)/5*pow(sinpr,2)*x*cospz*sint*pow(y,2)/pow(r,4);
		}

		Real diffz(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS
			
			sinpz = sin(M_PI*(z+1)/2);
			return	Real(1)/5*pow(sinpr,2)*y*sinpz*M_PI*sint*x/r2;
		}

		Real difft(Real x, Real y, Real z, Real t) const
		{
			if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
				return 0;
			SETUPLOCALS
			
			cospz = cos(M_PI*(z+1)/2);

			return 	-Real(1)/5*pow(sinpr,2)*y*cospz*cost*x/r2
					+Real(1)/5*pow(sinpr,2)*pow(x,2)*cospz*cost/r2;
			
		}

		Real diff2x(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS
			
			cospz = cos(M_PI*(z+1)/2);

			return	-Real( 4)/5*cos(2*M_PI*r)*pow(M_PI,2)*pow(x,3)*y*cospz*sint/pow(r,4)
					+Real(12)/5*pow(sinpr,2)*cospz*sint*y*x/pow(r,4)
					+Real( 4)/1*sinpr*y*cospz*sint*pow(x,3)*cospr*M_PI/pow(r,5)
					-Real(12)/5*sinpr*x*cospz*sint*y*cospr*M_PI/pow(r,3);
					-Real(16)/5*pow(sinpr,2)*y*cospz*sint*pow(x,3)/pow(r,6);
		}

		Real diff2y(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS
			
			cospz = cos(M_PI*(z+1)/2);
	
			return	-Real( 4)/5*cos(2*M_PI*r)*pow(M_PI,2)*pow(y,3)*cospz*sint*x/pow(r,4)
					-Real(12)/5*sinpr*y*cospz*sint*x*cospr*M_PI/pow(r,3)
					+Real( 4)/1*sinpr*pow(y,3)*cospz*sint*x*cospr*M_PI/pow(r,5)
					+Real(12)/5*pow(sinpr,2)*y*cospz*sint*x/pow(r,4)
					-Real(16)/5*pow(sinpr,2)*x*cospz*sint*pow(y,3)/pow(r,6);
		}

		Real diff2z(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS
			
			cospz = cos(M_PI*(z+1)/2);
	
			return	Real(1)/10*pow(sinpr,2)*y*cospz*pow(M_PI,2)*sint*x/r2;
		}
};

class solY3 : public basics::function3D {
	public:
		Real val(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS
			
			cospz = cos(M_PI*(z+1)/2);

			return 	-Real(1)/5*pow(sinpr,2)*pow(y,2)*cospz*sint/r2
					+Real(1)/5*pow(sinpr,2)*pow(x,2)*cospz*sint/r2;
		}

		Real diffx(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS
			
			cospz = cos(M_PI*(z+1)/2);
		
			return	-Real(2)/5*sinpr*pow(y,2)*cospz*sint*x*cospr*M_PI/pow(r,3)
					+Real(2)/5*pow(sinpr*y,2)*cospz*sint*x/pow(r,4)
					+Real(2)/5*sinpr*pow(x,3)*cospz*sint*cospr*M_PI/pow(r,3)
					+Real(2)/5*pow(sinpr,2)*x*cospz*sint/r2
					-Real(2)/5*pow(sinpr,2)*pow(x,3)*cospz*sint/pow(r,4);
		}

		Real diffy(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS
			
			cospz = cos(M_PI*(z+1)/2);
		
			return	-Real(2)/5*sinpr*pow(y,3)*cospz*sint*cospr*M_PI/pow(r,3)
					-Real(2)/5*pow(sinpr,2)*y*cospz*sint/r2
					+Real(2)/5*pow(sinpr,2)*pow(y,3)*cospz*sint/pow(r,4)
					+Real(2)/5*sinpr*pow(x,2)*cospz*sint*y*cospr*M_PI/pow(r,3)
					-Real(2)/5*pow(sinpr*x,2)*cospz*sint*y/pow(r,4);
		}

		Real diffz(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS
			
			sinpz = sin(M_PI*(z+1)/2);
			return	 Real(1)/10*pow(sinpr*y,2)*sinpz*M_PI*sint/r2
					-Real(1)/10*pow(sinpr*x,2)*sinpz*M_PI*sint/r2;
		}

		Real difft(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS
			
			cospz = cos(M_PI*(z+1)/2);

			return 	-Real(1)/5*pow(sinpr,2)*pow(y,2)*cospz*cost/r2
					+Real(1)/5*pow(sinpr,2)*pow(x,2)*cospz*cost/r2;
		}

		Real diff2x(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS
			
			cospz = cos(M_PI*(z+1)/2);
		
			return	-Real(2)/5*cos(2*M_PI*r)*pow(M_PI*x*y,2)*cospz*sint/pow(r,4)
					+Real(2)/1*sinpr*pow(x*y,2)*cospz*sint*cospr*M_PI/pow(r,5)
					-Real(2)/5*sinpr*pow(y,2)*cospz*sint*cospr*M_PI/pow(r,3)
					-Real(8)/5*pow(sinpr*y*x,2)*cospz*sint/pow(r,6)
					+Real(2)/5*pow(sinpr*x,2)*cospz*sint/pow(r,4)
					+Real(2)/5*pow(cospr*M_PI*x*x,2)*cospz*sint/pow(r,4)
					+Real(2)/1*sinpr*pow(x,2)*cospz*sint*cospr*M_PI/pow(r,3)
					-Real(2)/1*sinpr*pow(x,4)*cospz*sint*cospr*M_PI/pow(r,5)
					-Real(2)/5*pow(sinpr*x*x*M_PI,2)*cospz*sint/pow(r,4)
					+Real(2)/5*pow(sinpr,2)*cospz*sint/r2
					-Real(2)/1*pow(sinpr*x,2)*cospz*sint/pow(r,4)
					+Real(8)/5*pow(sinpr*x*x,2)*cospz*sint/pow(r,6);
		}
		
		Real diff2y(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS
			
			cospz = cos(M_PI*(z+1)/2);

			return	-Real(2)/5*pow(cospr*M_PI*y*y,2)*cospz*sint/pow(r,4)
					-Real(2)/1*sinpr*pow(y,2)*cospz*sint*cospr*M_PI/pow(r,3)
					+Real(2)/1*sinpr*pow(y,4)*cospz*sint*cospr*M_PI/pow(r,5)
					+Real(2)/5*pow(sinpr*M_PI*y*y,2)*cospz*sint/pow(r,4)
					-Real(2)/5*pow(sinpr,2)*cospz*sint/r2
					+Real(2)/1*pow(sinpr*y,2)*cospz*sint/pow(r,4)
					-Real(8)/5*pow(sinpr*y*y,2)*cospz*sint/pow(r,6)
					+Real(2)/5*pow(cospr*M_PI*x*y,2)*cospz*sint/pow(r,4)
					-Real(2)/1*sinpr*pow(x*y,2)*cospz*sint*cospr*M_PI/pow(r,5)
					+Real(2)/5*sinpr*pow(x,2)*cospz*sint*cospr*M_PI/pow(r,3)
					-Real(2)/5*pow(sinpr*x*y*M_PI,2)*cospz*sint/pow(r,4)
					+Real(8)/5*pow(sinpr*x*y,2)*cospz*sint/pow(r,6)
					-Real(2)/5*pow(sinpr*x,2)*cospz*sint/pow(r,4);
		}

		Real diff2z(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS
			
			sinpz = sin(M_PI*(z+1)/2);
	
			return 	 Real(1)/20*pow(sinpr*y*M_PI,2)*cospz*sint/r2
					-Real(1)/20*pow(sinpr*x*M_PI,2)*cospz*sint/r2;
		}
};

class solZ3 : public basics::function3D {
	public:
		Real val(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS
			
			sinpz = sin(M_PI*(z+1)/2);
			
			return Real(2)/5*sinpr*(2*M_PI*cospr+2*sinpr/r)*y*sinpz*sint/(M_PI*r);
		}

		Real diffx(Real x, Real y, Real z, Real t) const
		{
			SETUPLOCALS
			
			sinpz = sin(M_PI*(z+1)/2);
	
			return	 Real(2)/5*cospr*M_PI*x*(2*M_PI*cospr+2*sinpr/r)*y*sinpz*sint/(M_PI*r2)
					+Real(2)/5*sinpr*(-2*M_PI*sinpr*M_PI*x/r-2*sinpr*x/pow(r,3)+2*cospr*M_PI*x/r2)*y*sinpz*sint/(M_PI*r)
					-Real(2)/5*sinpr*(2*M_PI*cospr+2*sinpr/r)*y*sinpz*sint*x/(M_PI*pow(r,3));
		}
};

class solP3 : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return sinpr*sinpr*sinpz*sint;
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
	Real T = 1.f;
	if( argc > 3 )
		Dt = Real(atoi(argv[2]))/atoi(argv[3]);
	if( argc > 5 )
		T = Real(atoi(argv[4]))/atoi(argv[5]);
	Real Lz=1;
	bool lumped=true;
	if( argc > 6 )
		lumped = atoi(argv[6])==1?false:true;
	bool secondorder=true;
	bool preconditioned=true;
	bool rotational=false;
	bool p3D=false;
//    if( !loaded) {
		if( argc > 5 )
			T = Real(atoi(argv[4]))/atoi(argv[5]);
		if( argc > 7 )
			secondorder = atoi(argv[7])==1?false:true;
		
		if( argc > 8 )
			preconditioned = atoi(argv[8])==1?false:true;

		if( argc > 9 )
			rotational = atoi(argv[9])==1?true:false;

		if( argc > 10 )
			p3D = atoi(argv[10])==1?true:false;
//    }

	legendreLegendreW::poissonSolver SP(N,N,0,legendreLegendreW::poissonSolver::Homogenous,true,Lz);
	legendreLegendreW::poissonSolver SP2(N,4,0,legendreLegendreW::poissonSolver::Nonhomogenous,true,Lz);
//    legendreLegendreW::poissonSolver SP3(N,N,0,legendreLegendreW::poissonSolver::HomogenousLeft,true,Lz);

	basics::Vector gridGL("GL grid",N-1);
	basics::Vector grid("GL grid",N+1);
	basics::Vector weightGL("GL weight",N-1);
	utilities::GLL::GaussLobattoLegendreGrid(grid);
	utilities::GLL::GaussLegendreGrid(gridGL);
	utilities::GLL::GaussLegendreWeights(weightGL,gridGL);
	basics::Matrix GLL2G = utilities::GLL::interpolationMatrix(grid,gridGL);

	bigCircleGeometry circle(N,N,grid,SP.weightX());
//    quadraticGeometry circle(N,N,grid,SP.weightX());
	circle.setLz(Lz);

	Real nu;
	if( secondorder )
		nu = Real(3)/(2*Dt);
	else
		nu = Real(1)/Dt;

	extrudedConsistentPressureEvaluator pressure3D(circle,SP.Dx(),GLL2G,weightGL,SP.weightX(),
											 	   grid,gridGL,preconditioned,rank,size,
											 	   SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
	SEMLaplacianEvaluator eval3(circle,SP.Dx(),SP.weightX(),NULL,0,1,
								SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
	extrudedLaplacianEvaluator velocity3D1(eval3,SP,nu,&SP2.Ax(),preconditioned,rank,size,
										  SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
	extrudedLaplacianEvaluator velocity3D2(eval3,SP,nu,&SP2.Ax(),preconditioned,rank,size,
										  SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
	extrudedLaplacianEvaluator velocity3D3(eval3,SP,nu,&SP2.Ax(),preconditioned,rank,size,
										  SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
//    extrudedLaplacianEvaluator velocity3DN(eval3,SP3,Real(3)/(2*Dt),&SP2.Ax(),preconditioned,rank,size,
//                                           SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM);
	extrudedFEMConsistentPressurePreconditioner eval4(circle,SP.weightX(),weightGL,grid,gridGL,rank,size);

    basics::coarseGrid desc = circle.getDivisionInfo(size);
	u.add(*utilities::g_manager.aquireMatricesStackField("velocity 1",N+1,N+1,N+1,desc[rank].elements.size()));
	u.add(*utilities::g_manager.aquireMatricesStackField("velocity 2",N+1,N+1,N+1,desc[rank].elements.size()));
	u.add(*utilities::g_manager.aquireMatricesStackField("velocity 3",N+1,N+1,N+1,desc[rank].elements.size()));
	basics::matricesStack* ut = utilities::g_manager.aquireMatricesStack("buffer",N+1,N+1,N+1,desc[rank].elements.size());

	basics::Field3<basics::matricesStack>& F  = *utilities::g_manager.aquireMatricesStackField("source",N+1,N+1,N+1,desc[rank].elements.size());

	p.add(*utilities::g_manager.aquireMatricesStack("pressure 1",N-1,N-1,N-1,desc[rank].elements.size()));
	p.add(*utilities::g_manager.aquireMatricesStack("pressure 2",N-1,N-1,N-1,desc[rank].elements.size()));
	p.add(*utilities::g_manager.aquireMatricesStack("pressure 3",N-1,N-1,N-1,desc[rank].elements.size()));
	basics::matricesStack& pt = *utilities::g_manager.aquireMatricesStack("pressure temp",N-1,N-1,N-1,desc[rank].elements.size());
			
//    if( loaded )
//        loadState(u,p,n,Dt,rank,argv[2]);

	articleTest1X articleX;
	articleTest1Y articleY;
	articleTest1Z articleZ;
	articleTest1P articleP;
//    boxTest3DX articleX;
//    boxTest3DY articleY;
//    boxTest3DZ articleZ;
//    boxTest3DP articleP;
	source3 sourceX(articleX,articleP,0);
	source3 sourceY(articleY,articleP,1);
	source3 sourceZ(articleZ,articleP,2);

//    bugger bug;
//    circle.evaluate(p.get(0),gridGL,bug,0,desc[rank].elements);

	/* set up initial conditions */
	Real t0=0;
	circle.evaluate(u.get(0,-1),grid,articleX,articleY,articleZ,t0-Dt,desc[rank].elements);
	circle.evaluate(u.get(0, 0),grid,articleX,articleY,articleZ,t0,desc[rank].elements);
	circle.evaluate(p.get(0,-1),gridGL,articleP,t0-Dt,desc[rank].elements);
	circle.evaluate(p.get(0, 0),gridGL,articleP,t0,desc[rank].elements);

	if( rank == 0 ) {
		cout << "time step:\t\t"  	 << Dt << endl
			 <<	"final time:\t\t" 	 <<  T << endl
			 <<	"Lz:\t\t\t"  		 << Lz << endl
			 <<	"preconditioned:\t\t"<< (preconditioned?"yes":"no")  << endl
			 <<	"order :\t\t\t"	 	 << (secondorder?2:1) << endl
			 << "rotational:\t\t"  	 << (rotational?"yes":"no") << endl
			 << "lumped:\t\t\t"	 	 << (lumped?"yes":"no") << endl
			 << "3D evaluator:\t\t"  << (p3D?"yes":"no") << endl;
	}
		
	SEMInverseMassEvaluator invmass2D(circle,weightGL,GLL2G,rank,size);
	extrudedInverseMassEvaluator invmass3D(invmass2D,rank,size);

	MPIDotter pdotter(rank,size);

	/* setup keyboard */
	utilities::set_raw_tty();

	int elipticu = utilities::g_profiler.add("eliptic velocity");
	int elipticp = utilities::g_profiler.add("eliptic pressure");
	for( n;n<(T-t0)/Dt;++n ) {
		if( rank == 0 )
			std::cout << "Time: " << t0+(n+1)*Dt << std::endl;

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

//        F.clear();
//        F.Z() = -1;
		circle.evaluate(F,grid,sourceX,sourceY,sourceZ,t0+(n+1)*Dt,desc[rank].elements);

		if( secondorder ) {
			/* BDF2 */
			F.axpy(Real(2)/Dt,un);
			F.axpy(Real(-1)/(2*Dt),unm1);
		} else {
			/* BDF1*/
			F.axpy(Real(1)/Dt,un);
		}

		circle.mass(F,desc[rank].elements);

		/* extrapolate pressure */
		if( secondorder )
			utilities::extrapolate(pnp1,p,n,utilities::FIRST_ORDER);
		else
			pnp1 = pn;

		pressure3D.m_gradient.evaluate(unp1,pnp1);
		unp1 += F;

		/* delta */
		velocity3D1.evaluate(unm1,un,false,false);
		unp1 -= unm1;

		mask(unp1,circle,rank,size,velocity3D1.m_bc);
		F = unp1;
		dssum(unp1,circle,rank,size);
		utilities::g_profiler.start(elipticu);
		if( lumped )
			velocity3D1.solve(unp1,F,velocity3D1,velocity3D2,velocity3D3);
		else {
			velocity3D1.m_name = "velocity X";
			velocity3D1.solve(unp1.X(),F.X(),0);
			velocity3D1.m_name = "velocity Y";
			velocity3D1.solve(unp1.Y(),F.Y(),1);
			velocity3D1.m_name = "velocity Z";
			velocity3D1.solve(unp1.Z(),F.Z(),2);
		}
		utilities::g_profiler.pause(elipticu);
		unp1 += un;

		mask(unp1,circle,rank,size,velocity3D1.m_bc);
		pressure3D.m_divergence.evaluate(pt,unp1);
		invmass3D.evaluate(pnm1,pt);
		pt *= -nu;

		utilities::g_profiler.start(elipticp);
		if( p3D ) {
//            if( p3D == 2 )
//                eval.solve3(pt2,pt,eval.m_pre3);
//            else {
			{
				if( preconditioned )
					cout << "iterations pressure " << utilities::CGSolver::solve(pt,pressure3D,eval4,pdotter,1.e-10) << endl;
				else
					cout << "iterations pressure " << utilities::CGSolver::solve(pt,pressure3D,pdotter,1.e-10) << endl;
			}
		} else 
			pressure3D.solve(pt);

		utilities::g_profiler.pause(elipticp);

		pnp1 += pt;

		/* update velocity */
		pressure3D.m_gradient.evaluate(unm1,pt);
		dssum(unm1,circle,rank,size);
		circle.invMass(unm1,desc[rank].elements);
		mask(unm1,circle,rank,size,velocity3D1.m_bc);
		unp1.axpy(Real(1)/nu,unm1);
		
//        if( rotational )
//            pnp1 -= pnm1;
	}

	/* restore keyboard */
	utilities::set_normal_tty();

//    saveState(u,p,n,Dt,rank);

	velocity3D1.m_nu = 1;
	velocity3D1.m_nosum = &F.X();
	string err = errorReport(u.get(n),p.get(n),circle,articleX,articleY,articleZ,
							 articleP,grid,gridGL,t0+n*Dt,invmass3D,&velocity3D1);

	if( rank == 0 ) {
		cout << err;
		cout << utilities::g_profiler.report();
	}
#ifdef HAS_MPI
	MPI_Finalize();
#endif

	return 0; // AY-OH-KAY
}

