#ifndef FUNCTIONS_3D_H_
#define FUNCTIONS_3D_H_

#include <limits>

class evilStationaryTest3DX : public basics::function3D {
	public:
		Real val(Real x, Real y, Real z, Real t) const
		{
			return( sin(x)*sin(x)*sin(Real(2)*M_PI*y)*sin(z) );
		}

		Real diffx(Real x, Real y, Real z, Real t) const
		{
			return( sin(Real(2)*x)*sin(Real(2)*M_PI*y)*sin(z) );
		}

		Real diffy(Real x, Real y, Real z, Real t) const
		{
			return( Real(2)*M_PI*sin(x)*sin(x)*cos(Real(2)*M_PI*y)*sin(z) );
		}

		Real diffz(Real x, Real y, Real z, Real t) const
		{
			return( sin(x)*sin(x)*sin(Real(2)*M_PI*y)*cos(z) );
		}

		Real difft(Real x, Real y, Real z, Real t) const
		{
			return( 0 );
		}

		Real diff2x(Real x, Real y, Real z, Real t) const
		{
			return( Real(2)*cos(Real(2)*x)*sin(Real(2)*M_PI*y)*sin(z) );
		}

		Real diff2y(Real x, Real y, Real z, Real t) const
		{
			return( -Real(4)*M_PI*M_PI*sin(x)*sin(x)*sin(Real(2)*M_PI*y)*sin(z) );
		}

		Real diff2z(Real x, Real y, Real z, Real t) const
		{
			return( -sin(x)*sin(x)*sin(Real(2)*M_PI*y)*sin(z) );
		}

		static int id() 
		{
			return( 1 );
		}
};

class evilStationaryTest3DY : public basics::function3D {
	
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( sin(Real(2)*x)*sin(M_PI*y)*sin(M_PI*y)*cos(z)/M_PI );
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		return( Real(2)*cos(Real(2)*x)*sin(M_PI*y)*sin(M_PI*y)*cos(z)/M_PI );
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		return( sin(Real(2)*x)*sin(Real(2)*M_PI*y)*cos(z) );
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		return( -sin(Real(2)*x)*sin(M_PI*y)*sin(M_PI*y)*sin(z)/M_PI );
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		return( 0 );
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		return( -Real(4)*sin(Real(2)*x)*sin(M_PI*y)*sin(M_PI*y)*cos(z)/M_PI );
	}
	
	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		return( Real(2)*M_PI*sin(Real(2)*x)*cos(Real(2)*M_PI*y)*cos(z) );
	}
	
	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		return( -sin(Real(2)*x)*sin(M_PI*y)*sin(M_PI*y)*cos(z)/M_PI );
	}

	static int id() 
	{
		return( 1 );
	}
};


class evilStationaryTest3DZ : public basics::function3D {
	
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( sin(Real(2)*x)*sin(Real(2)*M_PI*y)*(cos(z)-sin(z)) );
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		return( Real(2)*cos(Real(2)*x)*sin(Real(2)*M_PI*y)*(cos(z)-sin(z)) );
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		return( Real(2)*M_PI*sin(Real(2)*x)*cos(Real(2)*M_PI*y)*(cos(z)-sin(z)) );
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		return( sin(Real(2)*x)*sin(Real(2)*M_PI*y)*(-sin(z)-cos(z)) );
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		return( 0 );
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		return( -Real(4)*sin(Real(2)*x)*sin(Real(2)*M_PI*y)*(cos(z)-sin(z)) );
	}
	
	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		return( -Real(4)*M_PI*M_PI*sin(Real(2)*x)*sin(Real(2)*M_PI*y)*(cos(z)-sin(z)) );
	}
	
	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		return( sin(Real(2)*x)*sin(Real(2)*M_PI*y)*(-cos(z)+sin(z)) );
	}

	static int id() 
	{
		return( 1 );
	}
};

class evilStationaryTest3DP : public basics::function3D {
	
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( 0 );
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
		return( 0 );
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

	static int id() 
	{
		return( 1 );
	}
};

class evilTest3DX : public basics::function3D {
	public:
		Real val(Real x, Real y, Real z, Real t) const
		{
			return( (Real(4)/5)*M_PI*pow(sin(M_PI*x),2)*
					 cos(M_PI*y)*sin(M_PI*y)*cos(M_PI*z)*sin(M_PI*z) );
		}

		Real diffx(Real x, Real y, Real z, Real t) const
		{
			return( (Real(8)/5)*pow(M_PI,2)*sin(M_PI*x)*cos(M_PI*y)*
					sin(M_PI*y)*cos(M_PI*z)*sin(M_PI*z)*cos(M_PI*x) );
		}

		Real diffy(Real x, Real y, Real z, Real t) const
		{
			return( -(Real(4)/5)*pow(M_PI,2)*(-1+pow(cos(M_PI*x),2))*cos(M_PI*z)*sin(M_PI*z)*(-1+2*pow(cos(M_PI*y),2)) );
		}

		Real diffz(Real x, Real y, Real z, Real t) const
		{
			return( -(Real(4)/5)*pow(M_PI,2)*pow(sin(M_PI*x),2)*cos(M_PI*y)*
					sin(M_PI*y)*pow(sin(M_PI*z),2)+(Real(4)/5)*pow(M_PI,2)*
					pow(sin(M_PI*x),2)*cos(M_PI*y)*sin(M_PI*y)*pow(cos(M_PI*z),2) );
		}

		Real difft(Real x, Real y, Real z, Real t) const
		{
			return( 0 );
		}

		Real diff2x(Real x, Real y, Real z, Real t) const
		{
			return( (Real(8)/5)*pow(M_PI,3)*pow(cos(M_PI*x),2)*cos(M_PI*y)*
					sin(M_PI*y)*cos(M_PI*z)*sin(M_PI*z)-(Real(8)/5)*pow(M_PI,3)*
					pow(sin(M_PI*x),2)*cos(M_PI*y)*sin(M_PI*y)*cos(M_PI*z)*sin(M_PI*z) );
		}

		Real diff2y(Real x, Real y, Real z, Real t) const
		{
			return( -(Real(16)/5)*pow(M_PI,3)*pow(sin(M_PI*x),2)*cos(M_PI*y)*sin(M_PI*y)*
					 cos(M_PI*z)*sin(M_PI*z) );
		}

		Real diff2z(Real x, Real y, Real z, Real t) const
		{
			return( -(Real(16)/5)*pow(M_PI,3)*pow(sin(M_PI*x),2)*cos(M_PI*y)*
					sin(M_PI*y)*cos(M_PI*z)*sin(M_PI*z) );
		}

		static int id() 
		{
			return( 1 );
		}
};

class evilTest3DY : public basics::function3D {
	
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( -(Real(2)/5)*M_PI*cos(M_PI*x)*sin(M_PI*x)*pow(sin(M_PI*y),2)*
				cos(M_PI*z)*sin(M_PI*z) );
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		return( (Real(2)/5)*pow(M_PI,2)*pow(sin(M_PI*x),2)*pow(sin(M_PI*y),2)*
				cos(M_PI*z)*sin(M_PI*z)-(Real(2)/5)*pow(M_PI,2)*pow(cos(M_PI*x),2)*
				pow(sin(M_PI*y),2)*cos(M_PI*z)*sin(M_PI*z) );
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		return( -(Real(4)/5)*pow(M_PI,2)*cos(M_PI*x)*sin(M_PI*x)*
				sin(M_PI*y)*cos(M_PI*z)*sin(M_PI*z)*cos(M_PI*y) );
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		return( (Real(2)/5)*pow(M_PI,2)*cos(M_PI*x)*sin(M_PI*x)*pow(sin(M_PI*y),2)*
				pow(sin(M_PI*z),2)-(Real(2)/5)*pow(M_PI,2)*cos(M_PI*x)*sin(M_PI*x)*
				pow(sin(M_PI*y),2)*pow(cos(M_PI*z),2) );
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		return( 0 );
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		return( (Real(8)/5)*pow(M_PI,3)*sin(M_PI*x)*pow(sin(M_PI*y),2)*cos(M_PI*z)*
				sin(M_PI*z)*cos(M_PI*x) );
	}
	
	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		return( -(Real(4)/5)*pow(M_PI,3)*cos(M_PI*x)*sin(M_PI*x)*pow(cos(M_PI*y),2)*
				cos(M_PI*z)*sin(M_PI*z)+(Real(4)/5)*pow(M_PI,3)*sin(M_PI*x)*
				pow(sin(M_PI*y),2)*cos(M_PI*z)*sin(M_PI*z)*cos(M_PI*x) );
	}
	
	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		return( (Real(8)/5)*pow(M_PI,3)*sin(M_PI*x)*pow(sin(M_PI*y),2)*
				cos(M_PI*z)*sin(M_PI*z)*cos(M_PI*x) );
	}

	static int id() 
	{
		return( 1 );
	}
};


class evilTest3DZ : public basics::function3D {
	
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( -(Real(2)/5)*M_PI*cos(M_PI*x)*sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*y)*
				pow(sin(M_PI*z),2) );
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		return( (Real(2)/5)*pow(M_PI,2)*pow(sin(M_PI*x),2)*cos(M_PI*y)*sin(M_PI*y)*
				pow(sin(M_PI*z),2)-(Real(2)/5)*pow(M_PI,2)*pow(cos(M_PI*x),2)*cos(M_PI*y)*
				sin(M_PI*y)*pow(sin(M_PI*z),2) );
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		return( (Real(2)/5)*pow(M_PI,2)*cos(M_PI*x)*sin(M_PI*x)*
				pow(sin(M_PI*y),2)*pow(sin(M_PI*z),2)-(Real(2)/5)*pow(M_PI,2)*
				cos(M_PI*x)*sin(M_PI*x)*pow(cos(M_PI*y),2)*pow(sin(M_PI*z),2) );
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		return( -(Real(4)/5)*pow(M_PI,2)*cos(M_PI*x)*sin(M_PI*x)*cos(M_PI*y)*
				sin(M_PI*y)*sin(M_PI*z)*cos(M_PI*z) );
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		return( 0 );
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		return( (Real(8)/5)*pow(M_PI,3)*sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*y)*
				pow(sin(M_PI*z),2)*cos(M_PI*x) );
	}
	
	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		return( (Real(8)/5)*pow(M_PI,3)*sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*y)*
				pow(sin(M_PI*z),2)*cos(M_PI*x) );
	}
	
	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		return( -(Real(4)/5)*pow(M_PI,3)*cos(M_PI*x)*sin(M_PI*x)*cos(M_PI*y)*
				sin(M_PI*y)*pow(cos(M_PI*z),2)+(Real(4)/5)*pow(M_PI,3)*sin(M_PI*x)*
				cos(M_PI*y)*sin(M_PI*y)*pow(sin(M_PI*z),2)*cos(M_PI*x) );
	}

	static int id() 
	{
		return( 1 );
	}
};

class evilTest3DP : public basics::function3D {
	
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) );
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		return( -sin(M_PI*x)*M_PI*sin(M_PI*y)*sin(M_PI*z) );
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		return( cos(M_PI*x)*cos(M_PI*y)*M_PI*sin(M_PI*z) );
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		return( cos(M_PI*x)*sin(M_PI*y)*cos(M_PI*z)*M_PI );
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		return( 0 );
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		return 0;
	}
	
	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		return 0;
	}
	
	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		return 0;
	}

	static int id() 
	{
		return( 1 );
	}
};


class rappTest3DX : public basics::function3D {
	public:
		Real val(Real x, Real y, Real z, Real t) const
		{
			return( 1.f-y*y+sin(x)*sin(x)*sin(Real(2)*M_PI*y)*sin(z)*sin(t) );
		}

		Real diffx(Real x, Real y, Real z, Real t) const
		{
			return( sin(Real(2)*x)*sin(Real(2)*M_PI*y)*sin(z)*sin(t) );
		}

		Real diffy(Real x, Real y, Real z, Real t) const
		{
			return( -2*y + Real(2)*M_PI*sin(x)*sin(x)*cos(Real(2)*M_PI*y)*sin(z)*sin(t));
		}

		Real diffz(Real x, Real y, Real z, Real t) const
		{
			return( sin(x)*sin(x)*sin(Real(2)*M_PI*y)*cos(z)*sin(t) );
		}

		Real difft(Real x, Real y, Real z, Real t) const
		{
			return( sin(x)*sin(x)*sin(Real(2)*M_PI*y)*sin(z)*cos(t) );
		}

		Real diff2x(Real x, Real y, Real z, Real t) const
		{
			return( Real(2)*cos(Real(2)*x)*sin(Real(2)*M_PI*y)*sin(z)*sin(t) );
		}

		Real diff2y(Real x, Real y, Real z, Real t) const
		{
			return( -2 + -Real(4)*M_PI*M_PI*sin(x)*sin(x)*sin(Real(2)*M_PI*y)*sin(z)*sin(t));
		}

		Real diff2z(Real x, Real y, Real z, Real t) const
		{
			return( -sin(x)*sin(x)*sin(Real(2)*M_PI*y)*sin(z)*sin(t) );
		}

		static int id() 
		{
			return( 1 );
		}
};

class rappTest3DY : public basics::function3D {
	
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( sin(Real(2)*x)*sin(M_PI*y)*sin(M_PI*y)*cos(z)*sin(t)/M_PI );
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		return( Real(2)*cos(Real(2)*x)*sin(M_PI*y)*sin(M_PI*y)*cos(z)*sin(t)/M_PI );
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		return( sin(Real(2)*x)*sin(Real(2)*M_PI*y)*cos(z)*sin(t) );
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		return( -sin(Real(2)*x)*sin(M_PI*y)*sin(M_PI*y)*sin(z)*sin(t)/M_PI );
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		return( sin(Real(2)*x)*sin(M_PI*y)*sin(M_PI*y)*cos(z)*cos(t)/M_PI );
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		return( -Real(4)*sin(Real(2)*x)*sin(M_PI*y)*sin(M_PI*y)*cos(z)*sin(t)/M_PI );
	}
	
	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		return( Real(2)*M_PI*sin(Real(2)*x)*cos(Real(2)*M_PI*y)*cos(z)*sin(t) );
	}
	
	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		return( -sin(Real(2)*x)*sin(M_PI*y)*sin(M_PI*y)*cos(z)*sin(t)/M_PI );
	}

	static int id() 
	{
		return( 1 );
	}
};


class rappTest3DZ : public basics::function3D {
	
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( sin(Real(2)*x)*sin(Real(2)*M_PI*y)*(cos(z)-sin(z))*sin(t) );
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		return( Real(2)*cos(Real(2)*x)*sin(Real(2)*M_PI*y)*(cos(z)-sin(z))*sin(t) );
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		return( Real(2)*M_PI*sin(Real(2)*x)*cos(Real(2)*M_PI*y)*(cos(z)-sin(z))*sin(t) );
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		return( sin(Real(2)*x)*sin(Real(2)*M_PI*y)*(-sin(z)-cos(z))*sin(t) );
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		return( sin(Real(2)*x)*sin(Real(2)*M_PI*y)*(cos(z)-sin(z))*cos(t) );
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		return( -Real(4)*sin(Real(2)*x)*sin(Real(2)*M_PI*y)*(cos(z)-sin(z))*sin(t) );
	}
	
	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		return( -Real(4)*M_PI*M_PI*sin(Real(2)*x)*sin(Real(2)*M_PI*y)*(cos(z)-sin(z))*sin(t) );
	}
	
	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		return( sin(Real(2)*x)*sin(Real(2)*M_PI*y)*(-cos(z)+sin(z))*sin(t) );
	}

	static int id() 
	{
		return( 1 );
	}
};

class rappTest3DP : public basics::function3D {
	
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( sin(x)*cos(M_PI*y)*sin(z)*sin(z)*sin(t) );
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		return( cos(x)*cos(M_PI*y)*sin(z)*sin(z)*sin(t) );
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		return( -M_PI*sin(x)*sin(M_PI*y)*sin(z)*sin(z)*sin(t) );
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		return( sin(x)*cos(M_PI*y)*sin(2*z)*sin(t) );
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		return( sin(x)*cos(M_PI*y)*sin(z)*sin(z)*cos(t) );
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		return( -sin(x)*cos(M_PI*y)*sin(z)*sin(z)*sin(t) );
	}
	
	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		return( -M_PI*M_PI*sin(x)*cos(M_PI*y)*sin(z)*sin(z)*sin(t) );
	}
	
	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		return( 2*sin(x)*cos(M_PI*y)*cos(2*z)*sin(t) );
	}

	static int id() 
	{
		return( 1 );
	}
};

class channelTest3DX : public basics::function3D {
	public:
		Real val(Real x, Real y, Real z, Real t) const
		{
			return( 1.f-y*y );
		}

		Real diffx(Real x, Real y, Real z, Real t) const
		{
			return( 0 );
		}

		Real diffy(Real x, Real y, Real z, Real t) const
		{
			return( -2.f*y );
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
			return -2;
		}

		Real diff2z(Real x, Real y, Real z, Real t) const
		{
			return( 0 );
		}

		static int id() 
		{
			return( 2 );
		}
};

class channelTest3DY : public basics::function3D {
	
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( 0 );
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
		return( 0 );
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

	static int id() 
	{
		return( 2 );
	}
};


class channelTest3DZ : public basics::function3D {
	
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( 0 );
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
		return( 0 );
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

	static int id() 
	{
		return( 2 );
	}
};

class channelTest3DP : public basics::function3D {
public:
	channelTest3DP(Real n_nu)
	{
		nu = n_nu;
	}
	
	Real val(Real x, Real y, Real z, Real t) const
	{
//        return( 2.f*x*nu );
		return 0.f;
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		return( 2.f*nu );
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
		return( 0 );
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

	static int id() 
	{
		return( 2 );
	}
protected:
	Real nu;	
};

#define AMPLITUDE 0.05

class channelPerturbedTest3DX : public basics::function3D {
	public:
		Real val(Real x, Real y, Real z, Real t) const
		{
			return( 1.f-y*y+AMPLITUDE*(sin(x)*sin(x)*sin(Real(2)*M_PI*y)*sin(z)*sin(t)) );
		}

		Real diffx(Real x, Real y, Real z, Real t) const
		{
			return( 0 );
		}

		Real diffy(Real x, Real y, Real z, Real t) const
		{
			return( -2.f*y );
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

		static int id() 
		{
			return( 1 );
		}
};

class channelPerturbedTest3DY : public basics::function3D {
	
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( 0.f );//sin(Real(2)*x)*sin(M_PI*y)*sin(M_PI*y)*cos(z)*sin(t)/M_PI );
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
		return( 0 );
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

	static int id() 
	{
		return( 1 );
	}
};


class channelPerturbedTest3DZ : public basics::function3D {
	
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( 0.f );
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
		return( 0 );
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

	static int id() 
	{
		return( 1 );
	}
};

class channelPerturbedTest3DP : public basics::function3D {
	
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( 0 );
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
		return( 0 );
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

	static int id() 
	{
		return( 1 );
	}
};

class deformedTest3D : public basics::function3D {
	
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( (1-x*x)*(1-y*y)*(1-z*z) );
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		return( -2*x*(y*y-1)*(z*z-1) );
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		return( -2*y*(x*x-1)*(z*z-1) );
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		return( -2*z*(x*x-1)*(y*y-1) );
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		return( 0 );
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		return( -2*(y*y-1)*(z*z-1) );
	}
	
	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		return( -2*(x*x-1)*(z*z-1) );
	}
	
	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		return( -2*(x*x-1)*(y*y-1) );
	}

	static int id() 
	{
		return( 1 );
	}
};

#define OMEGA Real(1)

#define SETUPLOCALS \
		Real r2		= x*x+y*y; \
		Real r		= sqrt(r2); \
		Real cospr	= cos(OMEGA*M_PI*r); \
		Real sinpr	= sin(OMEGA*M_PI*r); \
		Real sinpz	= sin(2*M_PI*z); \
		Real cospz	= cos(2*M_PI*z); \
		Real sint	= sin(t); \
		Real cost	= cos(t);


class articleTest1X : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
			return 0;

		return( (Real(2)/5)*pow(sinpr,2)*y*sinpz*sint*x/r2 );
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return(  (Real(4)/5)*sinpr*y*sinpz*sint*pow(x,2)*cospr*OMEGA*M_PI/pow(r2,(Real(3)/2))
				-(Real(4)/5)*pow(sinpr,2)*y*sinpz*sint*pow(x,2)/pow(r2,2)
				+(Real(2)/5)*pow(sinpr,2)*y*sinpz*sint/r2 );
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return  (Real(4)/5)*sinpr*pow(y,2)*sinpz*sint*x*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
			   +(Real(2)/5)*pow(sinpr,2)*sinpz*sint*x/r2
			   -(Real(4)/5)*pow(sinpr,2)*pow(y,2)*sinpz*sint*x/pow(r2,2);
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return (Real(4)/5)*pow(sinpr,2)*y*cospz*M_PI*sint*x/r2;
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
			return 0;

		SETUPLOCALS

		return( (Real(2)/5)*pow(sinpr,2)*y*sinpz*cost*x/r2 );
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return(  (Real( 4)/5)*pow(cospr,2)*pow(OMEGA,2)*pow(M_PI,2)*pow(x,3)*y*sinpz*sint/pow(r2,2)
				-(Real( 4)/1)*sinpr*y*sinpz*sint*pow(x,3)*cospr*OMEGA*M_PI/pow(r2,Real(5)/2)
				+(Real(12)/5)*sinpr*y*sinpz*sint*x*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
				-(Real( 4)/5)*pow(sinpr,2)*y*sinpz*sint*pow(x,3)*pow(OMEGA,2)*pow(M_PI,2)/pow(r2,2)
				+(Real(16)/5)*pow(sinpr,2)*y*sinpz*sint*pow(x,3)/pow(r2,3)
				-(Real(12)/5)*pow(sinpr,2)*y*sinpz*sint*x/pow(r2,2) );
	}

	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return(  (Real( 4)/5)*pow(cospr,2)*pow(OMEGA,2)*pow(M_PI,2)*pow(y,3)*sinpz*sint*x/pow(r2,2)
				+(Real(12)/5)*sinpr*y*sinpz*sint*x*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
				-(Real( 4)/1)*sinpr*pow(y,3)*sinpz*sint*x*cospr*OMEGA*M_PI/pow(r2,Real(5)/2)
				-(Real( 4)/5)*pow(sinpr,2)*pow(y,3)*sinpz*sint*x*pow(OMEGA,2)*pow(M_PI,2)/pow(r2,2)
				-(Real(12)/5)*pow(sinpr,2)*y*sinpz*sint*x/pow(r2,2)
				+(Real(16)/5)*pow(sinpr,2)*pow(y,3)*sinpz*sint*x/pow(r2,3) );
	}

	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return( -(Real(8)/5)*pow(sinpr,2)*y*sinpz*pow(M_PI,2)*sint*x/r2 );
	}
};

class articleTest1Y : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
			return 0;

		SETUPLOCALS

		return(  (Real(1)/5)*pow(sinpr,2)*pow(y,2)*sinpz*sint/r2
				-(Real(1)/5)*pow(sinpr,2)*pow(x,2)*sinpz*sint/r2 );
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return  (Real(2)/5)*sinpr*pow(y,2)*sinpz*sint*cospr*OMEGA*M_PI*x/pow(r2,Real(3)/2)
			   -(Real(2)/5)*pow(sinpr,2)*pow(y,2)*sinpz*sint*x/pow(r2,2)
			   -(Real(2)/5)*sinpr*pow(x,3)*sinpz*sint*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
			   -(Real(2)/5)*pow(sinpr,2)*x*sinpz*sint/r2
			   +(Real(2)/5)*pow(sinpr,2)*pow(x,3)*sinpz*sint/pow(r2,2);
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return  (Real(2)/5)*sinpr*pow(y,3)*sinpz*sint*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
			   +(Real(2)/5)*pow(sinpr,2)*y*sinpz*sint/r2
			   -(Real(2)/5)*pow(sinpr,2)*pow(y,3)*sinpz*sint/pow(r2,2)
			   -(Real(2)/5)*sinpr*y*sinpz*sint*pow(x,2)*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
			   +(Real(2)/5)*pow(sinpr,2)*y*sinpz*sint*pow(x,2)/pow(r2,2);

	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return  (Real(2)/5)*pow(sinpr,2)*pow(y,2)*cospz*M_PI*sint/r2
			   -(Real(2)/5)*pow(sinpr,2)*pow(x,2)*cospz*M_PI*sint/r2;
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
			return 0;

		SETUPLOCALS

		return(  (Real(1)/5)*pow(sinpr,2)*pow(y,2)*sinpz*cost/r2
				-(Real(1)/5)*pow(sinpr,2)*pow(x,2)*sinpz*cost/r2 );
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return  (Real(2)/5)*pow(cospr,2)*pow(OMEGA,2)*pow(M_PI,2)*pow(x,2)*pow(y,2)*sinpz*sint/pow(r2,2)
			   -(Real(2)/1)*sinpr*pow(y,2)*sinpz*sint*pow(x,2)*cospr*OMEGA*M_PI/pow(r2,Real(5)/2)
			   +(Real(2)/5)*sinpr*pow(y,2)*sinpz*sint*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
			   -(Real(2)/5)*pow(sinpr,2)*pow(y,2)*sinpz*sint*pow(x,2)*pow(OMEGA,2)*pow(M_PI,2)/pow(r2,2)
			   +(Real(8)/5)*pow(sinpr,2)*pow(y,2)*sinpz*sint*pow(x,2)/pow(r2,3)
			   -(Real(2)/5)*pow(sinpr,2)*pow(y,2)*sinpz*sint/pow(r2,2)
			   -(Real(2)/5)*pow(cospr,2)*pow(OMEGA,2)*pow(M_PI,2)*pow(x,4)*sinpz*sint/pow(r2,2)
			   -(Real(2)/1)*sinpr*pow(x,2)*sinpz*sint*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
			   +(Real(2)/1)*sinpr*pow(x,4)*sinpz*sint*cospr*OMEGA*M_PI/pow(r2,Real(5)/2)
			   +(Real(2)/5)*pow(sinpr,2)*pow(x,4)*sinpz*sint*pow(OMEGA,2)*pow(M_PI,2)/pow(r2,2)
			   +(Real(2)/1)*pow(sinpr,2)*sinpz*sint*pow(x,2)/pow(r2,2)
			   -(Real(2)/5)*pow(sinpr,2)*sinpz*sint/r2
			   -(Real(8)/5)*pow(sinpr,2)*pow(x,4)*sinpz*sint/pow(r2,3);
	}

	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return  (Real(2)/5)*pow(cospr,2)*pow(OMEGA,2)*pow(M_PI,2)*pow(y,4)*sinpz*sint/pow(r2,2)
			   +(Real(2)/1)*sinpr*pow(y,2)*sinpz*sint*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
			   -(Real(2)/1)*sinpr*pow(y,4)*sinpz*sint*cospr*OMEGA*M_PI/pow(r2,Real(5)/2)
			   -(Real(2)/5)*pow(sinpr,2)*pow(y,4)*sinpz*sint*pow(OMEGA,2)*pow(M_PI,2)/pow(r2,2)
			   +(Real(2)/5)*pow(sinpr,2)*sinpz*sint/r2
			   -(Real(2)/1)*pow(sinpr,2)*pow(y,2)*sinpz*sint/pow(r2,2)
			   +(Real(8)/5)*pow(sinpr,2)*pow(y,4)*sinpz*sint/pow(r2,3)
			   -(Real(2)/5)*pow(cospr,2)*pow(OMEGA,2)*pow(M_PI,2)*pow(x,2)*pow(y,2)*sinpz*sint/pow(r2,2)
			   -(Real(2)/5)*sinpr*pow(x,2)*sinpz*sint*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
			   +(Real(2)/1)*sinpr*pow(y,2)*sinpz*sint*pow(x,2)*cospr*OMEGA*M_PI/pow(r2,Real(5)/2)
			   +(Real(2)/5)*pow(sinpr,2)*pow(y,2)*sinpz*sint*pow(x,2)*pow(OMEGA,2)*pow(M_PI,2)/pow(r2,2)
			   +(Real(2)/5)*pow(sinpr,2)*sinpz*sint*pow(x,2)/pow(r2,2)
			   -(Real(8)/5)*pow(sinpr,2)*pow(y,2)*sinpz*sint*pow(x,2)/pow(r2,3);
	}

	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return -(Real(4)/5)*pow(sinpr,2)*pow(y,2)*sinpz*pow(M_PI,2)*sint/r2
			   +(Real(4)/5)*pow(sinpr,2)*pow(x,2)*sinpz*pow(M_PI,2)*sint/r2;
	}
};

class articleTest1Z : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
			return 0;

		SETUPLOCALS

		return  (Real(1)/(10*r*M_PI))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz-1)*sint);
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return   (Real(1)/(10*r2))*(cospr*OMEGA*x*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz-1)*sint)
			    +(Real(1)/(10*M_PI*r))*(sinpr*(-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr*x/r-2*sinpr*x/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI*x/r2)*y*(cospz-1)*sint)
				-(Real(1)/(10*pow(r2,Real(3)/2)*M_PI))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz-1)*sint*x);
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return  (Real(1)/(10*r2))*(cospr*OMEGA*pow(y,2)*(2*M_PI*OMEGA*cospr+2*sinpr/r)*(cospz-1)*sint)
			   +(Real(1)/(10*r*M_PI))*(sinpr*(-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr*y/r-2*sinpr*y/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI*y/r2)*y*(cospz-1)*sint)
			   +(Real(1)/(10*r*M_PI))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*(cospz-1)*sint)
			   -(Real(1)/(10*pow(r2,Real(3)/2)*M_PI))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*pow(y,2)*(cospz-1)*sint);
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return -(Real(1)/(5*r))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*sinpz*sint);
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
			return 0;

		SETUPLOCALS

		return Real(1)/(10*r*M_PI)*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz-1)*cost);
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return  -(Real(1)/(10*pow(r2,Real(3)/2)))*(sinpr*pow(OMEGA,2)*M_PI*pow(x,2)*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz-1)*sint)
				-(Real(3)/(10*pow(r2,2)))*(cospr*OMEGA*pow(x,2)*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz-1)*sint)
				+(Real(1)/(10*r2))*(cospr*OMEGA*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz-1)*sint)
				+(Real(1)/(5*r2))*(cospr*OMEGA*x*(-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr*x/r-2*sinpr*x/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI*x/r2)*y*(cospz-1)*sint)
				+(Real(1)/(10*r*M_PI))*(sinpr*(-2*pow(M_PI,3)*pow(OMEGA,3)*cospr*pow(x,2)/r2-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr/r+6*sinpr*pow(x,2)/pow(r2,Real(5)/2)-6*cospr*OMEGA*M_PI*pow(x,2)/pow(r2,2)-2*sinpr/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI/r2)*y*(cospz-1)*sint)
				-(Real(1)/( 5*pow(r2,Real(3)/2)*M_PI))*(sinpr*(-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr*x/r-2*sinpr*x/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI*x/r2)*y*(cospz-1)*sint*x)
				+(Real(3)/(10*pow(r2,Real(5)/2)*M_PI))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz-1)*sint*pow(x,2))
				-(Real(1)/(10*pow(r2,Real(3)/2)*M_PI))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz-1)*sint);
	}

	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return -(Real(1)/(10*pow(r2,Real(3)/2)))*(sinpr*pow(OMEGA,2)*M_PI*pow(y,3)*(2*M_PI*OMEGA*cospr+2*sinpr/r)*(cospz-1)*sint)
			   -(Real(3)/(10*pow(r2,2)))*(cospr*OMEGA*pow(y,3)*(2*M_PI*OMEGA*cospr+2*sinpr/r)*(cospz-1)*sint)
			   +(Real(3)/(10*r2))*(cospr*OMEGA*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz-1)*sint)
			   +(Real(1)/( 5*r2))*(cospr*OMEGA*pow(y,2)*(-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr*y/r-2*sinpr*y/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI*y/r2)*(cospz-1)*sint)
			   +(Real(1)/(10*r*M_PI))*(sinpr*(-2*pow(M_PI,3)*pow(OMEGA,3)*cospr*pow(y,2)/r2-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr/r+6*sinpr*pow(y,2)/pow(r2,Real(5)/2)-6*cospr*OMEGA*M_PI*pow(y,2)/pow(r2,2)-2*sinpr/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI/r2)*y*(cospz-1)*sint)
			   +(Real(1)/( 5*r*M_PI))*(sinpr*(-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr*y/r-2*sinpr*y/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI*y/r2)*(cospz-1)*sint)
			   -(Real(1)/( 5*pow(r2,Real(3)/2)*M_PI))*(sinpr*(-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr*y/r-2*sinpr*y/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI*y/r2)*pow(y,2)*(cospz-1)*sint)
			   -(Real(3)/(10*pow(r2,Real(3)/2)*M_PI))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz-1)*sint)
			   +(Real(3)/(10*pow(r2,Real(5)/2)*M_PI))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*pow(y,3)*(cospz-1)*sint);
	}

	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return -(Real(2)/(5*r))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*cospz*M_PI*sint);
	}
};

#define RADIUS (2)

class articleTest1P : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		Real r = sqrt(x*x+y*y);
		Real sinpr = sin(2*M_PI/RADIUS*r);
		Real sinpz = sin(2*M_PI*z);
		Real sint = sin(t);

		return sinpr*sinpr*sinpz*sint;
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		Real r = sqrt(x*x+y*y);
		Real sinpr = sin(2*M_PI/RADIUS*r);
		Real cospr = cos(2*M_PI/RADIUS*r);
		Real sinpz = sin(2*M_PI*z);
		Real sint = sin(t);

		return 4*sinpr*sinpz*sint*cospr*x*M_PI/(RADIUS*r);
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		Real r = sqrt(x*x+y*y);
		Real sinpr = sin(2*M_PI/RADIUS*r);
		Real cospr = cos(2*M_PI/RADIUS*r);
		Real sinpz = sin(2*M_PI*z);
		Real sint = sin(t);

		return 4*sinpr*sinpz*sint*cospr*y*M_PI/(RADIUS*r);
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		Real r = sqrt(x*x+y*y);
		Real sinpr = sin(2*M_PI/RADIUS*r);
		Real cospz = cos(2*M_PI*z);
		Real sint = sin(t);

		return 2*M_PI*sinpr*sinpr*cospz*sint;
	}	

	Real difft(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return sinpr*sinpr*sinpz*cost;
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return  2*pow(cospr,2)*pow(OMEGA,2)*pow(M_PI,2)*pow(x,2)*sinpz*sint/r2
			   -2*pow(sinpr,2)*sinpz*sint*pow(OMEGA,2)*pow(M_PI,2)*pow(x,2)/r2
			   -2*sinpr*sinpz*sint*cospr*OMEGA*M_PI*pow(x,2)/pow(r2,Real(3)/2)
			   +2*sinpr*sinpz*sint*cospr*OMEGA*M_PI/r;
	}

	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return  2*pow(cospr,2)*pow(OMEGA,2)*pow(M_PI,2)*pow(y,2)*sinpz*sint/r2
			   -2*pow(sinpr,2)*sinpz*sint*pow(OMEGA,2)*pow(M_PI,2)*pow(y,2)/r2
			   -2*sinpr*sinpz*sint*cospr*OMEGA*M_PI*pow(y,2)/pow(r2,Real(3)/2)
			   +2*sinpr*sinpz*sint*cospr*OMEGA*M_PI/r;
	}
	
	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		return -4*sinpr*sinpr*sinpz*pow(M_PI,2)*sint;
	}
};

class articleTest1X2 : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
			return 0;

		sinpz = pow(z,2)*pow(1-z,2);

		return( (Real(2)/5)*pow(sinpr,2)*y*sinpz*sint*x/r2 );
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		sinpz = pow(z,2)*pow(1-z,2);

		return(  (Real(4)/5)*sinpr*y*sinpz*sint*pow(x,2)*cospr*OMEGA*M_PI/pow(r2,(Real(3)/2))
				-(Real(4)/5)*pow(sinpr,2)*y*sinpz*sint*pow(x,2)/pow(r2,2)
				+(Real(2)/5)*pow(sinpr,2)*y*sinpz*sint/r2 );
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		sinpz = pow(z,2)*pow(1-z,2);

		return  (Real(4)/5)*sinpr*pow(y,2)*sinpz*sint*x*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
			   +(Real(2)/5)*pow(sinpr,2)*sinpz*sint*x/r2
			   -(Real(4)/5)*pow(sinpr,2)*pow(y,2)*sinpz*sint*x/pow(r2,2);
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		cospz = 2*z*pow(1-z,2)-2*pow(z,2)*(1-z);

		return (Real(4)/5)*pow(sinpr,2)*y*cospz*M_PI*sint*x/r2;
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
			return 0;

		SETUPLOCALS
		
		sinpz = pow(z,2)*pow(1-z,2);

		return( (Real(2)/5)*pow(sinpr,2)*y*sinpz*cost*x/r2 );
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		sinpz = pow(z,2)*pow(1-z,2);

		return(  (Real( 4)/5)*pow(cospr,2)*pow(OMEGA,2)*pow(M_PI,2)*pow(x,3)*y*sinpz*sint/pow(r2,2)
				-(Real( 4)/1)*sinpr*y*sinpz*sint*pow(x,3)*cospr*OMEGA*M_PI/pow(r2,Real(5)/2)
				+(Real(12)/5)*sinpr*y*sinpz*sint*x*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
				-(Real( 4)/5)*pow(sinpr,2)*y*sinpz*sint*pow(x,3)*pow(OMEGA,2)*pow(M_PI,2)/pow(r2,2)
				+(Real(16)/5)*pow(sinpr,2)*y*sinpz*sint*pow(x,3)/pow(r2,3)
				-(Real(12)/5)*pow(sinpr,2)*y*sinpz*sint*x/pow(r2,2) );
	}

	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS

		sinpz = pow(z,2)*pow(1-z,2);

		return(  (Real( 4)/5)*pow(cospr,2)*pow(OMEGA,2)*pow(M_PI,2)*pow(y,3)*sinpz*sint*x/pow(r2,2)
				+(Real(12)/5)*sinpr*y*sinpz*sint*x*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
				-(Real( 4)/1)*sinpr*pow(y,3)*sinpz*sint*x*cospr*OMEGA*M_PI/pow(r2,Real(5)/2)
				-(Real( 4)/5)*pow(sinpr,2)*pow(y,3)*sinpz*sint*x*pow(OMEGA,2)*pow(M_PI,2)/pow(r2,2)
				-(Real(12)/5)*pow(sinpr,2)*y*sinpz*sint*x/pow(r2,2)
				+(Real(16)/5)*pow(sinpr,2)*pow(y,3)*sinpz*sint*x/pow(r2,3) );
	}

	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		sinpz = 2*pow(1-z,2)-8*z*(1-z)+2*pow(z,2);

		return( -(Real(8)/5)*pow(sinpr,2)*y*sinpz*sint*x/r2 );
	}
};

class articleTest1Y2 : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
			return 0;

		SETUPLOCALS
		
		sinpz = pow(z,2)*pow(1-z,2);

		return(  (Real(1)/5)*pow(sinpr,2)*pow(y,2)*sinpz*sint/r2
				-(Real(1)/5)*pow(sinpr,2)*pow(x,2)*sinpz*sint/r2 );
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		sinpz = pow(z,2)*pow(1-z,2);

		return  (Real(2)/5)*sinpr*pow(y,2)*sinpz*sint*cospr*OMEGA*M_PI*x/pow(r2,Real(3)/2)
			   -(Real(2)/5)*pow(sinpr,2)*pow(y,2)*sinpz*sint*x/pow(r2,2)
			   -(Real(2)/5)*sinpr*pow(x,3)*sinpz*sint*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
			   -(Real(2)/5)*pow(sinpr,2)*x*sinpz*sint/r2
			   +(Real(2)/5)*pow(sinpr,2)*pow(x,3)*sinpz*sint/pow(r2,2);
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		sinpz = pow(z,2)*pow(1-z,2);

		return  (Real(2)/5)*sinpr*pow(y,3)*sinpz*sint*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
			   +(Real(2)/5)*pow(sinpr,2)*y*sinpz*sint/r2
			   -(Real(2)/5)*pow(sinpr,2)*pow(y,3)*sinpz*sint/pow(r2,2)
			   -(Real(2)/5)*sinpr*y*sinpz*sint*pow(x,2)*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
			   +(Real(2)/5)*pow(sinpr,2)*y*sinpz*sint*pow(x,2)/pow(r2,2);

	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		cospz = 2*z*pow(1-z,2)-2*pow(z,2)*(1-z);

		return  (Real(2)/5)*pow(sinpr,2)*pow(y,2)*cospz*M_PI*sint/r2
			   -(Real(2)/5)*pow(sinpr,2)*pow(x,2)*cospz*M_PI*sint/r2;
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
			return 0;
		
		SETUPLOCALS
		
		sinpz = pow(z,2)*pow(1-z,2);

		return(  (Real(1)/5)*pow(sinpr,2)*pow(y,2)*sinpz*cost/r2
				-(Real(1)/5)*pow(sinpr,2)*pow(x,2)*sinpz*cost/r2 );
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		sinpz = pow(z,2)*pow(1-z,2);

		return  (Real(2)/5)*pow(cospr,2)*pow(OMEGA,2)*pow(M_PI,2)*pow(x,2)*pow(y,2)*sinpz*sint/pow(r2,2)
			   -(Real(2)/1)*sinpr*pow(y,2)*sinpz*sint*pow(x,2)*cospr*OMEGA*M_PI/pow(r2,Real(5)/2)
			   +(Real(2)/5)*sinpr*pow(y,2)*sinpz*sint*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
			   -(Real(2)/5)*pow(sinpr,2)*pow(y,2)*sinpz*sint*pow(x,2)*pow(OMEGA,2)*pow(M_PI,2)/pow(r2,2)
			   +(Real(8)/5)*pow(sinpr,2)*pow(y,2)*sinpz*sint*pow(x,2)/pow(r2,3)
			   -(Real(2)/5)*pow(sinpr,2)*pow(y,2)*sinpz*sint/pow(r2,2)
			   -(Real(2)/5)*pow(cospr,2)*pow(OMEGA,2)*pow(M_PI,2)*pow(x,4)*sinpz*sint/pow(r2,2)
			   -(Real(2)/1)*sinpr*pow(x,2)*sinpz*sint*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
			   +(Real(2)/1)*sinpr*pow(x,4)*sinpz*sint*cospr*OMEGA*M_PI/pow(r2,Real(5)/2)
			   +(Real(2)/5)*pow(sinpr,2)*pow(x,4)*sinpz*sint*pow(OMEGA,2)*pow(M_PI,2)/pow(r2,2)
			   +(Real(2)/1)*pow(sinpr,2)*sinpz*sint*pow(x,2)/pow(r2,2)
			   -(Real(2)/5)*pow(sinpr,2)*sinpz*sint/r2
			   -(Real(8)/5)*pow(sinpr,2)*pow(x,4)*sinpz*sint/pow(r2,3);
	}

	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		sinpz = pow(z,2)*pow(1-z,2);

		return  (Real(2)/5)*pow(cospr,2)*pow(OMEGA,2)*pow(M_PI,2)*pow(y,4)*sinpz*sint/pow(r2,2)
			   +(Real(2)/1)*sinpr*pow(y,2)*sinpz*sint*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
			   -(Real(2)/1)*sinpr*pow(y,4)*sinpz*sint*cospr*OMEGA*M_PI/pow(r2,Real(5)/2)
			   -(Real(2)/5)*pow(sinpr,2)*pow(y,4)*sinpz*sint*pow(OMEGA,2)*pow(M_PI,2)/pow(r2,2)
			   +(Real(2)/5)*pow(sinpr,2)*sinpz*sint/r2
			   -(Real(2)/1)*pow(sinpr,2)*pow(y,2)*sinpz*sint/pow(r2,2)
			   +(Real(8)/5)*pow(sinpr,2)*pow(y,4)*sinpz*sint/pow(r2,3)
			   -(Real(2)/5)*pow(cospr,2)*pow(OMEGA,2)*pow(M_PI,2)*pow(x,2)*pow(y,2)*sinpz*sint/pow(r2,2)
			   -(Real(2)/5)*sinpr*pow(x,2)*sinpz*sint*cospr*OMEGA*M_PI/pow(r2,Real(3)/2)
			   +(Real(2)/1)*sinpr*pow(y,2)*sinpz*sint*pow(x,2)*cospr*OMEGA*M_PI/pow(r2,Real(5)/2)
			   +(Real(2)/5)*pow(sinpr,2)*pow(y,2)*sinpz*sint*pow(x,2)*pow(OMEGA,2)*pow(M_PI,2)/pow(r2,2)
			   +(Real(2)/5)*pow(sinpr,2)*sinpz*sint*pow(x,2)/pow(r2,2)
			   -(Real(8)/5)*pow(sinpr,2)*pow(y,2)*sinpz*sint*pow(x,2)/pow(r2,3);
	}

	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		sinpz = 2*pow(1-z,2)-8*z*(1-z)+2*pow(z,2);

		return -(Real(4)/5)*pow(sinpr,2)*pow(y,2)*sinpz*pow(M_PI,2)*sint/r2
			   +(Real(4)/5)*pow(sinpr,2)*pow(x,2)*sinpz*pow(M_PI,2)*sint/r2;
	}
};

class articleTest1Z2 : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
			return 0;

		SETUPLOCALS
		
		cospz = -2*M_PI*(2*z*pow(1-z,2)-2*pow(z,2)*(1-z));

		return  (Real(1)/(10*r*M_PI))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz)*sint);
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		cospz = -2*M_PI*(2*z*pow(1-z,2)-2*pow(z,2)*(1-z));

		return   (Real(1)/(10*r2))*(cospr*OMEGA*x*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz)*sint)
			    +(Real(1)/(10*M_PI*r))*(sinpr*(-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr*x/r-2*sinpr*x/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI*x/r2)*y*(cospz)*sint)
				-(Real(1)/(10*pow(r2,Real(3)/2)*M_PI))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz)*sint*x);
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		cospz = -2*M_PI*(2*z*pow(1-z,2)-2*pow(z,2)*(1-z));

		return  (Real(1)/(10*r2))*(cospr*OMEGA*pow(y,2)*(2*M_PI*OMEGA*cospr+2*sinpr/r)*(cospz)*sint)
			   +(Real(1)/(10*r*M_PI))*(sinpr*(-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr*y/r-2*sinpr*y/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI*y/r2)*y*(cospz)*sint)
			   +(Real(1)/(10*r*M_PI))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*(cospz)*sint)
			   -(Real(1)/(10*pow(r2,Real(3)/2)*M_PI))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*pow(y,2)*(cospz)*sint);
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		sinpz = -2*M_PI*(2*pow(1-z,2)-8*z*(1-z)+2*pow(z,2));

		return -(Real(1)/(5*r))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*sinpz*sint);
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		if( abs(x) < 1.e-14 && abs(y) < 1.e-14 )
			return 0;

		SETUPLOCALS
		
		cospz = -2*M_PI*(2*z*pow(1-z,2)-2*pow(z,2)*(1-z));

		return Real(1)/(10*r*M_PI)*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz)*cost);
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		cospz = -2*M_PI*(2*z*pow(1-z,2)-2*pow(z,2)*(1-z));

		return  -(Real(1)/(10*pow(r2,Real(3)/2)))*(sinpr*pow(OMEGA,2)*M_PI*pow(x,2)*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz)*sint)
				-(Real(3)/(10*pow(r2,2)))*(cospr*OMEGA*pow(x,2)*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz)*sint)
				+(Real(1)/(10*r2))*(cospr*OMEGA*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz)*sint)
				+(Real(1)/(5*r2))*(cospr*OMEGA*x*(-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr*x/r-2*sinpr*x/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI*x/r2)*y*(cospz-1)*sint)
				+(Real(1)/(10*r*M_PI))*(sinpr*(-2*pow(M_PI,3)*pow(OMEGA,3)*cospr*pow(x,2)/r2-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr/r+6*sinpr*pow(x,2)/pow(r2,Real(5)/2)-6*cospr*OMEGA*M_PI*pow(x,2)/pow(r2,2)-2*sinpr/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI/r2)*y*(cospz)*sint)
				-(Real(1)/( 5*pow(r2,Real(3)/2)*M_PI))*(sinpr*(-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr*x/r-2*sinpr*x/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI*x/r2)*y*(cospz)*sint*x)
				+(Real(3)/(10*pow(r2,Real(5)/2)*M_PI))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz)*sint*pow(x,2))
				-(Real(1)/(10*pow(r2,Real(3)/2)*M_PI))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz)*sint);
	}

	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		cospz = -2*M_PI*(2*z*pow(1-z,2)-2*pow(z,2)*(1-z));

		return -(Real(1)/(10*pow(r2,Real(3)/2)))*(sinpr*pow(OMEGA,2)*M_PI*pow(y,3)*(2*M_PI*OMEGA*cospr+2*sinpr/r)*(cospz)*sint)
			   -(Real(3)/(10*pow(r2,2)))*(cospr*OMEGA*pow(y,3)*(2*M_PI*OMEGA*cospr+2*sinpr/r)*(cospz)*sint)
			   +(Real(3)/(10*r2))*(cospr*OMEGA*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz)*sint)
			   +(Real(1)/( 5*r2))*(cospr*OMEGA*pow(y,2)*(-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr*y/r-2*sinpr*y/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI*y/r2)*(cospz)*sint)
			   +(Real(1)/(10*r*M_PI))*(sinpr*(-2*pow(M_PI,3)*pow(OMEGA,3)*cospr*pow(y,2)/r2-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr/r+6*sinpr*pow(y,2)/pow(r2,Real(5)/2)-6*cospr*OMEGA*M_PI*pow(y,2)/pow(r2,2)-2*sinpr/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI/r2)*y*(cospz)*sint)
			   +(Real(1)/( 5*r*M_PI))*(sinpr*(-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr*y/r-2*sinpr*y/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI*y/r2)*(cospz)*sint)
			   -(Real(1)/( 5*pow(r2,Real(3)/2)*M_PI))*(sinpr*(-2*pow(M_PI,2)*pow(OMEGA,2)*sinpr*y/r-2*sinpr*y/pow(r2,Real(3)/2)+2*cospr*OMEGA*M_PI*y/r2)*pow(y,2)*(cospz)*sint)
			   -(Real(3)/(10*pow(r2,Real(3)/2)*M_PI))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*(cospz)*sint)
			   +(Real(3)/(10*pow(r2,Real(5)/2)*M_PI))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*pow(y,3)*(cospz)*sint);
	}

	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		SETUPLOCALS
		
		cospz = -2*M_PI*(-12+24*z);

		return -(Real(2)/(5*r))*(sinpr*(2*M_PI*OMEGA*cospr+2*sinpr/r)*y*cospz*M_PI*sint);
	}
};

class articleTest1 : public basics::function3D {
public:
	Real val(Real x, Real y, Real z, Real t) const
	{
		return -cos((Real(7)/8)*M_PI*(x*x+y*y))*(cos(M_PI*z)-cos(3*M_PI*z))*exp(z);
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		return (Real(7)/4)*sin((Real(7)/8)*M_PI*(x*x+y*y))*M_PI*x*(cos(M_PI*z)-cos(3*M_PI*z))*exp(z);
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		return (Real(7)/4)*sin((Real(7)/8)*M_PI*(x*x+y*y))*M_PI*y*(cos(M_PI*z)-cos(3*M_PI*z))*exp(z);
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		return -cos((Real(7)/8)*M_PI*(x*x+y*y))*(-sin(M_PI*z)*M_PI+3*sin(3*M_PI*z)*M_PI)*exp(z)-
			    cos((Real(7)/8)*M_PI*(x*x+y*y))*(cos(M_PI*z)-cos(3*M_PI*z))*exp(z);
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		return( 0 );
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		return (Real(49)/16)*cos((Real(7)/8)*M_PI*(x*x+y*y))*M_PI*M_PI*x*x*(cos(M_PI*z)-cos(3*M_PI*z))*exp(z)+(Real(7)/4)*sin((Real(7)/8)*M_PI*(x*x+y*y))*M_PI*(cos(M_PI*z)-cos(3*M_PI*z))*exp(z);
	}

	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		return (Real(49)/16)*cos((Real(7)/8)*M_PI*(x*x+y*y))*M_PI*M_PI*y*y*(cos(M_PI*z)-cos(3*M_PI*z))*exp(z)+(Real(7)/4)*sin((Real(7)/8)*M_PI*(x*x+y*y))*M_PI*(cos(M_PI*z)-cos(3*M_PI*z))*exp(z);
	}

	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		return -cos((Real(7)/8)*M_PI*(x*x+y*y))*(-cos(M_PI*z)*M_PI*M_PI+9*cos(3*M_PI*z)*M_PI*M_PI)*exp(z)-2*cos((Real(7)/8)*M_PI*(x*x+y*y))*(-sin(M_PI*z)*M_PI+3*sin(3*M_PI*z)*M_PI)*exp(z)-cos((Real(7)/8)*M_PI*(x*x+y*y))*(cos(M_PI*z)-cos(3*M_PI*z))*exp(z);
	}
};

class boxTest3DX : public basics::function3D {
	public:
		Real val(Real x, Real y, Real z, Real t) const
		{
			return( M_PI/5*sin(M_PI*x)*sin(M_PI*x)*sin(Real(2)*M_PI*y)*sin(2*M_PI*z)*sin(t) );
		}

		Real diffx(Real x, Real y, Real z, Real t) const
		{
			return( 2*pow(M_PI,2)/5*sin(M_PI*x)*cos(M_PI*x)*sin(Real(2)*M_PI*y)*sin(2*M_PI*z)*sin(t) );
		}

		Real diffy(Real x, Real y, Real z, Real t) const
		{
			return( 2*pow(M_PI,2)/5*sin(M_PI*x)*sin(M_PI*x)*cos(Real(2)*M_PI*y)*sin(2*M_PI*z)*sin(t) );
		}

		Real diffz(Real x, Real y, Real z, Real t) const
		{
			return( 2*pow(M_PI,2)/5*sin(M_PI*x)*sin(M_PI*x)*sin(Real(2)*M_PI*y)*cos(2*M_PI*z)*sin(t) );
		}

		Real difft(Real x, Real y, Real z, Real t) const
		{
			return( M_PI/5*sin(M_PI*x)*sin(M_PI*x)*sin(Real(2)*M_PI*y)*sin(2*M_PI*z)*cos(t) );
		}

		Real diff2x(Real x, Real y, Real z, Real t) const
		{
			return( 2*pow(M_PI,3)/5*sin(t)*sin(Real(2)*M_PI*y)*sin(2*M_PI*z)*(cos(M_PI*x)*cos(M_PI*x)-sin(M_PI*x)*sin(M_PI*x)) );
		}

		Real diff2y(Real x, Real y, Real z, Real t) const
		{
			return( -4*pow(M_PI,3)/5*sin(M_PI*x)*sin(M_PI*x)*sin(Real(2)*M_PI*y)*sin(2*M_PI*z)*sin(t) );
		}

		Real diff2z(Real x, Real y, Real z, Real t) const
		{
			return( -4*pow(M_PI,3)/5*sin(M_PI*x)*sin(M_PI*x)*sin(Real(2)*M_PI*y)*sin(2*M_PI*z)*sin(t) );
		}

		static int id() 
		{
			return( 1 );
		}
};

class boxTest3DY : public basics::function3D {
	
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( -M_PI/10*sin(2*M_PI*x)*sin(M_PI*y)*sin(M_PI*y)*sin(2*M_PI*z)*sin(t) );
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		return( -pow(M_PI,2)/5*cos(2*M_PI*x)*sin(M_PI*y)*sin(M_PI*y)*sin(2*M_PI*z)*sin(t) );
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		return( -pow(M_PI,2)/5*sin(2*M_PI*x)*sin(M_PI*y)*cos(M_PI*y)*sin(2*M_PI*z)*sin(t) );
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		return( -pow(M_PI,2)/5*sin(2*M_PI*x)*sin(M_PI*y)*sin(M_PI*y)*cos(2*M_PI*z)*sin(t) );
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		return( -M_PI/10*sin(2*M_PI*x)*sin(M_PI*y)*sin(M_PI*y)*sin(2*M_PI*z)*cos(t) );
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		return( 2*pow(M_PI,3)/5*sin(2*M_PI*x)*sin(M_PI*y)*sin(M_PI*y)*sin(2*M_PI*z)*sin(t) );
	}
	
	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		return( pow(M_PI,3)/5*sin(2*M_PI*x)*sin(t)*sin(2*M_PI*z)*(sin(M_PI*y)*sin(M_PI*y)-cos(M_PI*y)*cos(M_PI*y)) );
	}
	
	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		return( 2*pow(M_PI,3)/5*sin(2*M_PI*x)*sin(M_PI*y)*sin(M_PI*y)*sin(2*M_PI*z)*sin(t) );
	}

	static int id() 
	{
		return( 1 );
	}
};


class boxTest3DZ : public basics::function3D {
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( -M_PI/10*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(M_PI*z)*sin(M_PI*z)*sin(t));
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		return( -pow(M_PI,2)/5*cos(2*M_PI*x)*sin(2*M_PI*y)*sin(M_PI*z)*sin(M_PI*z)*sin(t));
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		return( -pow(M_PI,2)/5*sin(2*M_PI*x)*cos(2*M_PI*y)*sin(M_PI*z)*sin(M_PI*z)*sin(t));
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	 {
		return( -pow(M_PI,2)/5*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(M_PI*z)*cos(M_PI*z)*sin(t) );
	}

	Real difft(Real x, Real y, Real z, Real t) const
	{
		return( -M_PI/10*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(M_PI*z)*sin(M_PI*z)*cos(t));
	}

	Real diff2x(Real x, Real y, Real z, Real t) const
	{
		return( 2*pow(M_PI,3)/5*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(M_PI*z)*sin(M_PI*z)*sin(t));
	}
	
	Real diff2y(Real x, Real y, Real z, Real t) const
	{
		return( 2*pow(M_PI,3)/5*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(M_PI*z)*sin(M_PI*z)*sin(t));
	}
	
	Real diff2z(Real x, Real y, Real z, Real t) const
	{
		return( pow(M_PI,3)/5*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(t)*(sin(M_PI*z)*sin(M_PI*z)-cos(M_PI*z)*cos(M_PI*z)) );
	}

	static int id() 
	{
		return( 1 );
	}
};

class boxTest3DP : public basics::function3D {
	
	Real val(Real x, Real y, Real z, Real t) const
	{
		return( cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*sin(t) );
	}

	Real diffx(Real x, Real y, Real z, Real t) const
	{
		return( -M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*sin(t) );
	}

	Real diffy(Real x, Real y, Real z, Real t) const
	{
		return( M_PI*cos(M_PI*x)*cos(M_PI*y)*sin(M_PI*z)*sin(t) );
	}

	Real diffz(Real x, Real y, Real z, Real t) const
	{
		return( M_PI*cos(M_PI*x)*sin(M_PI*y)*cos(M_PI*z)*sin(t) );
	}

	static int id() 
	{
		return( 1 );
	}
};

#endif

