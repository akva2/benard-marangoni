#include "uzawa.h"
#include "sem.h"

using namespace mnl;
using namespace std;

SEMCahouetChabardPreconditioner::SEMCahouetChabardPreconditioner(basics::geometryStack& geometry,
    const basics::Matrix& D,
    const basics::Matrix& GLL2G,
    const basics::Vector& weightGL,
    const basics::Vector& weight,
    const basics::Vector& grid,
    const basics::Vector& gridGL) :
  m_eval(geometry,D,GLL2G,weightGL,weight),
  m_invMass(geometry,weightGL,GLL2G),
  m_pre(geometry,weight,weightGL,grid,gridGL,GLL2G,D,0)
{
}

void SEMCahouetChabardPreconditioner::evaluate(basics::matrixStack& res, const basics::matrixStack& p) const
{
  basics::matrixStack* buffer = utilities::g_manager.clone(p);

  *buffer = p;
  utilities::CGSolver::solve(*buffer,m_eval,buffer->getDotter(),1.e-8);
  *buffer *= Real(3)/(2*m_Dt);
  m_invMass.evaluate(res,p);
  res += *buffer;

  utilities::g_manager.unlock(buffer);
}

SEMUzawaEvaluator::SEMUzawaEvaluator(basics::geometryStack& geometry,
    const basics::Matrix& D,
    const basics::Matrix& GLL2G,
    const basics::Vector& grid_gll,
    const basics::Vector& weight_gll,
    const basics::Vector& weight_gl,
    SEMLaplacianEvaluator& laplacian,
    bool  preconditioned,
    int rank, int size) :
  m_divergence(geometry,D,GLL2G,weight_gl,NULL,rank,size),
  m_laplacian(laplacian),
  m_gradient(geometry,D,GLL2G,weight_gl,NULL,rank,size),
  m_geometry(geometry), m_pre(NULL), m_rank(rank), m_size(size)
{
  m_G = geometry.getInterpolatedGeometryDerivatives(GLL2G);

  m_division = m_geometry.getDivisionInfo(size);
  m_W = utilities::g_manager.aquireMatrixStackField("working field",D.rows(),D.cols(),m_division[m_rank].elements.size());

  if( preconditioned ) {
    m_pre = new SEMFEMLaplacianPreconditioner(geometry,grid_gll,weight_gll,0,rank,size);
    m_pre->m_nosum = m_laplacian.m_nosum = utilities::g_manager.aquireMatrixStack("unsummed residual",weight_gl.length()+2,weight_gl.length()+2,m_division[m_rank].elements.size());
  } else {
    m_laplacian.m_nosum = NULL;
  }
}

SEMUzawaEvaluator::~SEMUzawaEvaluator()
{
  if( m_pre ) {
    utilities::g_manager.unlock(m_pre->m_nosum);
    delete m_pre;
  }
  for( int i=0;i<m_G.size();++i )
    delete m_G[i];
  utilities::g_manager.unlock(m_W);
}

void SEMUzawaEvaluator::evaluate(basics::matrixStack& res, const basics::matrixStack& p) const
{
  m_gradient.evaluate(*m_W,p);
  m_geometry.maskField(*m_W);
  MPIGeometryDotter<basics::geometryStack> dotter(m_geometry,m_division[m_rank].elements,m_rank,m_size);
  if( m_pre ) {
    *m_laplacian.m_nosum = m_W->X();
    m_geometry.dssum(m_W->X(),m_rank,m_size);
    utilities::CGSolver::solve(m_W->X(),m_laplacian,const_cast<SEMFEMLaplacianPreconditioner&>(*m_pre),dotter,1.e-10);
    *m_laplacian.m_nosum = m_W->Y();
    m_geometry.dssum(m_W->Y(),m_rank,m_size);
    utilities::CGSolver::solve(m_W->Y(),m_laplacian,const_cast<SEMFEMLaplacianPreconditioner&>(*m_pre),dotter,1.e-10);
  } else {
    m_geometry.dssum(*m_W,m_rank,m_size);
    utilities::CGSolver::solve(m_W->X(),m_laplacian,dotter,1.e-10);
    utilities::CGSolver::solve(m_W->Y(),m_laplacian,dotter,1.e-10);
  }
  m_divergence.evaluate(res,*m_W);
}

extrudedCahouetChabardPreconditioner::
extrudedCahouetChabardPreconditioner(const extrudedInverseMassEvaluator& invMass,
    extrudedConsistentPressureEvaluator& consistent,
    Real nu) :
  m_invMass(invMass), m_consistent(consistent), m_nu(nu)
{
}

void extrudedCahouetChabardPreconditioner::evaluate(basics::matricesStack& res,
    const basics::matricesStack& p) const
{
  basics::matricesStack* buffer = utilities::g_manager.clone(p);

  *buffer = p;
  m_consistent.solve(*buffer);
  *buffer *= m_nu;
  m_invMass.evaluate(res,p);
  res += *buffer;

  utilities::g_manager.unlock(buffer);
}

extrudedUzawaEvaluator::extrudedUzawaEvaluator(basics::geometryStack& geometry,
    const basics::Matrix& GLL2G,
    const basics::Vector& weight_gl,
    const legendreLegendreW::poissonSolver& SP,
    extrudedLaplacianEvaluator& laplacian,
    int rank, int size,
    extrudedLaplacianEvaluator* laplacian2) :
  m_uzawa2D   (geometry,SP.Dx(),GLL2G,SP.gridX(),SP.weightX(),weight_gl,laplacian.m_eval,false,rank,size),
  m_gradient  (m_uzawa2D.m_gradient),
  m_laplacian (laplacian),
  m_laplacian2(laplacian2),
  m_divergence(m_uzawa2D.m_divergence),
  m_geometry  (geometry),
  m_rank		(rank),
  m_size		(size)
{
}

void extrudedUzawaEvaluator::evaluate(basics::matricesStack& res, const basics::matricesStack& p) const
{
  basics::Field3<basics::matricesStack>* W = 
    utilities::g_manager.aquireMatricesStackField("temp field X",p[0].rows()+2,
        p[0].cols()+2,p[0].matrices()+2,
        p.size());
  basics::matricesStack* Wt = 
    utilities::g_manager.aquireMatricesStack("temp field temp",p[0].rows()+2,
        p[0].cols()+2,p[0].matrices()+2,
        p.size());

  m_gradient.evaluate(*W,p);
  *Wt = W->X();
  mask(*W,m_geometry,m_rank,m_size,m_laplacian.m_bc);
  dssum(W->X(),m_geometry,m_rank,m_size,m_laplacian.m_bc);
  const_cast<extrudedLaplacianEvaluator&>(m_laplacian).solve(W->X(),*Wt);
  *Wt = W->Y();
  dssum(W->Y(),m_geometry,m_rank,m_size,m_laplacian.m_bc);
  const_cast<extrudedLaplacianEvaluator&>(m_laplacian).solve(W->Y(),*Wt);
  *Wt = W->Z();
  dssum(W->Z(),m_geometry,m_rank,m_size,m_laplacian.m_bc);
  if( m_laplacian2 )
    const_cast<extrudedLaplacianEvaluator&>(*m_laplacian2).solve(W->Z(),*Wt);
  else
    const_cast<extrudedLaplacianEvaluator&>(m_laplacian).solve(W->Z(),*Wt);
  mask(*W,m_geometry,m_rank,m_size,m_laplacian.m_bc);
  m_divergence.evaluate(res,*W);

  res -= res.sum()/res.length();

  utilities::g_manager.unlock(W);
  utilities::g_manager.unlock(Wt);
}

void extrudedUzawaEvaluator::source(basics::Field3<basics::matricesStack>& u, 
    const basics::Vector&     grid, const basics::function3D& ux, 
    const basics::function3D& uy,   const basics::function3D& uz, 
    const basics::function3D& p,    Real t, int elem)
{
  for( int i=0;i<u.X().size();++i )
    source(u,grid,ux,uy,uz,p,t,i,elem+i);
}

void extrudedUzawaEvaluator::source(basics::Field3<basics::matricesStack>& u, 
    const basics::Vector&     grid, const basics::function3D& ux, 
    const basics::function3D& uy,   const basics::function3D& uz, 
    const basics::function3D& p,    Real t, int i, int elem)
{
  for( int alpha=0;alpha<grid.length();++alpha ) {
    for( int beta=0;beta<grid.length();++beta ) {
      pair<Real,Real> point = m_geometry[elem].getGH().evaluate(grid[beta],
          grid[alpha],grid);
      for( int gamma=0;gamma<grid.length();++gamma ) {
        if( abs(point.first) < std::numeric_limits<Real>::epsilon() && 
            abs(point.second) < std::numeric_limits<Real>::epsilon() ) {
          u.X()[i][gamma][alpha][beta] = 0;
          u.Y()[i][gamma][alpha][beta] = 0;
          u.Z()[i][gamma][alpha][beta] = 0;
        } else {
          u.X()[i][gamma][alpha][beta] = -(ux.diff2x(point.first,point.second,grid[gamma],t)+
              ux.diff2y(point.first,point.second,grid[gamma],t)+
              ux.diff2z(point.first,point.second,grid[gamma],t))
            +p.diffx(point.first,point.second,grid[gamma],t);
          u.Y()[i][gamma][alpha][beta] = -(uy.diff2x(point.first,point.second,grid[gamma],t)+
              uy.diff2y(point.first,point.second,grid[gamma],t)+
              uy.diff2z(point.first,point.second,grid[gamma],t))
            +p.diffy(point.first,point.second,grid[gamma],t);
          u.Z()[i][gamma][alpha][beta] = -(uz.diff2x(point.first,point.second,grid[gamma],t)+
              uz.diff2y(point.first,point.second,grid[gamma],t)+
              uz.diff2z(point.first,point.second,grid[gamma],t))
            +p.diffz(point.first,point.second,grid[gamma],t);
        }
      }
    }
  }
}

