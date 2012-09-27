#include "laplacian.h"
#include "poissonsolver-fe.h"
#include "mnl/hdf5.h"
#include "sem.h"

#ifdef HAS_MPI
#include "mpi.h"
#endif

using namespace linearElements;
using namespace mnl;
using namespace std;

deformedLaplacianEvaluator::
deformedLaplacianEvaluator(const basics::matrixStack& G,
    const basics::Matrix& D,
    const basics::Vector& weight,
    basics::Matrix& Uxi,
    basics::Matrix& Ueta,
    basics::Matrix& t1,
    basics::Matrix& t2,
    basics::Matrix& t3) :
  m_G(G), m_D(D), m_weight(weight),
  m_Uxi(Uxi), m_Ueta(Ueta), m_t1(t1), m_t2(t2), m_t3(t3)
{
}

void deformedLaplacianEvaluator::evaluate(basics::Matrix& res,
    const basics::Matrix& u, bool bMask) const
{
  basics::multTranspose(m_Uxi,m_D,u,'N','N');
  basics::multTranspose(m_Ueta,u,m_D,'N','T');

  basics::multPointwise(m_t1,m_G[0],m_Uxi);
  basics::multPointwise(m_t3,m_G[1],m_Ueta);
  m_t1 += m_t3;
  basics::multPointwise(m_t2,m_G[1],m_Uxi);
  basics::multPointwise(m_t3,m_G[2],m_Ueta);
  m_t2 += m_t3;

  massReference(m_t1,m_weight);
  massReference(m_t2,m_weight);

  basics::multTranspose(res,m_D,m_t1,'T','N');
  basics::multTranspose(res,m_t2,m_D,'N','N',mnlRealOne,mnlRealOne);
  if( bMask ) {
    res[0].clear();
    res[u.cols()-1].clear();
    res.clearRow(0);
    res.clearRow(u.rows()-1);
  }
}

deformed3DLaplacianEvaluator::deformed3DLaplacianEvaluator(const basics::matricesStack& G,
    const basics::Matrix& D,
    const basics::Vector& weight,
    basics::Matrices& Uxi,
    basics::Matrices& Ueta,
    basics::Matrices& Ugamma,
    basics::Matrices& t1,
    basics::Matrices& t2,
    basics::Matrices& t3) :
  m_G(G), m_D(D), m_weight(weight),
  m_Uxi(Uxi), m_Ueta(Ueta), m_Ugamma(Ugamma),
  m_t1(t1), m_t2(t2), m_t3(t3)
{
}

void deformed3DLaplacianEvaluator::evaluate(basics::Matrices& res,
    const basics::Matrices& u, bool bMask) const
{
  for( int l=0;l<u.matrices();++l ) {
    basics::multTranspose(m_Uxi[l],m_D,u[l],'N','N');
    basics::multTranspose(m_Ueta[l],u[l],m_D,'N','T');
  }
  basics::applyLocalGlobal(m_Ugamma,u,m_D,'N','T');

  basics::multPointwise(m_t1,m_Uxi,m_G[0]);
  basics::multPointwise(m_t2,m_Ueta,m_G[1]);
  basics::multPointwise(m_t3,m_Ugamma,m_G[2]);
  m_t1 += m_t2;
  m_t1 += m_t3;
  massReference(m_t1,m_weight);
  for( int l=0;l<u.matrices();++l )
    basics::multTranspose(res[l],m_D,m_t1[l],'T','N');

  basics::multPointwise(m_t1,m_Uxi,m_G[1]);
  basics::multPointwise(m_t2,m_Ueta,m_G[3]);
  basics::multPointwise(m_t3,m_Ugamma,m_G[4]);
  m_t1 += m_t2;
  m_t1 += m_t3;
  massReference(m_t1,m_weight);
  for( int l=0;l<u.matrices();++l )
    basics::multTranspose(res[l],m_t1[l],m_D,'N','N',mnlRealOne,mnlRealOne);

  basics::multPointwise(m_t1,m_Uxi,m_G[2]);
  basics::multPointwise(m_t2,m_Ueta,m_G[4]);
  basics::multPointwise(m_t3,m_Ugamma,m_G[5]);
  m_t1 += m_t2;
  m_t1 += m_t3;
  massReference(m_t1,m_weight);

  basics::applyLocalGlobal(res,m_t1,m_D,'N','N',0,mnlRealOne,mnlRealOne);
}

deformedTensoredLaplacianPreconditioner::deformedTensoredLaplacianPreconditioner(const basics::spectralElement2D& G,
    const legendreLegendreW::poissonSolver& SP,
    const Real nu) :
  m_SP(SP), m_nu(nu)
{
  pair<Real,Real> L = G.getAdjustedSize(SP.weightX(),SP.gridX());
  m_LyLx = L.second/L.first;
}

void deformedTensoredLaplacianPreconditioner::evaluate(basics::Matrix& res, const basics::Matrix& u) const
{
  legendreLegendreW::copyPreserveBorders(res,u);
  m_SP.solve(res,m_nu,1,m_LyLx);
}

SEMLaplacianEvaluator::SEMLaplacianEvaluator(basics::geometryStack& geometry,
    const basics::Matrix& D,
    const basics::Vector& weight,
    const std::vector<mnl::basics::matrixStack*>* G,
    int rank, int size, BC bc) :
  m_geometry(geometry),
  m_D(D), m_weight(weight), m_nu(0), m_mine(true), m_nosum(NULL), m_rank(rank), m_size(size), m_tag(0),
  m_bc(bc)
{
  m_division = geometry.getDivisionInfo(size);
  if( G ) {
    m_G = *G;
    m_mine = false;
  } else {
    for( int i=0;i<3;++i )
      m_G.push_back(new basics::matrixStack("geometry factors",D.rows(),
            D.cols(),geometry.size()));
    for( int i=0;i<m_geometry.size();++i) {
      const basics::matrixStack& GD = geometry[i].getGH().getGeometryDerivatives();
      const basics::Matrix& J = geometry[i].getGH().getJacobian();
      for( int j=0;j<D.cols();++j ) {
        for( int k=0;k<D.rows();++k ) {
          Real Ja = J[j][k];
          (*m_G[0])[i][j][k] = Real( 1)/Ja*(GD[1][j][k]*GD[1][j][k] + GD[3][j][k]*GD[3][j][k]);
          (*m_G[1])[i][j][k] = Real(-1)/Ja*(GD[0][j][k]*GD[1][j][k]  + GD[2][j][k]*GD[3][j][k]);
          (*m_G[2])[i][j][k] = Real( 1)/Ja*(GD[0][j][k]*GD[0][j][k]   + GD[2][j][k]*GD[2][j][k]);
        }
      }
    }
  }
  m_buffer = utilities::g_manager.aquireMatrixStack("laplacian buffer",D.rows(),D.cols(),m_division[m_rank].elements.size());
  m_Uxi  = utilities::g_manager.clone(*m_buffer);
  m_Ueta = utilities::g_manager.clone(*m_Uxi);
  m_t1   = utilities::g_manager.clone(*m_Uxi);
  m_t2   = utilities::g_manager.clone(*m_Uxi);
  m_t3   = utilities::g_manager.clone(*m_Uxi);
  m_view = new basics::componentView<basics::Matrix,basics::matrixStack>(m_G);
  for( int i=0;i<m_geometry.size();++i ) {
    m_evals.push_back(new deformedLaplacianEvaluator((*m_view)[i],
          m_D,m_weight,(*m_Uxi)[i],
          (*m_Ueta)[i],(*m_t1)[i],
          (*m_t2)[i],(*m_t3)[i]));
  }
}

SEMLaplacianEvaluator::~SEMLaplacianEvaluator()
{
  delete m_view;
  if( m_mine ) {
    for( unsigned int i=0;i<m_G.size();++i )
      delete m_G[i];
  }
  for( int i=0;i<m_evals.size();++i )
    delete m_evals[i];

  utilities::g_manager.unlock(m_Uxi);
  utilities::g_manager.unlock(m_Ueta);
  utilities::g_manager.unlock(m_t1);
  utilities::g_manager.unlock(m_t2);
  utilities::g_manager.unlock(m_t3);
  utilities::g_manager.unlock(m_buffer);
}

void SEMLaplacianEvaluator::evaluate(basics::matrixStack& res,
    const basics::matrixStack& u, bool bMask, bool doSum) const
{
  //    utilities::g_profiler.start(3);
  int max=m_division[m_rank].elements.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    int elem = m_division[m_rank].elements[i];
    m_evals[elem]->evaluate(res[i],u[i],false);
    if( m_nu > 1.e-14 ) {
      basics::multPointwise((*m_buffer)[i],u[i],(*m_geometry.m_mass)[elem]);
      res[elem].axpy(m_nu,(*m_buffer)[i]);
    }
  }
  if( m_bc == HOMOGENOUS_DIRICHLET && bMask )
    m_geometry.mask(res);
  if( doSum && m_nosum )
    *m_nosum = res;
  if( doSum ) {
    if( m_bc == PERIODIC )
      m_geometry.periodicDssum(res);
    else 
      m_geometry.dssum(res);
  }
}

SEMFEMLaplacianPreconditioner::SEMFEMLaplacianPreconditioner(basics::geometryStack& geometry,
    const basics::Vector& grid,
    const basics::Vector& weight,
    const Real nu,
    int rank, int size,
    SEMLaplacianEvaluator::BC bc) :
  m_geometry(geometry), m_nu(nu), m_A(NULL), m_nosum(NULL),
  m_rank(rank), m_size(size), m_tag(0), m_bc(bc)
{
  m_group = geometry.getCoarseGroups(size);
  basics::coarseGrid group = geometry.getCoarseGroups(1)[0];
  /* based on the first FULL OVERLAP part */
  int n=0;
  while( n < group.size() && group[n].type != basics::coarseDescriptor::FULL_OVERLAP )
    ++n;

  basics::Vector Lx("Lx",group[n].size1);
  basics::Vector Ly("Ly",group[n].size2);

  for( int j=0;j<group[n].size1;++j )
    Lx[j] = m_geometry[group[n].elements[j]].getAdjustedSize(weight,grid).first;
  for( int j=0;j<group[n].size2;++j )
    Ly[j] = m_geometry[group[n].elements[group[n].size1-1+
      j*group[n].size1]].getAdjustedSize(weight,grid).second;

  m_SP = new poissonSolver(poissonSolver::Homogenous,grid,Lx,Ly);

  int n2 = 0;
  while( n2 < group.size() && group[n2].type == basics::coarseDescriptor::FULL_OVERLAP )
    ++n2;

  basics::Vector Lx2("Lx2",group[n2].size1);
  basics::Vector Ly2("Ly2",group[n2].size2);

  for( int j=0;j<group[n2].size1;++j )
    Lx2[j] = m_geometry[group[n2].elements[j]].getAdjustedSize(weight,grid).first;
  for( int j=0;j<group[n2].size2;++j )
    Ly2[j] = m_geometry[group[n2].elements[group[n2].size1-1+
      j*group[n2].size1]].getAdjustedSize(weight,grid).second;

  m_SP2 = new poissonSolver(m_bc==SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN?
      poissonSolver::Nonhomogenous : poissonSolver::HomogenousLeft,
      grid,Lx2,Ly2);

  for( int c=0;c<group.size();++c ) {
    Real lx=0;
    Real ly=0;
    Real area=0;

    if( group[c].type2 == basics::coarseDescriptor::FLIP_XY ) {
      for( int j=0;j<group[c].size1;++j )
        lx += m_geometry[group[c].elements[j]].getSize(weight,grid).second;
      for( int j=0;j<group[c].size2;++j )
        ly += m_geometry[group[c].elements[group[c].size1-1+
          j*group[c].size1]].getSize(weight,grid).first;
    } else {
      for( int j=0;j<group[c].size1;++j )
        lx += m_geometry[group[c].elements[j]].getSize(weight,grid).first;
      for( int j=0;j<group[c].size2;++j )
        ly += m_geometry[group[c].elements[group[c].size1-1+
          j*group[c].size1]].getSize(weight,grid).second;
    }
    for( int i=0;i<group[c].elements.size();++i )
      area += m_geometry[group[c].elements[i]].getVolume(weight);

    Real scale = sqrt(area/(lx*ly));
    lx *= scale;
    ly *= scale;

    m_Lx.push_back(lx);
    m_Ly.push_back(ly);
  }
  m_coarse  = m_geometry.getCoarseBuffer(grid.length(),grid.length(),
      m_group[m_rank],"coarse",false,
      m_bc==SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN);
  m_coarse2 = m_geometry.getCoarseBuffer(grid.length(),grid.length(),
      m_group[m_rank],"coarse2",false,
      m_bc==SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN);

  if( m_geometry.hasCoarseSolver() ) {
    basics::coarseDescriptor desc = m_geometry.getRestrictionGridInfo();
    if( desc.size1 > 0 ) {
      m_work = utilities::g_manager.aquireMatrixStack("working buffer",desc.size1,desc.size1,desc.size3);
      if( bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
        m_A    = new basics::Matrix("helmholtz",desc.size2,desc.size2);
      if( bc == SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN )
        m_A    = new basics::Matrix("helmholtz",desc.size4,desc.size4);

      basics::Matrix* B    = utilities::g_manager.clone(*m_A);
      m_LG = utilities::g_manager.aquireMatrixStack("local to global",desc.size1,desc.size1,desc.size3);

      basics::Vector coarseGrid("coarse",desc.size1);
      std::string hdf("AB5.hdf5");
      if( desc.size1 == 3 )
        hdf = "bigcircle.hdf5";

      HDF5::HDF5Reader reader(hdf);
      utilities::GLL::GaussLobattoLegendreGrid(coarseGrid);
      coarseGrid += 1;
      coarseGrid /= 2;
      if( bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET ) {
        reader.read(*m_A,"A");
        reader.read(*B,"B");
        reader.read(*m_LG,"LG");
      }
      if( bc == SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN ) {
        reader.read(*m_A,"An");
        reader.read(*B,"Bn");
        reader.read(*m_LG,"LGn");
      }
      m_A->axpy(m_nu,*B);
      m_restricted = utilities::g_manager.aquireVector("restricted",m_A->rows());
      if( bc == SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN && m_nu < 1.e-12 ) {
        basics::Matrix* A = m_A;
        m_A = new basics::Matrix(m_A->submatrix(utilities::Range::colon(0,desc.size4-2),utilities::Range::colon(0,desc.size4-2)));
        delete A;
        m_restricted1 = new basics::Vector("restricted minus one",m_restricted->length()-1,false,m_restricted->data());
      } else
        m_restricted1 = m_restricted;

      *m_A = m_A->choleskyFactorize();

      utilities::g_manager.unlock(B);

      setupRestrictionOperators(m_RT,geometry,grid,group[n].size1,group[n].size2,
          group[n2].size1,group[n2].size2,bc);
      m_restrict  = utilities::g_manager.aquireMatrix("restrict buffer2",(*m_coarse)[0].rows(),desc.size1);
    }
  } else
    m_restricted = NULL;
}

SEMFEMLaplacianPreconditioner::~SEMFEMLaplacianPreconditioner()
{
  delete m_SP;
  delete m_SP2;
  if( m_A )
    delete m_A;

  delete m_coarse;
  delete m_coarse2;

  if( m_restricted1 != m_restricted )
    delete m_restricted1;

  utilities::g_manager.unlock(m_LG);
  utilities::g_manager.unlock(m_restricted);
  utilities::g_manager.unlock(m_restrict);
  utilities::g_manager.unlock(m_work);

  for( int i=0;i<m_RT.size();++i )
    delete m_RT[i];
  m_RT.clear();
}

void SEMFEMLaplacianPreconditioner::evaluate(basics::matrixStack& res,
    const basics::matrixStack& u) const
{
  m_geometry.fineToCoarse(*m_coarse,u,m_group[m_rank],
      m_bc==SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN);
  /* local solves */
  int max=m_group[m_rank].size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    if( m_group[m_rank][i].type == basics::coarseDescriptor::FULL_OVERLAP )
      m_SP->solve((*m_coarse)[i],(*m_coarse2)[i],m_nu,m_Lx[m_group[m_rank][i].size3],
          m_Ly[m_group[m_rank][i].size3]);
    else
      m_SP2->solve((*m_coarse)[i],(*m_coarse2)[i],m_nu,m_Lx[m_group[m_rank][i].size3],
          m_Ly[m_group[m_rank][i].size3]);
  }
  m_geometry.dssumCoarse(*m_coarse,m_group[m_rank],m_rank,m_size,m_tag);

  if( m_geometry.hasCoarseSolver() ) { // coarse solve
    assert( m_nosum );
    m_geometry.fineToCoarseRestriction(*m_coarse2,res,*m_nosum,m_group[m_rank],
        m_bc==SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN,
        m_rank,m_size);

    /* restrict */
    m_geometry.getRestriction(*m_restricted,*m_work,*m_LG,*m_coarse2,*m_restrict,
        m_group[m_rank],m_RT);

    /* coarse solve */
    m_A->LLSolve(*m_restricted1);
    if( m_bc == SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN && m_nu < 1.e-12 )
      (*m_restricted)[m_restricted->length()-1] = 0;

    /* prolong */
    m_geometry.getProlongiation(*m_coarse2,*m_work,*m_LG,*m_restricted,*m_restrict,m_RT);

    *m_coarse += *m_coarse2;
  }
  m_geometry.coarseToFine(res,*m_coarse,m_group[m_rank],
      m_bc==SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN);

  if( m_bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
    m_geometry.mask(res,m_rank,m_size);
}

void SEMFEMLaplacianPreconditioner::
setupRestrictionOperators(vector<basics::Matrix*>& result,
    const basics::geometryStack& geometry,
    const basics::Vector& grid,
    int size1, int size2, 
    int size3, int size4,
    SEMLaplacianEvaluator::BC bc)
{
  std::vector<const basics::Vector*> foo;
  basics::coarseDescriptor desc = geometry.getRestrictionGridInfo();
  basics::Vector coarseGrid("coarse",desc.size1);
  utilities::GLL::GaussLobattoLegendreGrid(coarseGrid);
  coarseGrid += 1;
  coarseGrid /= 2;
  basics::Vector Lx("Lx",size1);
  basics::Vector Ly("Ly",size2);
  Lx = Real(1)/size1;
  Ly = Real(1)/size2;
  poissonSolver bar(poissonSolver::Homogenous,grid,Lx,Ly);
  basics::Vector Lx2("Lx",size3);
  basics::Vector Ly2("Ly",size4);
  Lx2 = Real(1)/size3;
  Ly2 = Real(1)/size4;
  poissonSolver bar2(poissonSolver::Nonhomogenous,grid,Lx2,Ly2);
  basics::Vector foo1 =  bar.gridX()[utilities::Range::colon(1,bar.gridX().length()-2)];
  basics::Vector foo2 =  bar.gridY()[utilities::Range::colon(1,bar.gridY().length()-2)];
  basics::Vector foo3 = bar2.gridX()[utilities::Range::colon(1,bar2.gridX().length()-2)];
  int starten=1;
  if( bc == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM ||
      bc == SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN )
    starten = 0;
  basics::Vector foo4 = bar2.gridY()[utilities::Range::colon(starten,bar2.gridY().length()-2)];
  foo.push_back(&foo1);
  foo.push_back(&foo2);
  foo.push_back(&foo3);
  foo.push_back(&foo4);
  basics::Vector temp(grid);
  temp += 1;
  temp /= 2;
  foo.push_back(&temp);
  geometry.setupRestrictionOperators(result,coarseGrid,foo);
}

SEMLaplacianDiagonalPreconditioner::SEMLaplacianDiagonalPreconditioner(const basics::geometryStack& geometry, 
    const Real nu, 
    const legendreLegendreW::poissonSolver& SP) :
  m_geometry(geometry), m_nu(nu), m_iHD("inverse diagonal helmholtz")
{
  m_iHD.add(SP.Dx(),m_geometry.size());

  basics::Matrix* D2 = utilities::g_manager.clone(SP.Dx());
  basics::Matrix* W  = utilities::g_manager.clone(SP.Dx());
  basics::Matrix* W1 = utilities::g_manager.clone(SP.Dx());
  basics::Matrix* W2 = utilities::g_manager.clone(SP.Dx());
  basics::Matrix* W3 = utilities::g_manager.clone(SP.Dx());

  basics::multPointwise(*D2,SP.Dx(),SP.Dx());

  for( int i=0;i<m_geometry.size();++i ) {
    const basics::matrixStack& G = m_geometry[i].getGH().getGeometryDerivatives();
    const basics::Matrix& J = m_geometry[i].getGH().getJacobian();
    for( int j=0;j<W1->cols();++j ) {
      for( int k=0;k<W1->rows();++k ) {
        (*W)[j][k] =   SP.weightX()[j]*SP.weightX()[k]/J[j][k];
        (*W1)[j][k] =   (*W)[j][k]*(pow(G[3][j][k],2)+
            pow(-G[1][j][k],2));
        (*W3)[j][k] =   (*W)[j][k]*(pow(-G[2][j][k],2)+
            pow(G[0][j][k],2));

        (*W2)[j][k] = 2*(*W)[j][k]*( G[3][j][k]*-G[2][j][k]+
            -G[1][j][k]* G[0][j][k]);
      }
    }

    basics::multTranspose(m_iHD[i],*D2,*W1,'T','N');
    basics::multTranspose(m_iHD[i],*W3,*D2,'N','N',mnlRealOne,mnlRealOne);
    for( int j=0;j<W1->cols();++j )
      for( int k=0;k<W1->rows();++k )
        m_iHD[i][j][k] += SP.Dx()[j][j]*SP.Dx()[k][k]*(*W2)[j][k] +
          m_nu*SP.weightX()[j]*SP.weightX()[k];

  }

  m_geometry.dssum(m_iHD);
  for( int i=0;i<m_geometry.size();++i )
    for( int j=0;j<m_iHD[i].cols();++j )
      for( int k=0;k<m_iHD[i].rows();++k )
        m_iHD[i][j][k] = Real(1)/m_iHD[i][j][k];

  utilities::g_manager.unlock(D2);
  utilities::g_manager.unlock(W);
  utilities::g_manager.unlock(W1);
  utilities::g_manager.unlock(W2);
  utilities::g_manager.unlock(W3);
}

void SEMLaplacianDiagonalPreconditioner::evaluate(basics::matrixStack& res, const basics::matrixStack& u) const
{
  for( int i=0;i<res.size();++i )
    basics::multPointwise(res[i],u[i],m_iHD[i]);
}

SEMTensoredLaplacianPreconditioner::SEMTensoredLaplacianPreconditioner(const basics::geometryStack& geometry,
    const Real nu,
    const legendreLegendreW::poissonSolver& SP) :
  m_eval(geometry,nu,SP)
{
  for( int i=0;i<geometry.size();++i ) {
    deformedTensoredLaplacianPreconditioner* pre =
      new deformedTensoredLaplacianPreconditioner(geometry[i],SP,nu);
    m_pre.push_back(pre);
  }
}

SEMTensoredLaplacianPreconditioner::~SEMTensoredLaplacianPreconditioner()
{
  for( int i=0;i<m_pre.size();++i )
    delete m_pre[i];
}

void SEMTensoredLaplacianPreconditioner::evaluate(basics::matrixStack& res, const basics::matrixStack& u) const
{
  m_eval.evaluate(res,u);
  for( int i=0;i<res.size();++i )
    m_pre[i]->evaluate(res[i],u[i]);
}

extrudedLaplacianEvaluator::extrudedLaplacianEvaluator(SEMLaplacianEvaluator& eval, 
    const legendreLegendreW::poissonSolver& SP,
    const Real nu, const basics::Matrix* Az,
    bool preconditioned,
    int rank, int size,
    SEMLaplacianEvaluator::BC bc) :
  m_eval(eval), m_SP(SP), m_nu(nu), m_Az(Az), m_nosum(NULL),
  m_rank(rank), m_size(size), m_bc(bc), m_stat3(NULL)
{
  m_division = eval.m_geometry.getDivisionInfo(size);

  for( int i=0;i<3;++i )
    m_stat.push_back(new utilities::iterationStatistics(SP.eigX().length(),size,m_scount,m_sdispl,m_rcount,m_rdispl));
  m_stat[0]->getDisplacements2(m_scount,m_sdispl,m_rcount,m_rdispl,m_rank,m_size,m_SP.Dx().length());
  m_counts = m_stat[0]->getCounts2(m_size);

  int max = 1;
#ifdef OPENMP
  max = omp_get_max_threads();
#endif
  for( int i=0;i<max;++i )
    m_evals.push_back(new SEMLaplacianEvaluator(m_eval.m_geometry,
          m_eval.m_D,
          m_eval.m_weight,
          &m_eval.m_G,
          0,1,m_eval.m_bc));

  if( preconditioned ) {
    basics::Vector grid("GLL grid",m_SP.Dx().cols());
    utilities::GLL::GaussLobattoLegendreGrid(grid);
    for( int l=0;l<m_SP.eigX().length();++l ) {
      Real nu2 = m_SP.eigX()[l]+m_nu;
#ifdef USE_FINE_PRECONDITIONER
      SEMTensoredLaplacianPreconditioner* foo
        = new SEMTensoredLaplacianPreconditioner(m_eval.m_geometry,SP.eigX()[l]+m_nu,SP,0,1);
#else
      SEMFEMLaplacianPreconditioner* foo 
        = new SEMFEMLaplacianPreconditioner(m_eval.m_geometry,grid,SP.weightX(),nu2,0,1,m_eval.m_bc);
      foo->m_tag = l+1;
#endif
      m_pre.push_back(foo);
    }
  }
}

extrudedLaplacianEvaluator::~extrudedLaplacianEvaluator()
{
  for( int i=0;i<m_pre.size();++i )
    delete m_pre[i];
  for( int i=0;i<m_evals.size();++i )
    delete m_evals[i];
  m_stat[0]->cleanDisplacements(m_scount,m_sdispl,m_rcount,m_rdispl);
  for( int i=0;i<3;++i )
    delete m_stat[i];
  if( m_stat3 ) {
    m_stat3->cleanDisplacements(m_scount3,m_sdispl3,m_rcount3,m_rdispl3);
    delete m_stat3;
  }
}

void extrudedLaplacianEvaluator::evaluate(basics::matricesStack& res, 
    const basics::matricesStack& u, bool bMask, bool doSum) const
{
  assert( m_Az );
  /* res = Az*u */
  int max=u.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i )
    applyLocalGlobal(res[i],u[i],*m_Az,'N','T');

#ifdef HAS_MPI
  vector<int*> scount, sdispl, rcount, rdispl;
  utilities::iterationStatistics stat(res[0].matrices(),m_size,
      scount,rcount,sdispl,rdispl);
  stat.getDisplacements(scount,sdispl,rcount,rdispl,m_rank,m_size,res[0][0].length());
  basics::matricesStack* buffer = 
    utilities::g_manager.aquireMatricesStack("temp",res[0].rows(),res[0].cols(),
        (rdispl[0][m_size-1]+rcount[0][m_size-1])/
        res[0][0].length(),res.size());
  vector<int> start = stat.getStarts(m_size);
  basics::matricesStack* buffer2 = utilities::g_manager.clone(*buffer);
  vector<basics::matrixStack*> stackRes = sendAndSetupStack(*buffer,res,SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN,m_rank,m_size,scount[0],sdispl[0],rcount[0],rdispl[0]);
  vector<basics::matrixStack*> stackU = sendAndSetupStack(*buffer2,u,SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN,m_rank,m_size,scount[0],sdispl[0],rcount[0],rdispl[0]);
  basics::matricesStack* temp  = utilities::g_manager.aquireMatricesStack("temp",res[0].rows(),res[0].cols(),stackRes.size(),m_eval.m_geometry.size());
#pragma omp parallel for schedule(static)
  for( int l=0;l<stackRes.size();++l ) {
#ifdef OPENMP
    int eval = omp_get_thread_num();
#else
    int eval = 0;
#endif
    m_evals[eval]->m_nu = 0;
    m_evals[eval]->m_nosum = NULL;
    /* from Z direction */
    m_evals[eval]->m_geometry.mass(*stackRes[l],false);
    m_evals[eval]->evaluate(temp->at(l),*stackU[l],false,false);

    stackRes[l]->axpy(m_eval.m_weight[l+start[m_rank]]*m_eval.m_geometry.m_Lz/2,temp->at(l));

    if( m_nu > 1.e-14 ) {
      temp->at(l) = *stackU[l];
      m_eval.m_geometry.mass(temp->at(l),false);
      stackRes[l]->axpy(m_nu*m_eval.m_weight[l+start[m_rank]]*
          m_eval.m_geometry.m_Lz/2,temp->at(l));
    }
  }
  sendStack(res,*buffer,scount[1],sdispl[1],rcount[1],rdispl[1],
      SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN,&stackRes);
  stat.cleanDisplacements(scount,sdispl,rcount,rdispl);
  for( int l=0;l<stackU.size();++l )
    delete stackU[l];

  utilities::g_manager.unlock(temp);
  utilities::g_manager.unlock(buffer);
  utilities::g_manager.unlock(buffer2);
#else
  basics::matricesStack* temp = utilities::g_manager.clone(u);
  max=u[0].matrices();
#pragma omp parallel for schedule(static)
  for( int l=0;l<max;++l ) {
#ifdef OPENMP
    int eval = omp_get_thread_num();
#else
    int eval = 0;
#endif
    m_evals[eval]->m_nu = 0;
    m_evals[eval]->m_nosum = NULL;

    /* from Z direction */
    m_evals[eval]->m_geometry.mass(res.at(l),false);
    m_evals[eval]->evaluate(temp->at(l),u.at(l),false,false);

    res.at(l).axpy(m_eval.m_weight[l]*m_eval.m_geometry.m_Lz/2,temp->at(l));
  }
  if( m_nu > 1.e-14 ) {
    *temp = u;
    m_eval.m_geometry.mass(*temp,false);
    res.axpy(m_nu,*temp);
  }
  utilities::g_manager.unlock(temp);
#endif

  if( bMask && m_bc != SEMLaplacianEvaluator::PERIODIC )
    mask(res,m_eval.m_geometry,m_rank,m_size,m_bc);

  if( m_bc == SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN &&
      m_nu < 1.e-12)
    res -= res.sum()/res.length();

  if( doSum && m_nosum )
    *m_nosum = res;

  if( doSum ) {
    if( m_bc == SEMLaplacianEvaluator::PERIODIC 		  || 
        m_bc == SEMLaplacianEvaluator::PERIODIC_DIRICHLET ||
        m_bc == SEMLaplacianEvaluator::PERIODIC_DIRICHLET_BOTTOM )
      m_eval.m_geometry.periodicDssum(res);
    else
      dssum(res,m_eval.m_geometry,m_rank,m_size);
  }
}

int applyEigenMatrix(basics::matricesStack& res, const basics::matricesStack& u,
    const basics::Matrix& Q, const SEMLaplacianEvaluator::BC bc,
    const char trans)
{
  int skip=((bc==SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET||bc==SEMLaplacianEvaluator::PERIODIC_DIRICHLET)?2:0);
  skip=(bc==SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM?1:skip);
  skip=((bc==SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM||bc==SEMLaplacianEvaluator::PERIODIC_DIRICHLET_BOTTOM)?1:skip);
  int max=u.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i )
    applyLocalGlobal(res[i],u[i],Q,'N',trans,skip);

  return skip;
}

void extrudedLaplacianEvaluator::solve(basics::matricesStack& res, basics::matricesStack& u, int c)
{
  basics::matricesStack* buffer=NULL;
  basics::matricesStack* buffer2=NULL;
  if( m_pre.size() )
    buffer = utilities::g_manager.clone(u);

  /* pre-transform */
  if( m_pre.size() )
    applyEigenMatrix(*buffer,u,*m_SP.Qx(),m_bc,'N');
  int skip=applyEigenMatrix(u,res,*m_SP.Qx(),m_bc,'N');

  Real tol = 1.e-10;
  vector<basics::matrixStack*> stack;
  vector<basics::matrixStack*> pre;
  m_counts = m_stat[c]->getCounts2(m_size);
  //    vector<int> start = m_stat[c]->getStarts(m_size);
#ifdef HAS_MPI
  m_stat[c]->getDisplacements2(m_scount,m_sdispl,m_rcount,m_rdispl,m_rank,m_size,res[0][0].length());

  basics::matricesStack* temp = 
    utilities::g_manager.aquireMatricesStack("temp",res[0].rows(),res[0].cols(),
        (m_rdispl[0][m_size-1]+m_rcount[0][m_size-1])/(res[0][0].length()),res.size());

  stack = sendAndSetupStack(*temp,u,m_bc,m_rank,m_size,m_scount[0],m_sdispl[0],
      m_rcount[0],m_rdispl[0],&m_counts.second);
  if( m_pre.size() ) {
    buffer2 = utilities::g_manager.clone(*temp);
    pre  = sendAndSetupStack(*buffer2,*buffer,m_bc,m_rank,m_size,
        m_scount[0],m_sdispl[0],m_rcount[0],
        m_rdispl[0],&m_counts.second);
  }
#else
  for( int l=0;l<m_counts.second[m_rank].size();++l ) {
    int k=m_counts.second[m_rank][l]+(skip>0?1:0);
    stack.push_back(&u.at(k));
    if( m_pre.size() )
      pre.push_back(&buffer->at(k));
  }
#endif
  m_stat[c]->reset();
  int max=stack.size();
  basics::geometryStack::geometryDotter dotter = 
    m_eval.m_geometry.getDotter(m_evals[0]->m_bc==SEMLaplacianEvaluator::PERIODIC?true:false);
#pragma omp parallel for schedule(dynamic)
  for( int l=0;l<max;++l ) {
    int iter;
#ifdef OPENMP
    int eval = omp_get_thread_num();
#else
    int eval = 0;
#endif
    int plane = m_counts.second[m_rank][l];
    //        int plane = start[m_rank]+l;
    Real nu = m_SP.eigX()[plane]+m_nu;
    m_evals[eval]->m_nu = nu;
    if( m_pre.size() ) {
      m_evals[eval]->m_nosum = 
        ((SEMFEMLaplacianPreconditioner*)m_pre[plane])->m_nosum 
        = pre[l];
      iter = utilities::CGSolver::solve(*stack[l],*m_evals[eval],
          *m_pre[plane],
          dotter,tol);
    } else {
      m_evals[eval]->m_nosum = NULL;
      iter = utilities::CGSolver::solve(*stack[l],*m_evals[eval],
          dotter,tol);
    }

    m_stat[c]->get()[plane] = iter;
    cout << "iterations " << m_name << " (" 
      << (skip>0?1:0)+plane << ") " << iter << endl;
  }
#ifdef HAS_MPI
  sendStack(u,*temp,m_scount[1],m_sdispl[1],m_rcount[1],
      m_rdispl[1],m_bc,&stack,&m_counts.second);
  m_stat[c]->exchange();
  for( int l=0;l<pre.size();++l )
    delete pre[l];
  pre.clear();
  utilities::g_manager.unlock(temp);
#endif

  /* post-transform */
  applyEigenMatrix(res,u,*m_SP.Qx(),m_bc,'T');

  utilities::g_manager.unlock(buffer);
  utilities::g_manager.unlock(buffer2);
}

pair<int,int> getComponent(int l, int skipx, int skipy, int skipz,
    int sizex, int sizey, int sizez)
{
  pair<int,int> result;
  result.first=l;
  result.second=0;

  if( result.first >= sizex) {
    result.first -= sizex;
    ++result.second;
  } 
  if( result.first >= sizey ) {
    result.first -= sizey;
    ++result.second;
  }

  if( result.second == 0 )
    result.first += (skipx>0?1:0);
  if( result.second == 1 )
    result.first += (skipy>0?1:0);
  if( result.second == 2 )
    result.first += (skipz>0?1:0);

  return( result );
}

void extrudedLaplacianEvaluator::solve(basics::Field3<basics::matricesStack>& res,
    basics::Field3<basics::matricesStack>& u,
    extrudedLaplacianEvaluator& X,
    extrudedLaplacianEvaluator& Y,
    extrudedLaplacianEvaluator& Z)
{
  if( !m_stat3 ) {
    int len=X.m_SP.eigX().length()+Y.m_SP.eigX().length()+Z.m_SP.eigX().length();
    m_stat3 = new utilities::iterationStatistics(len,X.m_size,m_scount3,m_sdispl3,m_rcount3,m_rdispl3);
    m_stat3->getDisplacements2(m_scount3,m_sdispl3,m_rcount3,m_rdispl3,X.m_rank,X.m_size,X.m_SP.Dx().length());
    m_counts3 = m_stat3->getCounts2(m_size);
  }
  vector<extrudedLaplacianEvaluator*> solver;
  solver.push_back(&X);
  solver.push_back(&Y);
  solver.push_back(&Z);
  basics::Field3<basics::matricesStack>* buffer=NULL;
  basics::matricesStack* buffer2=NULL;
  basics::matricesStack* buffer3=NULL;
  if( m_pre.size() )
    buffer = utilities::g_manager.clone(u);

  /* pre-transform */
  int skip[3];
  for( int i=0;i<3;++i ) {
    if( m_pre.size() )
      applyEigenMatrix((*buffer)[i],u[i],*solver[i]->m_SP.Qx(),solver[i]->m_bc,'N');
    skip[i]=applyEigenMatrix(u[i],res[i],*solver[i]->m_SP.Qx(),solver[i]->m_bc,'N');
  }

  vector<basics::matrixStack*> stack;
  vector<basics::matrixStack*> pre;
  m_counts3 = m_stat3->getCounts2(m_size);
  //    vector<int> start = m_stat3->getStarts(m_size);

#ifdef HAS_MPI
  int len = res.X()[0][0].length();
  m_stat3->getDisplacements2(m_scount3,m_sdispl3,m_rcount3,m_rdispl3,m_rank,m_size,len);

  basics::matricesStack* temp = 
    utilities::g_manager.aquireMatricesStack("temp",res[0][0].rows(),res[0][0].cols(),
        (m_rdispl3[0][m_size-1]+m_rcount3[0][m_size-1])/len,res.X().size());
  basics::matricesStack* temp2 = 
    utilities::g_manager.aquireMatricesStack("temp2",res[0][0].rows(),res[0][0].cols(),
        (m_sdispl3[0][m_size-1]+m_scount3[0][m_size-1])/len,res.X().size());
  stack = sendAndSetupStack(*temp,*temp2,u,
      solver[0]->m_bc,solver[1]->m_bc,solver[2]->m_bc,
      m_rank,m_size,m_scount3[0],m_sdispl3[0],
      m_rcount3[0],m_rdispl3[0],&m_counts3.second);
  if( m_pre.size() ) {
    buffer2 = utilities::g_manager.clone(*temp);
    buffer3 = utilities::g_manager.clone(*temp2);
    pre  = sendAndSetupStack(*buffer2,*buffer3,*buffer,
        solver[0]->m_bc,solver[1]->m_bc,solver[2]->m_bc,
        m_rank,m_size,m_scount3[0],m_sdispl3[0],m_rcount3[0],
        m_rdispl3[0],&m_counts3.second);
  }
#else
  for( int l=0;l<m_counts3.second[m_rank].size();++l ) {
    pair<int,int> plane = getComponent(m_counts3.second[m_rank][l],
        skip[0],skip[1],skip[2],
        solver[0]->m_SP.eigX().length(),
        solver[1]->m_SP.eigX().length(),
        solver[2]->m_SP.eigX().length());
    stack.push_back(&(u[plane.second].at(plane.first)));
    if( m_pre.size() )
      pre.push_back(&((*buffer)[plane.second].at(plane.first)));
  }
#endif
  Real tol = 1.e-8;
  m_stat3->reset();
  int max=stack.size();
  basics::geometryStack::geometryDotter dotter = 
    m_eval.m_geometry.getDotter(m_evals[0]->m_bc==SEMLaplacianEvaluator::PERIODIC?true:false);
#pragma omp parallel for schedule(dynamic)
  for( int l=0;l<max;++l ) {
    int iter;
#ifdef OPENMP
    int eval = omp_get_thread_num();
#else
    int eval = 0;
#endif
    int plan = m_counts3.second[m_rank][l];
    //        int plan = start[m_rank]+l;
    pair<int,int> plane = getComponent(plan,0,0,0,
        solver[0]->m_SP.eigX().length(),
        solver[1]->m_SP.eigX().length(),
        solver[2]->m_SP.eigX().length());
    Real nu = solver[plane.second]->m_SP.eigX()[plane.first]+m_nu;
    m_evals[eval]->m_nu = nu;
    if( m_pre.size() ) {
      m_evals[eval]->m_nosum = 
        ((SEMFEMLaplacianPreconditioner*)solver[plane.second]->m_pre[plane.first])->m_nosum 
        = pre[l];
      iter = utilities::CGSolver::solve(*stack[l],*m_evals[eval],
          *solver[plane.second]->m_pre[plane.first],
          dotter,tol);
    } else {
      m_evals[eval]->m_nosum = NULL;
      iter = utilities::CGSolver::solve(*stack[l],*m_evals[eval],
          dotter,tol);
    }

    //        m_stat3->get()[start[m_rank]+l] = iter;
    m_stat3->get()[m_counts3.second[m_rank][l]] = iter;
    cout << "iterations component " << plane.second << " (" 
      << plane.first+(skip[plane.second]?1:0) << ") " << iter << endl;
  }
#ifdef HAS_MPI
  sendStack(u,*temp2,*temp,m_scount3[1],m_sdispl3[1],m_rcount3[1],
      m_rdispl3[1],solver[0]->m_bc,solver[1]->m_bc,solver[2]->m_bc,
      &stack,&m_counts3.second);
  m_stat3->exchange();
  for( int l=0;l<pre.size();++l )
    delete pre[l];
  pre.clear();
  utilities::g_manager.unlock(temp);
  utilities::g_manager.unlock(temp2);
#endif

  /* post-transform */
  for( int i=0;i<3;++i )
    applyEigenMatrix(res[i],u[i],*solver[i]->m_SP.Qx(),solver[i]->m_bc,'T');

  utilities::g_manager.unlock(buffer);
  utilities::g_manager.unlock(buffer2);
  utilities::g_manager.unlock(buffer3);
}

extrudedFEMLaplacianPreconditioner::
extrudedFEMLaplacianPreconditioner(basics::geometryStack& geometry,
    const basics::Vector& weight,
    const basics::Vector& grid,
    Real nu, int rank, int size,
    SEMLaplacianEvaluator::BC bc) :
  m_geometry(geometry), m_rank(rank), m_size(size), m_nu(nu), m_bc(bc),
  m_SP(NULL), m_SP2(NULL), m_coarse(NULL), m_coarse2(NULL), m_work3(NULL),
  m_restricted1(NULL)
{
  m_group = geometry.getCoarseGroups(size);
  basics::coarseGrid group = geometry.getCoarseGroups(1)[0];
  /* based on the first FULL OVERLAP part */
  int n=0;
  while( n < group.size() && group[n].type != basics::coarseDescriptor::FULL_OVERLAP )
    ++n;

  basics::Vector Lx("Lx",group[n].size1);
  basics::Vector Ly("Ly",group[n].size2);

  for( int j=0;j<group[n].size1;++j )
    Lx[j] = m_geometry[group[n].elements[j]].getAdjustedSize(weight,grid).first;
  for( int j=0;j<group[n].size2;++j )
    Ly[j] = m_geometry[group[n].elements[group[n].size1-1+
      j*group[n].size1]].getAdjustedSize(weight,grid).second;

  poissonSolver::BC bcp = poissonSolver::Homogenous;
  poissonSolver::BC bcp2 = poissonSolver::HomogenousLeft;
  bool neumann=false;
  if( m_bc == SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN ) {
    bcp2 = bcp = poissonSolver::Nonhomogenous;
    neumann = true;
  }
  if( m_bc == SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM ||
      m_bc == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM)
  {
    bcp = poissonSolver::HomogenousLeft;
    if( m_bc == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM ) {
      bcp2 = poissonSolver::Nonhomogenous;
      neumann = true;
    }
  }
  m_SP = new poissonSolver(poissonSolver::Homogenous,grid,Lx,Ly,false,true,bcp);

  int n2 = 0;
  while( n2 < group.size() && 
      group[n2].type == basics::coarseDescriptor::FULL_OVERLAP )
    ++n2;

  basics::Vector Lx2("Lx2",group[n2].size1);
  basics::Vector Ly2("Ly2",group[n2].size2);

  for( int j=0;j<group[n2].size1;++j )
    Lx2[j] = m_geometry[group[n2].elements[j]].getAdjustedSize(weight,grid).first;
  for( int j=0;j<group[n2].size2;++j )
    Ly2[j] = m_geometry[group[n2].elements[group[n2].size1-1+
      j*group[n2].size1]].getAdjustedSize(weight,grid).second;

  m_SP2 = new poissonSolver(bcp2,grid,Lx2,Ly2,false,true,bcp);

  for( int c=0;c<group.size();++c ) {
    Real lx=0;
    Real ly=0;
    Real area=0;

    if( group[c].type2 == basics::coarseDescriptor::FLIP_XY ) {
      for( int j=0;j<group[c].size1;++j )
        lx += m_geometry[group[c].elements[j]].getSize(weight,grid).second;
      for( int j=0;j<group[c].size2;++j )
        ly += m_geometry[group[c].elements[group[c].size1-1+
          j*group[c].size1]].getSize(weight,grid).first;
    } else {
      for( int j=0;j<group[c].size1;++j )
        lx += m_geometry[group[c].elements[j]].getSize(weight,grid).first;
      for( int j=0;j<group[c].size2;++j )
        ly += m_geometry[group[c].elements[group[c].size1-1+
          j*group[c].size1]].getSize(weight,grid).second;
    }
    for( int i=0;i<group[c].elements.size();++i )
      area += m_geometry[group[c].elements[i]].getVolume(weight);

    Real scale = sqrt(area/(lx*ly));
    lx *= scale;
    ly *= scale;

    m_Lx.push_back(lx);
    m_Ly.push_back(ly);
  }
  m_coarse  = m_geometry.getCoarseBuffer(grid.length(),grid.length(),grid.length(),m_group[m_rank],"coarse",false,neumann);
  m_coarse2 = m_geometry.getCoarseBuffer(grid.length(),grid.length(),grid.length(),m_group[m_rank],"coarse2",false,neumann);
  m_work3 = m_geometry.getCoarseBuffer(grid.length(),grid.length(),grid.length(),group,"coarse3",false,neumann);

  if( m_geometry.hasCoarseSolver() ) {
    basics::coarseDescriptor desc = m_geometry.getRestrictionGridInfo();
    if( desc.size1 > 0 ) {
      m_work  = utilities::g_manager.aquireMatricesStack("working buffer",
          desc.size1,
          desc.size1,
          desc.size1,
          desc.size3);
      m_work2 = m_geometry.getCoarseBuffer(grid.length(),grid.length(),
          desc.size1,m_group[m_rank],
          "working buffer 2",false,neumann);
      m_work3 = m_geometry.getCoarseBuffer(grid.length(),grid.length(),
          desc.size1,group,
          "working buffer 3",false,neumann);

      HDF5::HDF5Reader reader(m_geometry.getOperatorFile());
      m_LG = utilities::g_manager.aquireMatricesStack("local to global",
          desc.size1,
          desc.size1,
          desc.size1,
          desc.size3);
      if( bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET)  {
        basics::Matrix* A = new basics::Matrix("helmholtz",desc.size2,desc.size2);
        reader.read(*m_LG,"LG3d");
        reader.read(*A,"A3");
        m_A = new basics::Matrix(*A);
        reader.read(*A,"B3");
        m_A->axpy(nu,*A); // helmholtz factor
        delete A;
        m_restricted = utilities::g_manager.aquireVector("restricted",desc.size2);
      } 
      if( bc == SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM) {
        basics::Matrix* A = new basics::Matrix("helmholtz",2*desc.size2,
            2*desc.size2);
        reader.read(*m_LG,"LG32d");
        reader.read(*A,"A32");
        m_A = new basics::Matrix(*A);
        reader.read(*A,"B32");
        m_A->axpy(nu,*A); // helmholtz factor
        delete A;
        m_restricted = utilities::g_manager.aquireVector("restricted",desc.size2*2);
      }
      if( bc == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM) {
        basics::Matrix* A = new basics::Matrix("helmholtz",2*desc.size4,
            2*desc.size4);
        reader.read(*m_LG,"LG32n");
        reader.read(*A,"A32n");
        m_A = new basics::Matrix(*A);
        reader.read(*A,"B32n");
        m_A->axpy(nu,*A); // helmholtz factor
        delete A;
        m_restricted = utilities::g_manager.aquireVector("restricted",desc.size4*2);
      }
      if( bc == SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN ) {
        basics::Matrix* A = new basics::Matrix("helmholtz",desc.size4*desc.size1,desc.size4*desc.size1);
        reader.read(*m_LG,"LG3n");
        reader.read(*A,"A3n");
        m_A = new basics::Matrix(*A);
        reader.read(*A,"B3n");
        m_A->axpy(nu,*A); // helmholtz factor
        delete A;
        m_restricted = utilities::g_manager.aquireVector("restricted",m_A->rows());
        if( m_nu < 1.e-12 ) {
          m_restricted1 = new basics::Vector("restricted minus one",m_restricted->length()-1,false,m_restricted->data());
          A = m_A;
          m_A = new basics::Matrix(m_A->submatrix(utilities::Range::colon(0,A->rows()-2),utilities::Range::colon(0,A->cols()-2)));
          delete A;
        }
      }

      *m_A = m_A->choleskyFactorize();

      SEMFEMLaplacianPreconditioner::
        setupRestrictionOperators(m_RT,geometry,grid,
            group[ n].size1,group[ n].size2,
            group[n2].size1,group[n2].size2,bc);
      m_restrict  = utilities::g_manager.aquireMatrix("restrict buffer2",(*m_coarse)[0].rows(),desc.size1);
    }
  }
}

extrudedFEMLaplacianPreconditioner::~extrudedFEMLaplacianPreconditioner()
{
  delete m_SP;
  delete m_SP2;

  delete m_coarse;
  delete m_coarse2;
}

void extrudedFEMLaplacianPreconditioner::evaluate(basics::matricesStack& res,
    const basics::matricesStack& u) const
{
  bool neumann = (m_bc == SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN ||
      m_bc == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM);
  m_geometry.fineToCoarse(*m_coarse,u,m_group[m_rank],neumann);
  /* local solves */
  int max=m_group[m_rank].size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    if( m_group[m_rank][i].type == basics::coarseDescriptor::FULL_OVERLAP )
      m_SP->solve((*m_coarse)[i],(*m_coarse2)[i],
          m_Lx[m_group[m_rank][i].size3],
          m_Ly[m_group[m_rank][i].size3],
          m_geometry.m_Lz,m_nu);
    else
      m_SP2->solve((*m_coarse)[i],(*m_coarse2)[i],
          m_Lx[m_group[m_rank][i].size3],
          m_Ly[m_group[m_rank][i].size3],
          m_geometry.m_Lz,m_nu);
  }
  dssumCoarse(*m_coarse,res[0].rows(),m_geometry,m_rank,m_size,neumann);

  if( m_geometry.hasCoarseSolver() ) { // coarse solve
    assert( m_nosum );
    m_geometry.fineToCoarseRestriction(*m_coarse2,res,*m_nosum,
        m_group[m_rank],neumann,m_rank,m_size);
    /* restrict */
    m_geometry.gatherAndRestrict(*m_restricted,*m_work,*m_work2,
        *m_work3,*m_LG,*m_coarse2,*m_restrict,m_RT,
        m_rank,m_size);
    /* coarse solve */
    if( m_rank == 0 ) {
      if( m_bc != SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN || m_nu > 1.e-12)
        m_A->LLSolve(*m_restricted);
      else {
        m_A->LLSolve(*m_restricted1);
        (*m_restricted)[m_restricted->length()-1] = 0;
      }
    }
    /* prolong */
    m_geometry.prolongAndScatter(*m_coarse2,*m_work,*m_work2,*m_work3,
        *m_LG,*m_restricted,*m_restrict,m_RT,
        m_rank,m_size);

    *m_coarse += *m_coarse2;
  }
  m_geometry.coarseToFine(res,*m_coarse,m_group[m_rank],neumann);
  mask(res,m_geometry,m_rank,m_size,m_bc);
  //    if( m_bc == SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN && m_nu < 1.e-12 )
  //        res -= res.sum()/res.length();
}

extrudedLaplacianPreconditioner::extrudedLaplacianPreconditioner(extrudedLaplacianEvaluator& eval,
    const basics::geometryStack3D& geometry) :
  m_eval(eval), m_geometry(geometry)
{
}

void extrudedLaplacianPreconditioner::evaluate(mnl::basics::matricesStack& res,
    const mnl::basics::matricesStack& u) const
{
  basics::matricesStack* temp  = m_geometry.localToGlobalZ(u,m_eval.m_rank,m_eval.m_size);
  basics::matricesStack* temp2 = m_geometry.localToGlobalZ(*m_eval.m_nosum,m_eval.m_rank,
      m_eval.m_size,false,true);
  m_eval.solve(*temp,*temp2);
  m_geometry.globalToLocalZ(res,temp,m_eval.m_rank,m_eval.m_size);
  utilities::g_manager.unlock(temp2);
  mask(res,m_geometry,m_eval.m_rank,m_eval.m_size,m_eval.m_bc);
}

SEM3DLaplacianEvaluator::SEM3DLaplacianEvaluator(const basics::geometryStack3D& geometry,
    const basics::Vector& weight,
    const basics::Matrix& D, 
    Real nu, int rank, int size,
    SEMLaplacianEvaluator::BC bc) :
  m_geometry(geometry), m_rank(rank), m_bc(bc),
  m_size(size), m_D(D), m_weight(weight), m_nu(nu), m_nosum(NULL)
{
  m_division = geometry.getDivisionInfo(size);
  vector<basics::matricesStack*> derivatives = m_geometry.getReferenceGeometryDerivatives();
  for( int i=0;i<6;++i )
    m_G.push_back(new basics::matricesStack("geometry factors",D.rows(),D.cols(),
          D.rows(),m_geometry.size()));
  mnl::basics::componentView<mnl::basics::Matrices,mnl::basics::matricesStack> view2(derivatives);
  for( int i=0;i<m_geometry.size();++i) {
    const basics::Matrices& J = geometry[i].getGH().getJacobian();	
    for( int l=0;l<J.matrices();++l ) {
      for( int j=0;j<J.cols();++j ) {
        for( int k=0;k<J.rows();++k ) {
          Real Ja 	= J[l][j][k];
          Real xix    = (*derivatives[0])[i][l][j][k];
          Real xiy    = (*derivatives[1])[i][l][j][k];
          Real xiz    = (*derivatives[2])[i][l][j][k];
          Real etax   = (*derivatives[3])[i][l][j][k];
          Real etay   = (*derivatives[4])[i][l][j][k];
          Real etaz   = (*derivatives[5])[i][l][j][k];
          Real gammax = (*derivatives[6])[i][l][j][k];
          Real gammay = (*derivatives[7])[i][l][j][k];
          Real gammaz = (*derivatives[8])[i][l][j][k];

          (*m_G[0])[i][l][j][k] = 1/Ja*(pow(xix,2)+pow(xiy,2)+pow(xiz,2));
          (*m_G[1])[i][l][j][k] = 1/Ja*(xix*etax+xiy*etay+xiz*etaz);
          (*m_G[2])[i][l][j][k] = 1/Ja*(xix*gammax+xiy*gammay+xiz*gammaz);
          (*m_G[3])[i][l][j][k] = 1/Ja*(pow(etax,2)+pow(etay,2)+pow(etaz,2));
          (*m_G[4])[i][l][j][k] = 1/Ja*(etax*gammax+etay*gammay+etaz*gammaz);
          (*m_G[5])[i][l][j][k] = 1/Ja*(pow(gammax,2)+pow(gammay,2)+pow(gammaz,2));
        }
      }
    }
  }

  int N = geometry[0].getGH().getJacobian().rows();
  m_Uxi    = utilities::g_manager.aquireMatricesStack("Uxi",N,N,N,
      m_division[m_rank].elements.size());
  m_Ueta   = utilities::g_manager.clone(*m_Uxi);
  m_Ugamma = utilities::g_manager.clone(*m_Uxi);
  m_t1     = utilities::g_manager.clone(*m_Uxi);
  m_t2     = utilities::g_manager.clone(*m_Uxi);
  m_t3     = utilities::g_manager.clone(*m_Uxi);

  m_view = new mnl::basics::componentView<mnl::basics::Matrices,
         mnl::basics::matricesStack>(m_G);

  for( int i=0;i<m_division[m_rank].elements.size();++i ) {
    int eval=m_division[m_rank].elements[i];
    m_evals.push_back(new deformed3DLaplacianEvaluator((*m_view)[eval],
          m_D,m_weight,
          (*m_Uxi)[i],(*m_Ueta)[i],
          (*m_Ugamma)[i],(*m_t1)[i],
          (*m_t2)[i],(*m_t3)[i]));
  }
}

SEM3DLaplacianEvaluator::~SEM3DLaplacianEvaluator()
{
  delete m_view;
  for( unsigned int i=0;i<m_G.size();++i )
    delete m_G[i];
  for( int i=0;i<m_evals.size();++i )
    delete m_evals[i];
}

void SEM3DLaplacianEvaluator::evaluate(basics::matricesStack& res,
    const basics::matricesStack& u, bool bMask, bool doSum) const
{
  int max=m_division[m_rank].elements.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    int eval = m_division[m_rank].elements[i];
    m_evals[i]->evaluate(res[i],u[i],false);
    if( m_nu > 1.e-14 ) {
      basics::multPointwise((*m_Uxi)[i],u[i],(*m_geometry.m_mass)[eval]);
      res[i].axpy(m_nu,(*m_Uxi)[i]);
    }
  }
  if( bMask )
    mask(res,m_geometry,m_rank,m_size,m_bc);
  if( m_nosum && doSum)
    *m_nosum = res;
  if( doSum )
    dssum(res,m_geometry,m_rank,m_size);
}

SEM3DFEMLaplacianPreconditioner::
SEM3DFEMLaplacianPreconditioner(basics::geometryStack3D& geometry,
    const basics::Vector& weight,
    const basics::Vector& grid,
    Real nu, int rank, int size) :
  m_geometry(geometry), m_rank(rank), m_size(size), m_nu(nu),
  m_A(NULL), m_restrict(NULL), m_restricted(NULL), m_work(NULL),
  m_work2(NULL), m_coarse(NULL), m_coarse2(NULL)
{
  //    return;
  m_group = geometry.getCoarseGroups(size);
  basics::coarseGrid group = geometry.getCoarseGroups(1)[0];
  /* based on the first FULL OVERLAP part */
  int n=0;
  while( n < group.size() && group[n].type != basics::coarseDescriptor::FULL_OVERLAP )
    ++n;

  basics::Vector Lx("Lx",group[n].size1);
  basics::Vector Ly("Ly",group[n].size2);

  for( int j=0;j<group[n].size1;++j )
    Lx[j] = m_geometry[group[n].elements[j]].getAdjustedSize(weight,grid)[0];
  for( int j=0;j<group[n].size2;++j )
    Ly[j] = m_geometry[group[n].elements[group[n].size1-1+
      j*group[n].size1]].getAdjustedSize(weight,grid)[1];
  m_SP.push_back(new poissonSolver(poissonSolver::Homogenous,
        grid,Lx,Ly,false,true,
        poissonSolver::HomogenousLeft));
  m_SP.push_back(new poissonSolver(poissonSolver::Homogenous,
        grid,Lx,Ly,false,true,
        poissonSolver::Nonhomogenous));
  m_SP.push_back(new poissonSolver(poissonSolver::Homogenous,
        grid,Lx,Ly,false,true,
        poissonSolver::HomogenousRight));

  int n2 = 0;
  while( n2 < group.size() && group[n2].type == basics::coarseDescriptor::FULL_OVERLAP )
    ++n2;

  basics::Vector Lx2("Lx2",group[n2].size1);
  basics::Vector Ly2("Ly2",group[n2].size2);

  for( int j=0;j<group[n2].size1;++j )
    Lx2[j] = m_geometry[group[n2].elements[j]].getAdjustedSize(weight,grid)[0];
  for( int j=0;j<group[n2].size2;++j )
    Ly2[j] = m_geometry[group[n2].elements[group[n2].size1-1+
      j*group[n2].size1]].getAdjustedSize(weight,grid)[1];

  m_SP.push_back(new poissonSolver(poissonSolver::Nonhomogenous,
        grid,Lx2,Ly2,false,true,
        poissonSolver::HomogenousLeft));
  m_SP.push_back(new poissonSolver(poissonSolver::Nonhomogenous,
        grid,Lx2,Ly2,false,true,
        poissonSolver::Nonhomogenous));
  m_SP.push_back(new poissonSolver(poissonSolver::Nonhomogenous,
        grid,Lx2,Ly2,false,true,
        poissonSolver::HomogenousRight));

  for( int c=0;c<group.size();++c ) {
    Real lx=0;
    Real ly=0;
    Real lz=0;
    Real volume=0;

    for( int j=0;j<group[c].size1;++j ) {
      vector<Real> sizes = m_geometry[group[c].elements[j]].getSize(weight,grid);
      if( group[c].type2 == basics::coarseDescriptor::FLIP_XY )
        lx += sizes[1];
      else
        lx += sizes[0];
    }
    for( int j=0;j<group[c].size2;++j ) {
      vector<Real> sizes = m_geometry[group[c].elements[group[c].size1-1+
        j*group[c].size1]].getSize(weight,grid);
      if( group[c].type2 == basics::coarseDescriptor::FLIP_XY )
        ly += sizes[0];
      else
        ly += sizes[1];
    }
    for( int i=0;i<group[c].elements.size();++i ) {
      lz += m_geometry[group[c].elements[i]].getSize(weight,grid)[2];
      volume += m_geometry[group[c].elements[i]].getVolume(weight);
    }
    lz /= group[c].elements.size();

    Real scale = pow(volume/(lx*ly*lz),Real(1)/3);
    lx *= scale;
    ly *= scale;
    lz *= scale;

    m_Lx.push_back(lx);
    m_Ly.push_back(ly);
    m_Lz.push_back(lz);
  }
  m_coarse  = m_geometry.getCoarseBuffer(grid.length(),grid.length(),
      grid.length(),m_group[m_rank],"coarse");
  m_coarse2 = m_geometry.getCoarseBuffer(grid.length(),grid.length(),
      grid.length(),m_group[m_rank],"coarse2");

  if( m_geometry.hasCoarseSolver() ) {
    basics::coarseDescriptor desc = m_geometry.getRestrictionGridInfo();
    if( desc.size1 > 0 ) {
      m_work  = utilities::g_manager.aquireMatricesStack("working buffer",desc.size1,desc.size1,desc.size1,desc.size3);
      m_work2 = m_geometry.getCoarseBuffer(grid.length(),grid.length(),desc.size1,m_group[m_rank],"working buffer 2");
      basics::Matrix* A = new basics::Matrix("helmholtz",desc.size2,desc.size2);

      HDF5::HDF5Reader reader(m_geometry.getOperatorFile());
      reader.read(*A,"A3");
      m_LG = utilities::g_manager.aquireMatricesStack("local to global",desc.size1,desc.size1,desc.size1,desc.size3);
      reader.read(*m_LG,"LG3Dd");
      m_A = new basics::Matrix(*A);
      reader.read(*A,"B3");
      m_A->axpy(nu,*A); // helmholtz factor
      *m_A = m_A->choleskyFactorize();

      delete A;

      m_restricted = utilities::g_manager.aquireVector("restricted",desc.size2);
      m_restrict  = utilities::g_manager.aquireMatrix("restrict buffer2",(*m_coarse)[0].rows(),desc.size1);
    }
  } else
    m_restricted = NULL;
}

SEM3DFEMLaplacianPreconditioner::~SEM3DFEMLaplacianPreconditioner()
{
  for( int i=0;i<m_SP.size();++i )
    delete m_SP[i];

  delete m_coarse;
  delete m_coarse2;
  delete m_A;

  delete m_work2;

  utilities::g_manager.unlock(m_restricted);
  utilities::g_manager.unlock(m_restrict);
  utilities::g_manager.unlock(m_work);
}

void SEM3DFEMLaplacianPreconditioner::evaluate(basics::matricesStack& res,
    const basics::matricesStack& u) const
{
  //    utilities::g_profiler.start(2);
  m_geometry.fineToCoarse(*m_coarse,u,m_group[m_rank]);
  /* local solves */
  int max=m_group[m_rank].size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    int solver=1;
    if( i > max-32 )
      solver = 2;
    if( i < 32 )
      solver = 0;
    if( m_group[m_rank][i].type == basics::coarseDescriptor::FULL_OVERLAP )
      m_SP[solver]->solve((*m_coarse)[i],(*m_coarse2)[i],
          m_Lx[m_group[m_rank][i].size3],
          m_Ly[m_group[m_rank][i].size3],
          m_Lz[m_group[m_rank][i].size3],m_nu);
    else {
      m_SP[solver+3]->solve((*m_coarse)[i],(*m_coarse2)[i],
          m_Lx[m_group[m_rank][i].size3],
          m_Ly[m_group[m_rank][i].size3],
          m_Lz[m_group[m_rank][i].size3],m_nu);
    }
  }
  m_geometry.dssumCoarse(*m_coarse,m_group[m_rank],m_rank,m_size);

  if( m_geometry.hasCoarseSolver() ) { // coarse solve
    assert( m_nosum );
    m_geometry.fineToCoarseRestriction(*m_coarse2,res,*m_nosum,m_group[m_rank],m_rank,m_size);

    /* restrict */
    m_geometry.getRestriction(*m_restricted,*m_work,*m_work2,*m_LG,*m_coarse2,*m_restrict,m_group[m_rank]);

    /* coarse solve */
    m_A->LLSolve(*m_restricted);

    /* prolong */
    m_geometry.getProlongiation(*m_coarse2,*m_work,*m_work2,*m_LG,*m_restricted,*m_restrict);

    *m_coarse += *m_coarse2;
  }
  m_geometry.coarseToFine(res,*m_coarse,m_group[m_rank]);
  mask(res,m_geometry,m_rank,m_size,m_bc);
  //    utilities::g_profiler.pause(2);
}

