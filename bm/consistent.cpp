#include "consistent.h"

#include "mnl/hdf5.h"

#include "poissonsolver-lgw.h"

using namespace mnl;
using namespace std;
using namespace legendreLegendreW;

deformedConsistentPressureEvaluator::deformedConsistentPressureEvaluator(const deformedDivergenceEvaluator& divergence,
    const deformedGradientEvaluator& gradient,
    const basics::spectralElement2D& elem,
    const basics::Vector& weight)  :
  m_divergence(divergence), m_gradient(gradient), m_element(elem), m_weight(weight)
{
  m_W = utilities::g_manager.aquireField2D("buffer",weight.length(),weight.length());
}

deformedConsistentPressureEvaluator::~deformedConsistentPressureEvaluator()
{
  utilities::g_manager.unlock(m_W);
}

void deformedConsistentPressureEvaluator::evaluate(basics::Matrix& res,
    const basics::Matrix& p) const
{
  m_gradient.evaluate(*m_W,p);
  m_element.invMass(m_W->X(),m_weight);
  m_element.invMass(m_W->Y(),m_weight);
  m_divergence.evaluate(res,*m_W);
}

SEMConsistentPressureEvaluator::SEMConsistentPressureEvaluator(const basics::geometryStack& geometry,
    const basics::Matrix& D,
    const basics::Matrix& GLL2G,
    const basics::Vector& weightGL,
    const basics::Vector& weight,
    basics::matrixStack* J,
    const Real nu,
    SEMLaplacianEvaluator::BC bc,
    int rank, int size) :
  m_weight(weight),
  m_geometry(geometry), 
  m_divergence(geometry,D,GLL2G,weightGL,NULL,rank,size), 
  m_gradient(geometry,D,GLL2G,weightGL,NULL,rank,size),
  m_GLL2G(GLL2G), m_E(NULL), m_nu(nu), m_mine(true),
  m_rank(rank), m_size(size), m_deflated(false),
  m_restricted(NULL), m_restricted1(NULL), m_bc(bc)
{
  if( J ) {
    m_J = J;
    m_mine = false;
  }
  else {
    vector<basics::matrixStack*> G = geometry.getInterpolatedGeometryDerivatives(GLL2G);
    m_J = utilities::g_manager.aquireMatrixStack("interpolated jacobian",GLL2G.rows(),GLL2G.rows(),geometry.size());
    for( int i=0;i<m_J->size();++i ) {
      for( int j=0;j<(*m_J)[i].cols();++j )
        for( int k=0;k<(*m_J)[i].rows();++k )
          (*m_J)[i][j][k] = (((*G[0])[i][j][k]*(*G[3])[i][j][k])
              -((*G[2])[i][j][k]*(*G[1])[i][j][k]))*weightGL[j]*weightGL[k];
    }
    for( int i=0;i<G.size();++i )
      delete G[i];
  }
  m_desc     = m_geometry.getCoarseGroups(size);
  m_division = geometry.getDivisionInfo(size);
  basics::coarseDescriptor desc = m_geometry.getRestrictionGridInfo();
  m_buf  	= utilities::g_manager.aquireMatrixStack("buffer",weightGL.length(),weightGL.length(),m_division[m_rank].elements.size());
  m_buf2 	= utilities::g_manager.clone(*m_buf);
  m_buf3 	= utilities::g_manager.aquireMatrixStack("buffer",D.rows(),weightGL.length(),m_division[m_rank].elements.size());
  m_buf4 	= utilities::g_manager.clone(*m_buf);
  m_buf5 	= utilities::g_manager.clone(*m_buf);
  m_W 	= utilities::g_manager.aquireMatrixStackField("buffer",D.rows(),D.cols(),m_division[m_rank].elements.size());
  if( 0 && desc.size3 > 0 ) {
    m_restricted  = utilities::g_manager.aquireVector("restricted",desc.size3);
    m_restricted1 = new basics::Vector("restricted minus one",m_restricted->length()-1,false,m_restricted->data());
    int skip=1;
    if( m_nu > 1.e-12 )
      skip = 0;
    m_E = new basics::Matrix("coarse E",desc.size3-skip,desc.size3-skip);

    basics::matrixStack* buf3 = utilities::g_manager.clone(*m_buf2);
    for( int i=0;i<m_E->cols();++i ) {
      m_restricted->clear();
      (*m_restricted)[i] = 1;
      getProlongiation(*m_buf,*m_restricted,m_desc[m_rank],m_nu<1.e-12);
      evaluateSimple(*buf3,*m_buf);
      getRestriction(*m_restricted,*buf3,m_desc[m_rank],m_rank,m_size);
      if( m_nu > 1.e-12 )
        (*m_E)[i] = *m_restricted;
      else
        (*m_E)[i] = *m_restricted1;
    }
    utilities::g_manager.unlock(buf3);

    *m_E = m_E->choleskyFactorize();
  }

}

SEMConsistentPressureEvaluator::~SEMConsistentPressureEvaluator()
{
  delete m_E;
  if( m_mine )
    utilities::g_manager.unlock(m_J);
  delete m_restricted1;
  utilities::g_manager.unlock(m_restricted);
  utilities::g_manager.unlock(m_buf);
  utilities::g_manager.unlock(m_buf2);
  utilities::g_manager.unlock(m_buf3);
  utilities::g_manager.unlock(m_buf4);
  utilities::g_manager.unlock(m_buf5);
  utilities::g_manager.unlock(m_W);
}

void SEMConsistentPressureEvaluator::evaluate(basics::matrixStack& res,
    const basics::matrixStack& p) const
{
  evaluateSimple(res,p);
  if( m_deflated ) {
    getRestriction(*m_restricted,res,m_desc[m_rank],m_rank,m_size);
    if( m_nu > 1.e-12 )
      m_E->LLSolve(*m_restricted);
    else
      m_E->LLSolve(*m_restricted1);
    getProlongiation(*m_buf4,*m_restricted,m_desc[m_rank],m_nu<1.e-12);
    evaluateSimple(*m_buf5,*m_buf4);
    res -= *m_buf5;
    orthogonalize(res,m_desc[m_rank]);
  }
  if( m_nu < 1.e-10 )
    res -= res.sum()/res.length();
}

void SEMConsistentPressureEvaluator::evaluateSimple(basics::matrixStack& res, 
    const basics::matrixStack& p) const
{
  m_gradient.evaluate(*m_W,p);
  if( m_bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET ) {
    m_geometry.dssum(*m_W);
    m_geometry.maskField(*m_W);
    m_geometry.invMass(*m_W,m_division[m_rank].elements);
  }
  if( m_bc == SEMLaplacianEvaluator::PERIODIC ) {
    m_geometry.periodicDssum(*m_W);
    m_geometry.invMassP(m_W->X());
    m_geometry.invMassP(m_W->Y());
  }

  if( m_bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
    m_geometry.maskField(*m_W,m_rank,m_size);
  m_divergence.evaluate(res,*m_W);
  if( m_nu > 1.e-12 ) {
    int max=p.size();
#pragma omp parallel for schedule(static)
    for( int i=0;i<max;++i ) {
      basics::multPointwise((*m_buf)[i],p[i],(*m_J)[m_division[m_rank].elements[i]]);
      basics::multTranspose((*m_buf3)[i],m_GLL2G,(*m_buf)[i],'T','N');
      basics::multTranspose(m_W->X()[i],(*m_buf3)[i],m_GLL2G,'N','N');
    }
    if( m_bc == SEMLaplacianEvaluator::PERIODIC ) {
      m_geometry.periodicDssum(m_W->X());
      m_geometry.periodicDssum(m_W->Y());
      m_geometry.invMassP(m_W->X(),m_division[m_rank].elements);
      m_geometry.invMassP(m_W->Y(),m_division[m_rank].elements);
    } else {
      m_geometry.dssum(*m_W,m_rank,m_size);
      m_geometry.invMass(*m_W,m_division[m_rank].elements);
    }
    if( m_bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
      m_geometry.mask(m_W->X(),m_rank,m_size);

#pragma omp parallel for schedule(static)
    for( int i=0;i<max;++i ) {
      basics::multTranspose((*m_buf3)[i],m_W->X()[i],m_GLL2G,'N','T');
      basics::multTranspose((*m_buf)[i],m_GLL2G,(*m_buf3)[i],'N','N');
      basics::multPointwise((*m_buf2)[i],(*m_buf)[i],(*m_J)[m_division[m_rank].elements[i]]);
    }
    res.axpy(m_nu,*m_buf2);
  }
}

int SEMConsistentPressureEvaluator::solve(basics::matrixStack& res, 
    utilities::Evaluator<basics::matrixStack>* pre)
{
  basics::matrixStack* buf   = utilities::g_manager.clone(res);
  basics::matrixStack* buf2  = utilities::g_manager.clone(res);

  /* find gN */
  getRestriction(*m_restricted,res,m_desc[m_rank],m_rank,m_size);
  if( m_nu > 1.e-12 )
    m_E->LLSolve(*m_restricted);
  else
    m_E->LLSolve(*m_restricted1);
  getProlongiation(*buf2,*m_restricted,m_desc[m_rank],m_nu<1.e-12);
  evaluateSimple(*buf,*buf2);
  *buf *= -1;
  *buf += res;
  orthogonalize(*buf,m_desc[m_rank]);

  /* solve ENpN = gN */
  int iter;
  MPIDotter dotter(m_rank,m_size);
  m_deflated = true;
  if( pre )
    iter = utilities::CGSolver::solve(*buf,*this,*pre,dotter,1.e-8);
  else
    iter = utilities::CGSolver::solve(*buf,*this,dotter,1.e-8);

  m_deflated = false;
  /* find g0 */
  evaluateSimple(*buf2,*buf);
  *buf2 *= -1;
  *buf2 += res;
  getRestriction(*m_restricted,*buf2,m_desc[m_rank],m_rank,m_size);

  /* solve E0p0 = g0 */
  if( m_nu > 1.e-12 )
    m_E->LLSolve(*m_restricted);
  else
    m_E->LLSolve(*m_restricted1);

  /* find the resulting p */
  getProlongiation(res,*m_restricted,m_desc[m_rank],m_nu<1.e-12);
  res += *buf;

  utilities::g_manager.unlock(buf);
  utilities::g_manager.unlock(buf2);

  return( iter );
}

SEMTensoredConsistentPressurePreconditioner::
SEMTensoredConsistentPressurePreconditioner(const basics::geometryStack& geometry,
    const basics::Vector& weight,
    const basics::Vector& weightGL,
    const basics::Vector& grid,
    const basics::Vector& gridGL,
    const basics::Matrix& GLL2G,
    const basics::Matrix& D,
    const Real nu, int rank, int size) :
  m_nu(nu), m_grid(grid), m_geometry(geometry), m_rank(rank), m_size(size)
{
  m_desc = geometry.getCoarseGroups(1)[0];
  for( int i=0;i<m_desc.size();++i ) {
    m_SP.push_back(new poissonSolver(m_desc[i].size1*GLL2G.rows(),m_desc[i].size2*GLL2G.rows(),0,poissonSolver::Nonhomogenous,false,-1));

    basics::Vector Lx("Lx",m_desc[i].size1);
    basics::Vector Ly("Ly",m_desc[i].size2);

    /* find element sizes */
    if( m_desc[i].type2 == basics::coarseDescriptor::FLIP_XY ) {
      for( int j=0;j<m_desc[i].size1;++j )
        Lx[j] = geometry[m_desc[i].elements[j]].getAdjustedSize(weight,grid).second;
      for( int j=0;j<m_desc[i].size2;++j )
        Ly[j] = geometry[m_desc[i].elements[m_desc[i].size1-1+
          j*m_desc[i].size1]].getAdjustedSize(weight,grid).first;
    } else {
      for( int j=0;j<m_desc[i].size1;++j )
        Lx[j] = geometry[m_desc[i].elements[j]].getAdjustedSize(weight,grid).first;
      for( int j=0;j<m_desc[i].size2;++j )
        Ly[j] = geometry[m_desc[i].elements[m_desc[i].size1-1+
          j*m_desc[i].size1]].getAdjustedSize(weight,grid).second;
    }

    /* construct B^-1 */
    basics::Matrix BxM = constructInvMass(Lx,weight);
    basics::Matrix ByM = constructInvMass(Ly,weight);

    /* construct D and Bt */
    basics::Matrix DxM  = constructDivergence(weightGL,GLL2G,D,m_desc[i].size1);
    basics::Matrix DyM  = constructDivergence(weightGL,GLL2G,D,m_desc[i].size2);
    basics::Matrix BtxM = constructIntMass(Lx,weightGL,GLL2G);
    basics::Matrix BtyM = constructIntMass(Ly,weightGL,GLL2G);

    /* construct E and consistent mass */
    basics::Matrix tempx = basics::multTranspose(DxM,BxM,'N','N');
    basics::Matrix Ex = basics::multTranspose(tempx,DxM,'N','T');
    tempx = basics::multTranspose(BtxM,BxM,'N','N');
    basics::Matrix CMx = basics::multTranspose(tempx,BtxM,'N','T');
    basics::Matrix tempy = basics::multTranspose(DyM,ByM,'N','N');
    basics::Matrix Ey = basics::multTranspose(tempy,DyM,'N','T');
    tempy = basics::multTranspose(BtyM,ByM,'N','N');
    basics::Matrix CMy = basics::multTranspose(tempy,BtyM,'N','T');

    m_SP.back()->setOperatorX(Ex,CMx);
    m_SP.back()->setOperatorY(Ey,CMy);
    m_SP.back()->setGridX(gridGL);
  }

  m_desc   = geometry.getCoarseGroups(size)[rank];
  m_coarse = geometry.getCoarseBuffer(gridGL.length(),gridGL.length(),m_desc,"coarse pressure",true);


  m_desc   = geometry.getCoarseGroups(size)[rank];
  m_coarse = geometry.getCoarseBuffer(gridGL.length(),gridGL.length(),
      m_desc,"coarse pressure",true);
}

SEMTensoredConsistentPressurePreconditioner::~SEMTensoredConsistentPressurePreconditioner()
{
  for( int i=0;i<m_SP.size();++i )
    delete m_SP[i];
  delete m_coarse;
}

void SEMTensoredConsistentPressurePreconditioner::evaluate(basics::matrixStack& res,
    const basics::matrixStack& p) const
{
  m_geometry.fineToCoarseL2(*m_coarse,p,m_desc);
  int max=m_coarse->size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i )
    m_SP[m_desc[i].size3]->solve((*m_coarse)[i]);
  m_geometry.coarseToFineL2(res,*m_coarse,m_desc);

  SEMConsistentPressureEvaluator::orthogonalize(res,m_desc);
}

basics::Matrix SEMTensoredConsistentPressurePreconditioner::constructInvMass(const basics::Vector& L,
    const basics::Vector& weight,
    SEMLaplacianEvaluator::BC bc)
{
  int N = L.length()*weight.length()-(L.length()-1);
  basics::Matrix B("B",N,N);
  for( int j=0;j<L.length();++j )
    for( int k=0;k<weight.length();++k )
      B[k+j*weight.length()-j][k+j*weight.length()-j] += weight[k]*L[j]/2;
  for( int k=0;k<B.rows();++k )
    B[k][k] = Real(1)/B[k][k];

  if( bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
    return B.submatrix(utilities::Range::colon(1,B.rows()-2),
        utilities::Range::colon(1,B.cols()-2));
  if( bc == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM ||
      bc == SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM )
    return B.submatrix(utilities::Range::colon(1,B.rows()-1),
        utilities::Range::colon(1,B.cols()-1));
  return B;
}

  basics::Matrix 
SEMTensoredConsistentPressurePreconditioner::constructIntMass(const basics::Vector& L,
    const basics::Vector& weightGL,
    const basics::Matrix& GLL2G,
    SEMLaplacianEvaluator::BC bc)
{
  basics::Matrix Bt("temp 1D intmass",GLL2G.rows()*L.length(),(GLL2G.rows()+2)*L.length()-(L.length()-1));

  for( int n=0;n<L.length();++n ) {
    for( int j=0;j<GLL2G.cols();++j )
      for( int k=0;k<GLL2G.rows();++k ) {
        Bt[n*GLL2G.cols()+j-n][n*GLL2G.rows()+k] = weightGL[k]*GLL2G[j][k]*L[n]/2;
      }
  }

  if( bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
    return Bt.submatrix(utilities::Range::colon(0,Bt.rows()-1),
        utilities::Range::colon(1,Bt.cols()-2));
  if( bc == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM ||
      bc == SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM )
    return Bt.submatrix(utilities::Range::colon(0,Bt.rows()-1),
        utilities::Range::colon(1,Bt.cols()-1));

  return Bt;
}

  basics::Matrix 
SEMTensoredConsistentPressurePreconditioner::constructDivergence(const basics::Vector& weightGL,
    const basics::Matrix& GLL2G,
    const basics::Matrix& D,
    int elem,
    SEMLaplacianEvaluator::BC bc)
{
  /* construct D */
  basics::Matrix Dt("temp 1D divergence",GLL2G.rows()*elem,(GLL2G.rows()+2)*elem-(elem-1));

  basics::Matrix D1 = basics::multTranspose(GLL2G,D,'N','N');

  for( int n=0;n<elem;++n ) {
    for( int j=0;j<D1.cols();++j )
      for( int k=0;k<D1.rows();++k )
        Dt[n*D1.cols()+j-n][n*D1.rows()+k] = weightGL[k]*D1[j][k];
  }
  if( bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
    return Dt.submatrix(utilities::Range::colon(0, Dt.rows()-1),
        utilities::Range::colon(1, Dt.cols()-2));
  if( bc == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM  ||
      bc == SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM )
    return Dt.submatrix(utilities::Range::colon(0, Dt.rows()-1),
        utilities::Range::colon(1, Dt.cols()-1));

  return Dt;
}

basics::Matrix SEMTensoredConsistentPressurePreconditioner::constructRestriction(int N, int elem)
{
  basics::Matrix result("restriction",1,N*elem);
  result = 1;

  return( result );
}

basics::Matrix SEMTensoredConsistentPressurePreconditioner::constructE0(const basics::Matrix& E,
    const basics::Matrix& R)
{
  basics::Matrix result("coarse E",R.rows(),R.rows());

  basics::Vector temp("temp",R.rows());

  for( int i=0;i<R.rows();++i ) {
    temp.clear();
    temp[i] = 1;
    basics::Vector temp2 = multTranspose(R,temp,'T');
    temp2 = E*temp2;
    temp = basics::multTranspose(R,temp2,'N');
    result[i] = temp;
  }

  return( result );
}

void SEMTensoredConsistentPressurePreconditioner::constructEn(basics::Matrix& E, int N, int elem)
{
  basics::Matrix R = constructRestriction(N,elem);
  basics::Matrix E0 = constructE0(E,R);
  E0.invert();
  basics::Matrix temp = basics::multTranspose(R,E,'N','N');
  basics::Matrix temp3 = basics::multTranspose(E0,temp,'N','N');
  basics::Matrix bar = basics::multTranspose(R,temp3,'T','N');
  basics::Matrix CM = basics::multTranspose(E,bar,'N','N');
  E -= CM;
}

SEMFEMConsistentPressurePreconditioner::
SEMFEMConsistentPressurePreconditioner(basics::geometryStack& geometry,
    const basics::Vector& weight,
    const basics::Vector& weightGL,
    const basics::Vector& grid,
    const basics::Vector& gridGL,
    const Real nu, int rank, int size) :
  m_geometry(geometry), m_rank(rank), m_size(size),
  m_GLL2G("prolongation",gridGL.length(),grid.length()), m_nu(nu)
{
  m_GLL2G = utilities::GLL::interpolationMatrix(grid,gridGL);
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

  m_SP = new linearElements::poissonSolver(linearElements::poissonSolver::Homogenous,
      gridGL,Lx,Ly,true);

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

  m_SP2 = new linearElements::poissonSolver(linearElements::poissonSolver::Nonhomogenous,
      gridGL,Lx2,Ly2,true);
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
  m_coarse  = m_geometry.getCoarseBuffer(gridGL.length(),gridGL.length(),m_group[m_rank],"coarse",true);
  m_coarse2 = m_geometry.getCoarseBuffer(gridGL.length(),gridGL.length(),m_group[m_rank],"coarse2",true);
  m_coarse3 = m_geometry.getCoarseBuffer(grid.length(),grid.length(),m_group[m_rank],"coarse3",false,true);

  if( m_geometry.hasCoarseSolver() ) {
    m_velTemp = utilities::g_manager.aquireMatrixStack("temp velocity",grid.length(),gridGL.length(),geometry.size());

    basics::coarseDescriptor desc = m_geometry.getRestrictionGridInfo();
    if( desc.size1 > 0 ) {
      m_work = utilities::g_manager.aquireMatrixStack("working buffer",desc.size1,desc.size1,desc.size3);
      basics::Matrix* A = new basics::Matrix("helmholtz",desc.size4,desc.size4);
      basics::Matrix* B = new basics::Matrix("helmholtz",desc.size4,desc.size4);

      HDF5::HDF5Reader reader(m_geometry.getOperatorFile());
      reader.read(*A,"An");
      reader.read(*B,"Bn");
      m_LG = utilities::g_manager.aquireMatrixStack("local to global",desc.size1,desc.size1,desc.size3);
      reader.read(*m_LG,"LGn");
      if( m_nu > 1.e-12 ) {
        m_A = new basics::Matrix(*A);
        m_A->axpy(m_nu,*B);
      } else
        m_A = new basics::Matrix(A->submatrix(utilities::Range::colon(0,desc.size4-2),utilities::Range::colon(0,desc.size4-2)));
      *m_A = m_A->choleskyFactorize();

      delete A;
      delete B;

      m_restricted = utilities::g_manager.aquireVector("restricted",desc.size4);
      if( m_nu < 1.e-12 )
        m_restricted1 = new basics::Vector("restricted minus one",m_restricted->length()-1,false,m_restricted->data());
      else
        m_restricted1 = m_restricted;
      m_restrict  = utilities::g_manager.aquireMatrix("restrict buffer2",(*m_coarse3)[0].rows(),desc.size1);
      m_velBuf  = utilities::g_manager.aquireMatrixStack("velocity buffer",grid.length(),grid.length(),geometry.size());
      m_velBuf2 = utilities::g_manager.aquireMatrixStack("velocity buffer",grid.length(),grid.length(),geometry.size());

      SEMFEMLaplacianPreconditioner::
        setupRestrictionOperators(m_RT,geometry,grid,
            group[ n].size1,group[ n].size2,
            group[n2].size1,group[n2].size2,
            SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN);

    }
  } else
    m_restricted = NULL;

}

SEMFEMConsistentPressurePreconditioner::~SEMFEMConsistentPressurePreconditioner()
{
  utilities::g_manager.unlock(m_restrict);
  utilities::g_manager.unlock(m_restricted);
  utilities::g_manager.unlock(m_work);
  utilities::g_manager.unlock(m_velBuf);
  utilities::g_manager.unlock(m_velBuf2);
  utilities::g_manager.unlock(m_LG);
  utilities::g_manager.unlock(m_velTemp);

  delete m_coarse;
  delete m_coarse2;
  delete m_coarse3;
  delete m_A;
  if( m_nu < 1.e-12 )
    delete m_restricted1;
  delete m_SP;
  delete m_SP2;
}

void SEMFEMConsistentPressurePreconditioner::evaluate(basics::matrixStack& res,
    const basics::matrixStack& p) const
{
  m_geometry.fineToCoarseL2(*m_coarse,p,m_group[m_rank]);
  /* local solves */
  int max=m_group[m_rank].size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    Real Lx = m_Lx[m_group[m_rank][i].size3];
    Real Ly = m_Ly[m_group[m_rank][i].size3];
    if( m_group[m_rank][i].type == basics::coarseDescriptor::FULL_OVERLAP )
      m_SP->solve((*m_coarse)[i],(*m_coarse2)[i],m_nu*Ly*Lx,1,Ly/Lx);
    else
      m_SP2->solve((*m_coarse)[i],(*m_coarse2)[i],m_nu*Ly*Lx,1,Ly/Lx);
  }
  m_geometry.coarseToFineL2(res,*m_coarse,m_group[m_rank]);

  if( m_geometry.hasCoarseSolver() ) { // coarse solve
    max=p.size();
#pragma omp parallel for schedule(static)
    for( int i=0;i<max;++i ) {
      basics::multTranspose((*m_velTemp)[i],m_GLL2G,p[i],'T','N');
      basics::multTranspose((*m_velBuf)[i],(*m_velTemp)[i],m_GLL2G,'N','N');
    }
    m_geometry.fineToCoarseRestriction(*m_coarse3,*m_velBuf2,*m_velBuf,
        m_group[m_rank],true,m_rank,m_size);

    /* restrict */
    m_geometry.getRestriction(*m_restricted,*m_work,*m_LG,*m_coarse3,*m_restrict,
        m_group[m_rank],m_RT);

    /* coarse solve */
    m_A->LLSolve(*m_restricted1);

    if( m_nu < 1.e-12 )
      (*m_restricted)[m_restricted->length()-1] = 0;

    /* prolong */
    m_geometry.getProlongiation(*m_coarse3,*m_work,*m_LG,*m_restricted,*m_restrict,m_RT);
    m_geometry.coarseToFine(*m_velBuf,*m_coarse3,m_group[m_rank],true);
    max=m_velBuf->size();
#pragma omp parallel for schedule(static)
    for( int i=0;i<max;++i ) {
      basics::multTranspose((*m_velTemp)[i],(*m_velBuf)[i],m_GLL2G,'N','T');
      basics::multTranspose(res[i],m_GLL2G,(*m_velTemp)[i],'N','N',mnlRealOne,mnlRealOne);
    }
  }
  if( m_nu < 1.e-10 )
    res -= res.sum()/res.length();
}

extrudedConsistentPressureEvaluator::
extrudedConsistentPressureEvaluator(basics::geometryStack& geometry,
    const basics::Matrix& D,
    const basics::Matrix& GLL2G,
    const basics::Vector& weightGL,
    const basics::Vector& weight,
    const basics::Vector& grid,
    const basics::Vector& gridGL,
    bool preconditioned, int rank, 
    int size, SEMLaplacianEvaluator::BC bc,
    int elems) :
  m_eval(geometry,D,GLL2G,weightGL,weight),
  m_SP((GLL2G.rows()-1)*elems,(GLL2G.rows()-1)*elems,0,legendreLegendreW::poissonSolver::Nonhomogenous,false,-1),
  m_gradient(m_eval.m_gradient,rank,size), m_divergence(m_eval.m_divergence,rank,size),
  m_pre3(NULL), m_rank(rank), m_size(size), m_restricted(NULL), m_restricted1(NULL), 
  m_bc(bc), m_stat(gridGL.length()*elems,size,m_scount,m_sdispl,m_rcount,m_rdispl), m_W(NULL), m_E(NULL),
  m_deflated(false)
{
  /* setup consistent mass matrix and desired eigenvalues */
  basics::Vector L("temp",elems);
  L = geometry.m_Lz/elems;
  basics::Matrix Dz = 
    SEMTensoredConsistentPressurePreconditioner::constructDivergence(weightGL,GLL2G,D,elems);
  basics::Matrix B  = 
    SEMTensoredConsistentPressurePreconditioner::constructInvMass(L,weight,m_bc);
  basics::Matrix B2  = 
    SEMTensoredConsistentPressurePreconditioner::constructInvMass(L,weight);
  basics::Matrix Bz = 
    SEMTensoredConsistentPressurePreconditioner::constructIntMass(L,weightGL,GLL2G,m_bc);

  /* setup 'laplacian' */
  basics::Matrix temp = basics::multTranspose(B2,Dz,'N','T');
  basics::Matrix E = basics::multTranspose(Dz,temp,'N','N');

  /* setup 'mass' */
  basics::Matrix temp2 = basics::multTranspose(B,Bz,'N','T');
  basics::Matrix Bs = basics::multTranspose(Bz,temp2,'N','N');

  m_SP.setOperatorX(E,Bs);
  m_SP.setGridX(gridGL);
  m_SP.setWeightX(weightGL);

  /* now push them evaluators to the vector */
  m_evals.push_back(&m_eval);
  int max = 1;
#ifdef OPENMP
  max = omp_get_max_threads();
#endif
  SEMLaplacianEvaluator::BC bc2 = SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET;
  if( m_bc == SEMLaplacianEvaluator::PERIODIC_DIRICHLET_BOTTOM )
    bc2 = SEMLaplacianEvaluator::PERIODIC;
  m_eval.m_bc = bc2;
  for( int i=1;i<max;++i )
    m_evals.push_back(new SEMConsistentPressureEvaluator(geometry,D,GLL2G,weightGL,weight,NULL,0,bc2));
  m_stat.getDisplacements2(m_scount,m_sdispl,m_rcount,m_rdispl,m_rank,m_size,E.length());
  m_counts = m_stat.getCounts2(m_size);
  if( preconditioned ) {
    for( int l=0;l<m_SP.eigX().length();++l ) {
      if( m_SP.eigX()[l] > 1000 )
        m_pre.push_back(new SEMInverseMassEvaluator(geometry,
              weightGL,
              GLL2G));
      else
        m_pre.push_back(new SEMFEMConsistentPressurePreconditioner(geometry,
              weight,
              weightGL,
              grid,
              gridGL,
              m_SP.eigX()[l]));
    }
    //        m_pre3 = new extrudedTensoredConsistentPressurePreconditioner(geometry,
    //                                                                        weight,
    //                                                                        weightGL, 
    //                                                                        grid,
    //                                                                        gridGL,
    //                                                                        GLL2G, D);
  }
  m_desc     = geometry.getCoarseGroups(m_size);
  m_division = geometry.getDivisionInfo(m_size);
  m_W = utilities::g_manager.aquireMatricesStackField("temp field X",gridGL.length()+2,
      gridGL.length()+2,gridGL.length()+2,
      m_division[m_rank].elements.size());

  if( 0 && size == 1 && geometry.hasCoarseSolver() ) {
    basics::coarseDescriptor desc = geometry.getRestrictionGridInfo();
    m_buf1 	= utilities::g_manager.aquireMatricesStack("buffer",weightGL.length(),weightGL.length(),weightGL.length(),m_division[m_rank].elements.size());
    m_buf2 	= utilities::g_manager.clone(*m_buf1);
    if( desc.size3 > 0 ) {
      m_restricted  = utilities::g_manager.aquireVector("restricted",desc.size3);
      m_restricted1 = new basics::Vector("restricted minus one",m_restricted->length()-1,false,m_restricted->data());
      m_E = new basics::Matrix("coarse E",desc.size3-1,desc.size3-1);

      for( int i=0;i<m_E->cols();++i ) {
        m_restricted->clear();
        (*m_restricted)[i] = 1;
        SEMConsistentPressureEvaluator::getProlongiation(*m_buf1,*m_restricted,m_desc[m_rank]);
        evaluateSimple(*m_buf2,*m_buf1);
        SEMConsistentPressureEvaluator::getRestriction(*m_restricted,*m_buf2,m_desc[m_rank],m_rank,m_size);
        (*m_E)[i] = *m_restricted1;
      }

      *m_E = m_E->choleskyFactorize();
    }
  }
}

extrudedConsistentPressureEvaluator::~extrudedConsistentPressureEvaluator()
{
  for( int i=0;i<m_pre.size();++i )
    delete m_pre[i];
  for( int i=1;i<m_evals.size();++i )
    delete m_evals[i];
  utilities::g_manager.unlock(m_W);
  delete m_E;
  delete m_pre3;
  utilities::g_manager.unlock(m_restricted);
  delete m_restricted1;
  m_stat.cleanDisplacements(m_scount,m_sdispl,m_rcount,m_rdispl);
}

void extrudedConsistentPressureEvaluator::evaluate(basics::matricesStack& res,
    const basics::matricesStack& p) const
{
  evaluateSimple(res,p);
  if( m_deflated ) {
    SEMConsistentPressureEvaluator::getRestriction(*m_restricted,res,m_desc[m_eval.m_rank],m_eval.m_rank,m_eval.m_size);
    m_E->LLSolve(*m_restricted1);
    SEMConsistentPressureEvaluator::getProlongiation(*m_buf1,*m_restricted,m_desc[m_eval.m_rank]);
    evaluateSimple(*m_buf2,*m_buf1);
    res -= *m_buf2;
    SEMConsistentPressureEvaluator::orthogonalize(res,m_desc[0]);
  }
  applyFilter(res);
}

void extrudedConsistentPressureEvaluator::evaluateSimple(basics::matricesStack& res, 
    const basics::matricesStack& p) const
{
  m_gradient.evaluate(*m_W,p);
  dssum(*m_W,m_eval.m_geometry,m_rank,m_size);
  m_eval.m_geometry.invMass(*m_W,m_division[m_rank].elements);
  mask(*m_W,m_eval.m_geometry,m_rank,m_size,m_bc);
  m_divergence.evaluate(res,*m_W);
}

int extrudedConsistentPressureEvaluator::solve3(basics::matricesStack& res,
    basics::matricesStack& p,
    extrudedTensoredConsistentPressurePreconditioner* pre)
{
  basics::matricesStack* buf   = utilities::g_manager.clone(res);
  basics::matricesStack* buf2  = utilities::g_manager.clone(res);

  /* find gN */
  SEMConsistentPressureEvaluator::getRestriction(*m_restricted,res,m_desc[m_rank],m_rank,m_size);
  m_E->LLSolve(*m_restricted1);
  SEMConsistentPressureEvaluator::getProlongiation(*buf2,*m_restricted,m_desc[m_rank]);
  evaluateSimple(*buf,*buf2);
  *buf *= -1;
  *buf += res;
  SEMConsistentPressureEvaluator::orthogonalize(*buf,m_desc[m_rank]);

  /* solve ENpN = gN */
  int iter;
  MPIDotter dotter(m_rank,m_size);
  m_deflated = true;
  if( pre )
    cout << "iterations pressure " << (iter = utilities::CGSolver::solve(*buf,*this,*pre,dotter,1.e-10)) << endl;
  else
    cout << "iterations pressure " << (iter = utilities::CGSolver::solve(*buf,*this,dotter,1.e-10)) << endl;
  m_deflated = false;

  /* find g0 */
  evaluateSimple(*buf2,*buf);
  *buf2 *= -1;
  *buf2 += res;
  SEMConsistentPressureEvaluator::getRestriction(*m_restricted,*buf2,m_desc[m_rank],m_rank,m_size);

  /* solve E0p0 = g0 */
  m_E->LLSolve(*m_restricted1);

  /* find the resulting p */
  SEMConsistentPressureEvaluator::getProlongiation(res,*m_restricted,m_desc[m_rank]);
  res += *buf;

  utilities::g_manager.unlock(buf);
  utilities::g_manager.unlock(buf2);

  return( iter );
}

void extrudedConsistentPressureEvaluator::solve(basics::matricesStack& res)
{
  basics::matricesStack* p = utilities::g_manager.clone(res);
  /* pretransform */
  int max=p->size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i )
    applyLocalGlobal((*p)[i],res[i],*m_SP.Qx(),'N','N',0);

  vector<basics::matrixStack*> stack;
  m_counts = m_stat.getCounts2(m_size);
#ifdef HAS_MPI
  m_stat.getDisplacements2(m_scount,m_sdispl,m_rcount,m_rdispl,m_rank,m_size,res[0][0].length());
  basics::matricesStack* temp = 
    utilities::g_manager.aquireMatricesStack("temp",res[0].rows(),res[0].cols(),
        (m_rdispl[0][m_size-1]+m_rcount[0][m_size-1])/(res[0][0].length()),res.size());
  stack = sendAndSetupStack(*temp,*p,SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN,m_rank,m_size,
      m_scount[0],m_sdispl[0],m_rcount[0],m_rdispl[0],&m_counts.second);
#else
  for( int l=0;l<(*p)[0].matrices();++l ) 
    stack.push_back(&p->at(m_counts.second[m_rank][l]));
#endif

  Real tol = 1.e-10;
  max=stack.size();
  //    vector<int> start = m_stat.getStarts(m_size);
  m_stat.reset();
#pragma omp parallel for schedule(dynamic)
  for( int l=0;l<max;++l ) {
    int iter;
#ifdef OPENMP
    int eval = omp_get_thread_num();
#else
    int eval = 0;
#endif
    //        int plane = start[m_rank]+l;
    int plane = m_counts.second[m_rank][l];
    m_evals[eval]->m_nu = m_SP.eigX()[plane];
    if( m_pre.size() ) {
      iter = utilities::CGSolver::solve(*stack[l],*m_evals[eval],*m_pre[plane],stack[l]->getDotter(),tol);
    } else
      iter = utilities::CGSolver::solve(*stack[l],*m_evals[eval],stack[l]->getDotter(),tol);

    m_stat.get()[plane] = iter;
    cout << "iterations pressure (" << plane << ") " << iter << endl;
  }
#ifdef HAS_MPI
  sendStack(*p,*temp,m_scount[1],m_sdispl[1],m_rcount[1],m_rdispl[1],
      SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN,&stack,&m_counts.second);
  m_stat.exchange();
  utilities::g_manager.unlock(temp);
#endif
  max=p->size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i )
    applyLocalGlobal(res[i],(*p)[i],*m_SP.Qx(),'N','T');

  utilities::g_manager.unlock(p);
}

extrudedTensoredConsistentPressurePreconditioner::
extrudedTensoredConsistentPressurePreconditioner(const basics::geometryStack& geometry,
    const basics::Vector& weight,
    const basics::Vector& weightGL,
    const basics::Vector& grid,
    const basics::Vector& gridGL,
    const basics::Matrix& GLL2G,
    const basics::Matrix& D,
    int rank, int size) :
  m_grid(grid), m_geometry(geometry), m_rank(rank), m_size(size)
{
  m_desc = geometry.getCoarseGroups(1)[0];

  for( int i=0;i<m_desc.size();++i ) {
    m_SP.push_back(new poissonSolver(m_desc[i].size1*GLL2G.rows(),m_desc[i].size2*GLL2G.rows(),GLL2G.rows(),poissonSolver::Nonhomogenous,false,-1,-1,-1));

    basics::Vector Lx("Lx",m_desc[i].size1);
    basics::Vector Ly("Ly",m_desc[i].size2);

    /* find element sizes */
    if( m_desc[i].type2 == basics::coarseDescriptor::FLIP_XY ) {
      for( int j=0;j<m_desc[i].size1;++j )
        Lx[j] = geometry[m_desc[i].elements[j]].getAdjustedSize(weight,grid).second;
      for( int j=0;j<m_desc[i].size2;++j )
        Ly[j] = geometry[m_desc[i].elements[m_desc[i].size1-1+
          j*m_desc[i].size1]].getAdjustedSize(weight,grid).first;
    } else {
      for( int j=0;j<m_desc[i].size1;++j )
        Lx[j] = geometry[m_desc[i].elements[j]].getAdjustedSize(weight,grid).first;
      for( int j=0;j<m_desc[i].size2;++j )
        Ly[j] = geometry[m_desc[i].elements[m_desc[i].size1-1+
          j*m_desc[i].size1]].getAdjustedSize(weight,grid).second;
    }

    /* construct B^-1 */
    basics::Vector Lz("Lz",1);
    Lz = 2;
    basics::Matrix BxM = SEMTensoredConsistentPressurePreconditioner::constructInvMass(Lx,weight);
    basics::Matrix ByM = SEMTensoredConsistentPressurePreconditioner::constructInvMass(Ly,weight);
    basics::Matrix BzM = SEMTensoredConsistentPressurePreconditioner::constructInvMass(Lz,weight);

    /* construct D and Bt */
    basics::Matrix DxM = SEMTensoredConsistentPressurePreconditioner::constructDivergence(weightGL,GLL2G,D,m_desc[i].size1);
    basics::Matrix DyM = SEMTensoredConsistentPressurePreconditioner::constructDivergence(weightGL,GLL2G,D,m_desc[i].size2);
    basics::Matrix DzM = SEMTensoredConsistentPressurePreconditioner::constructDivergence(weightGL,GLL2G,D,1);

    basics::Matrix BtxM = SEMTensoredConsistentPressurePreconditioner::constructIntMass(Lx,weightGL,GLL2G);
    basics::Matrix BtyM = SEMTensoredConsistentPressurePreconditioner::constructIntMass(Ly,weightGL,GLL2G);
    basics::Matrix BtzM = SEMTensoredConsistentPressurePreconditioner::constructIntMass(Lz,weightGL,GLL2G);

    /* construct E and consistent mass */
    basics::Matrix tempx = basics::multTranspose(DxM,BxM,'N','N');
    basics::Matrix Ex = basics::multTranspose(tempx,DxM,'N','T');
    tempx = basics::multTranspose(BtxM,BxM,'N','N');
    basics::Matrix CMx = basics::multTranspose(tempx,BtxM,'N','T');
    m_SP.back()->setOperatorX(Ex,CMx);

    basics::Matrix tempy = basics::multTranspose(DyM,ByM,'N','N');
    basics::Matrix Ey = basics::multTranspose(tempy,DyM,'N','T');
    tempy = basics::multTranspose(BtyM,ByM,'N','N');
    basics::Matrix CMy = basics::multTranspose(tempy,BtyM,'N','T');
    m_SP.back()->setOperatorY(Ey,CMy);

    basics::Matrix tempz = basics::multTranspose(DzM,BzM,'N','N');
    basics::Matrix Ez = basics::multTranspose(tempz,DzM,'N','T');
    tempz = basics::multTranspose(BtzM,BzM,'N','N');
    basics::Matrix CMz = basics::multTranspose(tempz,BtzM,'N','T');
    m_SP.back()->setOperatorZ(Ez,CMz);

    m_SP.back()->setGridX(gridGL);
  }
  m_desc   = geometry.getCoarseGroups(size)[rank];
  m_coarse = geometry.getCoarseBuffer(gridGL.length(),gridGL.length(),gridGL.length(),m_desc,"coarse pressure",true);
}

extrudedTensoredConsistentPressurePreconditioner::~extrudedTensoredConsistentPressurePreconditioner()
{
  for( int i=0;i<m_SP.size();++i )
    delete m_SP[i];
  if ( m_coarse )
    delete m_coarse;
}

void extrudedTensoredConsistentPressurePreconditioner::evaluate(basics::matricesStack& res,
    const basics::matricesStack& p) const
{
  m_geometry.fineToCoarseL2(*m_coarse,p,m_desc);
  for( int i=0;i<m_coarse->size();++i )
    m_SP[m_desc[i].size3]->solve((*m_coarse)[i]);

  m_geometry.coarseToFineL2(res,*m_coarse,m_desc);
  SEMConsistentPressureEvaluator::orthogonalize(res,m_desc);
}

extrudedFEMConsistentPressurePreconditioner::
extrudedFEMConsistentPressurePreconditioner(basics::geometryStack& geometry,
    const basics::Vector& weight,
    const basics::Vector& weightGL,
    const basics::Vector& grid,
    const basics::Vector& gridGL,
    int rank, int size) :
  m_geometry(geometry), m_rank(rank), m_size(size),
  m_GLL2G("prolongation",gridGL.length(),grid.length()),
  m_velBuf(NULL), m_velBuf2(NULL), m_velTemp(NULL), m_velTemp2(NULL),
  m_A(NULL), m_work(NULL), m_work2(NULL), m_LG(NULL),
  m_restricted(NULL), m_restricted1(NULL), m_restrict(NULL)
{
  m_GLL2G = utilities::GLL::interpolationMatrix(grid,gridGL);
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

  m_SP = new linearElements::poissonSolver(linearElements::poissonSolver::Homogenous,
      gridGL,Lx,Ly,true,true);

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

  m_SP2 = new linearElements::poissonSolver(linearElements::poissonSolver::Nonhomogenous,
      gridGL,Lx2,Ly2,true,true);
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
  m_coarse  = m_geometry.getCoarseBuffer(gridGL.length(),gridGL.length(),gridGL.length(),m_group[m_rank],"coarse",true);
  m_coarse2 = m_geometry.getCoarseBuffer(gridGL.length(),gridGL.length(),gridGL.length(),m_group[m_rank],"coarse2",true);
  m_coarse3 = m_geometry.getCoarseBuffer(grid.length(),grid.length(),grid.length(),m_group[m_rank],"coarse3",false,true);

  if( m_geometry.hasCoarseSolver() ) {
    m_velTemp  = new basics::Matrix("temp velocity",grid.length(),gridGL.length());
    m_velTemp2 = new basics::Matrices("temp velocity",grid.length(),grid.length(),gridGL.length());

    basics::coarseDescriptor desc = m_geometry.getRestrictionGridInfo();
    if( desc.size1 > 0 ) {
      m_work  = utilities::g_manager.aquireMatricesStack("working buffer",desc.size1,desc.size1,desc.size1,desc.size3);
      m_work2 = m_geometry.getCoarseBuffer(grid.length(),grid.length(),desc.size1,m_group[m_rank],"working buffer 2",false,true);
      m_work3 = m_geometry.getCoarseBuffer(grid.length(),grid.length(),desc.size1,group,"working buffer 2",false,true);
      basics::Matrix* A = new basics::Matrix("helmholtz",desc.size4*desc.size1,desc.size4*desc.size1);

      m_LG = utilities::g_manager.aquireMatricesStack("local to global",desc.size1,desc.size1,desc.size1,desc.size3);
      HDF5::HDF5Reader reader(m_geometry.getOperatorFile());
      if( reader.read(*A,"A3n") ) {
        reader.read(*m_LG,"LG3n");
        m_A = new basics::Matrix(A->submatrix(utilities::Range::colon(0,A->rows()-2),utilities::Range::colon(0,A->cols()-2)));
        *m_A = m_A->choleskyFactorize();
      }

      delete A;

      m_restricted = utilities::g_manager.aquireVector("restricted",desc.size4*desc.size1);
      m_restricted1 = new basics::Vector("restricted minus one",m_restricted->length()-1,false,m_restricted->data());
      m_restrict  = utilities::g_manager.aquireMatrix("restrict buffer2",(*m_coarse3)[0].rows(),desc.size1);
      m_velBuf  = utilities::g_manager.aquireMatricesStack("velocity buffer",grid.length(),grid.length(),grid.length(),geometry.size());
      m_velBuf2 = utilities::g_manager.aquireMatricesStack("velocity buffer",grid.length(),grid.length(),grid.length(),geometry.size());

      SEMFEMLaplacianPreconditioner::
        setupRestrictionOperators(m_RT,geometry,grid,
            group[ n].size1,group[ n].size2,
            group[n2].size1,group[n2].size2,
            SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN);
    }
  } else
    m_restricted = NULL;
}

extrudedFEMConsistentPressurePreconditioner::~extrudedFEMConsistentPressurePreconditioner()
{
  utilities::g_manager.unlock(m_restrict);
  utilities::g_manager.unlock(m_restricted);
  utilities::g_manager.unlock(m_velBuf);
  utilities::g_manager.unlock(m_velBuf2);
  utilities::g_manager.unlock(m_LG);
  utilities::g_manager.unlock(m_work);

  delete m_coarse;
  delete m_coarse2;
  delete m_coarse3;
  delete m_A;
  delete m_velTemp;
  delete m_work2;
  delete m_restricted1;
}

void extrudedFEMConsistentPressurePreconditioner::evaluate(basics::matricesStack& res,
    const basics::matricesStack& p) const
{
  m_geometry.fineToCoarseL2(*m_coarse,p,m_group[m_rank]);
  /* local solves */
  int max=m_group[m_rank].size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    if( m_group[m_rank][i].type == basics::coarseDescriptor::FULL_OVERLAP ) {
      m_SP->solve((*m_coarse)[i],(*m_coarse2)[i],
          m_Lx[m_group[m_rank][i].size3],
          m_Ly[m_group[m_rank][i].size3],
          m_geometry.m_Lz);
    }
    else
      m_SP2->solve((*m_coarse)[i],(*m_coarse2)[i],
          m_Lx[m_group[m_rank][i].size3],
          m_Ly[m_group[m_rank][i].size3],
          m_geometry.m_Lz);
  }
  m_geometry.coarseToFineL2(res,*m_coarse,m_group[m_rank]);

  if( m_geometry.hasCoarseSolver() ) { // coarse solve
    for( int n=0;n<p.size();++n ) {
      for( int l=0;l<p[n].matrices();++l ) {
        basics::multTranspose(*m_velTemp,m_GLL2G,p[n][l],'T','N');
        basics::multTranspose((*m_velTemp2)[l],*m_velTemp,m_GLL2G,'N','N');
      }
      basics::applyLocalGlobal((*m_velBuf)[n],*m_velTemp2,m_GLL2G,'N','N');
    }
    m_geometry.fineToCoarseRestriction(*m_coarse3,*m_velBuf2,*m_velBuf,
        m_group[m_rank],true,m_rank,m_size);

    /* restrict */
    m_geometry.gatherAndRestrict(*m_restricted,*m_work,*m_work2,
        *m_work3,*m_LG,*m_coarse3,*m_restrict,m_RT,
        m_rank,m_size);
    /* coarse solve */
    if( m_rank == 0 ) {
      m_A->LLSolve(*m_restricted1);
      (*m_restricted)[m_restricted->length()-1] = 0;
    }

    /* prolong */
    m_geometry.prolongAndScatter(*m_coarse3,*m_work,*m_work2,*m_work3,
        *m_LG,*m_restricted,*m_restrict,m_RT,
        m_rank,m_size);
    m_geometry.coarseToFine(*m_velBuf,*m_coarse3,m_group[m_rank],true);
    max=res.size();
    for( int n=0;n<max;++n ) {
      basics::applyLocalGlobal(*m_velTemp2,(*m_velBuf)[n],m_GLL2G,'N','T');
      for( int l=0;l<m_velTemp2->matrices();++l ) {
        basics::multTranspose(*m_velTemp,(*m_velTemp2)[l],m_GLL2G,'N','T');
        basics::multTranspose(res[n][l],m_GLL2G,*m_velTemp,'N','N',mnlRealOne,mnlRealOne);
      }
    }
  }

  extrudedConsistentPressureEvaluator::applyFilter(res);
  //    res -= res.sum()/res.length();
}

SEM3DConsistentPressureEvaluator::SEM3DConsistentPressureEvaluator(basics::geometryStack3D& geometry,
    const basics::Matrix& D,
    const basics::Matrix& GLL2G,
    const basics::Vector& weightGL,
    const basics::Vector& weight,
    int rank, int size,
    SEMLaplacianEvaluator::BC bc) :
  m_divergence(geometry,D,GLL2G,weightGL,rank,size),
  m_gradient(geometry,D,GLL2G,weightGL,NULL,rank,size), 
  m_rank(rank), m_size(size), m_bc(bc)
{
  m_division = geometry.getDivisionInfo(size);
  m_W = utilities::g_manager.aquireMatricesStackField("temp field X",weightGL.length()+2,
      weightGL.length()+2,weightGL.length()+2,
      m_division[m_rank].elements.size());
}

SEM3DConsistentPressureEvaluator::~SEM3DConsistentPressureEvaluator()
{
  utilities::g_manager.unlock(m_W);
}

void SEM3DConsistentPressureEvaluator::evaluate(basics::matricesStack& res,
    const basics::matricesStack& p) const
{
  m_gradient.evaluate(*m_W,p);
  dssum(*m_W,m_divergence.m_geometry,m_rank,m_size);
  m_divergence.m_geometry.invMass(*m_W,m_division[m_rank].elements);
  mask(*m_W,m_divergence.m_geometry,m_rank,m_size,m_bc);
  m_divergence.evaluate(res,*m_W);
  Real sum = res.sum();
  int length=res.length();
#ifdef HAS_MPI
  Real mysum = sum;
  MPI_Allreduce(&mysum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  length *= m_size;
#endif
  res -= sum/length;
}

SEM3DFEMConsistentPressurePreconditioner::
SEM3DFEMConsistentPressurePreconditioner(basics::geometryStack3D& geometry,
    const basics::Vector& weight,
    const basics::Vector& weightGL,
    const basics::Vector& grid,
    const basics::Vector& gridGL,
    int rank, int size) :
  m_geometry(geometry), m_rank(rank), m_size(size),
  m_GLL2G("prolongation",gridGL.length(),grid.length())
{
  m_GLL2G = utilities::GLL::interpolationMatrix(grid,gridGL);
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

  m_SP.push_back(new linearElements::poissonSolver(linearElements::poissonSolver::Homogenous,
        gridGL,Lx,Ly,true,true));

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

  m_SP.push_back(new linearElements::poissonSolver(linearElements::poissonSolver::Nonhomogenous,
        gridGL,Lx2,Ly2,true,true));
  for( int c=0;c<group.size();++c ) {
    Real lx=0;
    Real ly=0;
    Real volume=0;

    if( group[c].type2 == basics::coarseDescriptor::FLIP_XY ) {
      for( int j=0;j<group[c].size1;++j )
        lx += m_geometry[group[c].elements[j]].getSize(weight,grid)[1];
      for( int j=0;j<group[c].size2;++j )
        ly += m_geometry[group[c].elements[group[c].size1-1+
          j*group[c].size1]].getSize(weight,grid)[0];
    } else {
      for( int j=0;j<group[c].size1;++j )
        lx += m_geometry[group[c].elements[j]].getSize(weight,grid)[0];
      for( int j=0;j<group[c].size2;++j )
        ly += m_geometry[group[c].elements[group[c].size1-1+
          j*group[c].size1]].getSize(weight,grid)[1];
    }
    Real lz=0;
    for( int i=0;i<group[c].elements.size();++i ) {
      volume += m_geometry[group[c].elements[i]].getVolume(weight);
      lz += m_geometry[group[c].elements[i]].getSize(weight,grid)[2];
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
  m_coarse  = m_geometry.getCoarseBuffer(gridGL.length(),gridGL.length(),gridGL.length(),m_group[m_rank],"coarse",true);
  m_coarse2 = m_geometry.getCoarseBuffer(gridGL.length(),gridGL.length(),gridGL.length(),m_group[m_rank],"coarse2",true);
  m_coarse3 = m_geometry.getCoarseBuffer(grid.length(),grid.length(),grid.length(),m_group[m_rank],"coarse3");

  if( m_geometry.hasCoarseSolver() ) {
    m_velTemp  = new basics::Matrix("temp velocity",grid.length(),gridGL.length());
    m_velTemp2 = new basics::Matrices("temp velocity",grid.length(),grid.length(),gridGL.length());

    basics::coarseDescriptor desc = m_geometry.getRestrictionGridInfo();
    if( desc.size1 > 0 ) {
      m_work  = utilities::g_manager.aquireMatricesStack("working buffer",desc.size1,desc.size1,desc.size1,desc.size3);
      m_work2 = m_geometry.getCoarseBuffer(grid.length(),grid.length(),desc.size1,m_group[m_rank],"working buffer 2");
      basics::Matrix* A = new basics::Matrix("helmholtz",desc.size4,desc.size4);

      m_LG = utilities::g_manager.aquireMatricesStack("local to global",desc.size1,desc.size1,desc.size1,desc.size3);
      HDF5::HDF5Reader reader(m_geometry.getOperatorFile());
      reader.read(*m_LG,"LG3Dn");
      reader.read(*A,"A3n");
      m_A = new basics::Matrix(A->submatrix(utilities::Range::colon(0,A->rows()-2),utilities::Range::colon(0,A->cols()-2)));
      *m_A = m_A->choleskyFactorize();

      delete A;

      m_restricted = utilities::g_manager.aquireVector("restricted",desc.size4*desc.size1);
      m_restricted1 = new basics::Vector("restricted minus one",m_restricted->length()-1,false,m_restricted->data());
      m_restrict  = utilities::g_manager.aquireMatrix("restrict buffer2",(*m_coarse3)[0].rows(),desc.size1);
      m_velBuf  = utilities::g_manager.aquireMatricesStack("velocity buffer",grid.length(),grid.length(),grid.length(),geometry.size());
      m_velBuf2 = utilities::g_manager.aquireMatricesStack("velocity buffer",grid.length(),grid.length(),grid.length(),geometry.size());
    }
  } else
    m_restricted = NULL;
}

SEM3DFEMConsistentPressurePreconditioner::~SEM3DFEMConsistentPressurePreconditioner()
{
  for( int i=0;i<m_SP.size();++i )
    delete m_SP[i];

  utilities::g_manager.unlock(m_restrict);
  utilities::g_manager.unlock(m_restricted);
  utilities::g_manager.unlock(m_velBuf);
  utilities::g_manager.unlock(m_velBuf2);
  utilities::g_manager.unlock(m_LG);
  utilities::g_manager.unlock(m_work);

  delete m_coarse;
  delete m_coarse2;
  delete m_coarse3;
  delete m_A;
  delete m_velTemp;
  delete m_velTemp2;
  delete m_work2;
  delete m_restricted1;
}

void SEM3DFEMConsistentPressurePreconditioner::evaluate(basics::matricesStack& res,
    const basics::matricesStack& p) const
{
  m_geometry.fineToCoarseL2(*m_coarse,p,m_group[m_rank]);
  /* local solves */
  int max=m_group[m_rank].size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    if( m_group[m_rank][i].type == basics::coarseDescriptor::FULL_OVERLAP ) {
      m_SP[0]->solve((*m_coarse)[i],(*m_coarse2)[i],
          m_Lx[m_group[m_rank][i].size3],
          m_Ly[m_group[m_rank][i].size3],
          m_Lz[m_group[m_rank][i].size3]);
    }
    else
      m_SP[1]->solve((*m_coarse)[i],(*m_coarse2)[i],
          m_Lx[m_group[m_rank][i].size3],
          m_Ly[m_group[m_rank][i].size3],
          m_Lz[m_group[m_rank][i].size3]);
  }
  m_geometry.coarseToFineL2(res,*m_coarse,m_group[m_rank]);

  if( m_geometry.hasCoarseSolver() ) { // coarse solve
    for( int n=0;n<p.size();++n ) {
      for( int l=0;l<p[n].matrices();++l ) {
        basics::multTranspose(*m_velTemp,m_GLL2G,p[n][l],'T','N');
        basics::multTranspose((*m_velTemp2)[l],*m_velTemp,m_GLL2G,'N','N');
      }
      basics::applyLocalGlobal((*m_velBuf)[n],*m_velTemp2,m_GLL2G,'N','N');
    }
    m_geometry.fineToCoarseRestriction(*m_coarse3,*m_velBuf2,*m_velBuf,m_group[m_rank],m_rank,m_size);

    /* restrict */
    m_geometry.getRestriction(*m_restricted,*m_work,*m_work2,*m_LG,*m_coarse3,*m_restrict,m_group[m_rank]);

    /* coarse solve */
    m_A->LLSolve(*m_restricted1);
    (*m_restricted)[m_restricted->length()-1] = 0;

    /* prolong */
    m_geometry.getProlongiation(*m_coarse3,*m_work,*m_work2,*m_LG,*m_restricted,*m_restrict);
    m_geometry.coarseToFine(*m_velBuf,*m_coarse3,m_group[m_rank]);
    max=m_velBuf->size();
    for( int n=0;n<max;++n ) {
      basics::applyLocalGlobal(*m_velTemp2,(*m_velBuf)[n],m_GLL2G,'N','T');
      for( int l=0;l<m_velTemp2->matrices();++l ) {
        basics::multTranspose(*m_velTemp,(*m_velTemp2)[l],m_GLL2G,'N','T');
        basics::multTranspose(res[n][l],m_GLL2G,*m_velTemp,'N','N',mnlRealOne,mnlRealOne);
      }
    }
  }

  res -= res.sum()/res.length();
}
