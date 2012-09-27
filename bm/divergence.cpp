#include "divergence.h"
#include "sem.h"

using namespace mnl;
using namespace std;

deformedDivergenceEvaluator::deformedDivergenceEvaluator(const basics::matrixStack& G,
    const basics::Matrix& Dt,
    const basics::Matrix& GLL2G,
    const basics::Vector& weight,
    basics::Matrix& W,
    basics::Matrix& Wt,
    basics::Matrix& buf) :
  m_G(G), m_weight(weight), m_GLL2G(GLL2G), m_Dt(Dt),
  m_W(W), m_Wt(Wt), m_buf(buf)
{
}

void deformedDivergenceEvaluator::evaluate(basics::Matrix& res,
    const basics::Field2<basics::Matrix>& u, bool mass) const
{
  basics::multTranspose(m_W,m_Dt,u.X(),'N','N');
  basics::multTranspose(m_Wt,m_W,m_GLL2G,'N','T');
  basics::multPointwise(res,m_Wt,m_G[3]);

  basics::multTranspose(m_W,m_GLL2G,u.X(),'N','N');
  basics::multTranspose(m_Wt,m_W,m_Dt,'N','T');
  basics::multPointwise(m_buf,m_Wt,m_G[2]);
  res -= m_buf;

  basics::multTranspose(m_W,m_Dt,u.Y(),'N','N');
  basics::multTranspose(m_Wt,m_W,m_GLL2G,'N','T');
  basics::multPointwise(m_buf,m_Wt,m_G[1]);
  res -= m_buf;

  basics::multTranspose(m_W,m_GLL2G,u.Y(),'N','N');
  basics::multTranspose(m_Wt,m_W,m_Dt,'N','T');
  basics::multPointwise(m_buf,m_Wt,m_G[0]);
  res += m_buf;

  if( mass )
    massReference(res,m_weight);
}

deformed3DDivergenceEvaluator::deformed3DDivergenceEvaluator(const basics::matricesStack& G,
    const basics::Matrix& Dt,
    const basics::Matrix& GLL2G,
    const basics::Vector& weight,
    basics::Matrix& W,
    basics::Matrices& Wt,
    basics::Matrices& buf) :
  m_G(G), m_weight(weight), m_GLL2G(GLL2G), m_Dt(Dt),
  m_W(W), m_Wt(Wt), m_buf(buf)
{
}

deformed3DDivergenceEvaluator::~deformed3DDivergenceEvaluator()
{
}

void deformed3DDivergenceEvaluator::evaluate(basics::Matrices& res,
    const basics::Field3<basics::Matrices>& u) const
{
  evaluateComponent(res,u.X(),0);
  evaluateComponent(res,u.Y(),1);
  evaluateComponent(res,u.Z(),2);
  massReference(res,m_weight);
}

void deformed3DDivergenceEvaluator::evaluateComponent(basics::Matrices& res,
    const basics::Matrices& field,
    int index) const
{
  for( int l=0;l<field.matrices();++l ) {
    basics::multTranspose(m_W,m_Dt,field[l],'N','N');
    basics::multTranspose(m_Wt[l],m_W,m_GLL2G,'N','T');
  }
  basics::applyLocalGlobal(m_buf,m_Wt,m_GLL2G,'N','T');
  basics::multPointwise(m_buf,m_G[index+0]);
  if( index )
    res += m_buf;
  else
    res = m_buf;

  for( int l=0;l<field.matrices();++l ) {
    basics::multTranspose(m_W,m_GLL2G,field[l],'N','N');
    basics::multTranspose(m_Wt[l],m_W,m_Dt,'N','T');
  }
  basics::applyLocalGlobal(m_buf,m_Wt,m_GLL2G,'N','T');
  basics::multPointwise(m_buf,m_G[index+3]);
  res += m_buf;

  for( int l=0;l<field.matrices();++l ) {
    basics::multTranspose(m_W,m_GLL2G,field[l],'N','N');
    basics::multTranspose(m_Wt[l],m_W,m_GLL2G,'N','T');
  }
  basics::applyLocalGlobal(m_buf,m_Wt,m_Dt,'N','T');
  basics::multPointwise(m_buf,m_G[index+6]);
  res += m_buf;
}


SEMDivergenceEvaluator::SEMDivergenceEvaluator(const basics::geometryStack& geometry,
    const basics::Matrix& D,
    const basics::Matrix& GLL2G,
    const basics::Vector& weight,
    const std::vector<basics::matrixStack*>* G,
    int rank, int size) :
  m_geometry(geometry), m_weight(weight), m_GLL2G(GLL2G),
  m_Dt("interpolate derivative matrix", D.rows()-2,D.cols()),
  m_rank(rank), m_size(size)
{
  basics::multTranspose(m_Dt,m_GLL2G,D,'N','N');
  if( G ) {
    m_G = *G;
    m_mine = false;
  }
  else  {
    m_G = geometry.getInterpolatedGeometryDerivatives(GLL2G);
    m_mine = true;
  }

  m_division = geometry.getDivisionInfo(size);
  m_W   = utilities::g_manager.aquireMatrixStack("temp W"  ,D.rows()-2,D.cols(),m_division[rank].elements.size());
  m_Wt  = utilities::g_manager.aquireMatrixStack("temp W1" ,D.rows()-2,D.cols()-2,m_division[rank].elements.size());
  m_buf = utilities::g_manager.aquireMatrixStack("temp buf",D.rows()-2,D.cols()-2,m_division[rank].elements.size());
  m_view = new basics::componentView<basics::Matrix,basics::matrixStack>(m_G);
  for( int i=0;i<geometry.size();++i )
    m_evals.push_back(new deformedDivergenceEvaluator((*m_view)[i],m_Dt,m_GLL2G,
          m_weight,(*m_W)[i],(*m_Wt)[i],
          (*m_buf)[i]));
}

SEMDivergenceEvaluator::~SEMDivergenceEvaluator()
{
  if( m_mine ) {
    for( int i=0;i<m_G.size();++i )
      delete m_G[i];
  }
  for( int i=0;i<m_evals.size();++i )
    delete m_evals[i];
  delete m_view;
  utilities::g_manager.unlock(m_W);
  utilities::g_manager.unlock(m_Wt);
  utilities::g_manager.unlock(m_buf);
}

void SEMDivergenceEvaluator::evaluate(basics::matrixStack& res, const basics::Field2<basics::matrixStack>& u) const
{
  int max=m_division[m_rank].elements.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    const basics::Field2<basics::Matrix> foo(const_cast<basics::Matrix*>(&u.X()[i]),
        const_cast<basics::Matrix*>(&u.Y()[i]));
    m_evals[i]->evaluate(res[i],foo);
  }
}

extrudedDivergenceEvaluator::extrudedDivergenceEvaluator(const SEMDivergenceEvaluator& eval,
    int rank, int size) :
  m_eval(eval), m_rank(rank), m_size(size)
{
  m_division = eval.m_geometry.getDivisionInfo(size);
}

extrudedDivergenceEvaluator::~extrudedDivergenceEvaluator()
{
}

void extrudedDivergenceEvaluator::evaluate(basics::matricesStack& res, 
    const basics::Field3<basics::matricesStack>& u) const
{
  basics::matricesStack* temp = 
    utilities::g_manager.aquireMatricesStack("temporary stack",
        res[0].rows(),res[0].cols(),
        u.X()[0].matrices(),u.X().size());
  //    int layers = m_eval.m_geometry.size()/192;
  Real Lz = m_eval.m_geometry.m_Lz;//layers;

  /* X and Y contributions */
  int max=res.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    for( int l=0;l<u.X()[0].matrices();++l ) {
      int eval = m_division[m_rank].elements[i];
      const basics::Field2<basics::Matrix> src(const_cast<basics::Matrix*>(&u.X()[i][l]),
          const_cast<basics::Matrix*>(&u.Y()[i][l]));
      m_eval.m_evals[eval]->evaluate((*temp)[i][l],src,false);
    }
    applyLocalGlobal(res[i],(*temp)[i],m_eval.m_GLL2G,'N','T');
  }

  /* Z contribution */
  const vector<basics::matrixStack*>& G = m_eval.m_G;
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    int elem = m_division[m_rank].elements[i];
    for( int l=0;l<(*temp)[i].matrices();++l ) {
      basics::multTranspose((*m_eval.m_W)[i],m_eval.m_GLL2G,u.Z()[i][l],'N','N');
      basics::multTranspose((*temp)[i][l],(*m_eval.m_W)[i],m_eval.m_GLL2G,'N','T');
      for( int j=0;j<res[i].cols();++j )
        for( int k=0;k<res[i].rows();++k )
          (*temp)[i][l][j][k] *= 	(*G[0])[elem][j][k]*(*G[3])[elem][j][k]-
            (*G[2])[elem][j][k]*(*G[1])[elem][j][k];
    }
    applyLocalGlobal(res[i],(*temp)[i],m_eval.m_Dt,'N','T',0,
        Real(2)/(Lz),mnlRealOne);
  }
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    for( int l=0;l<res[i].matrices();++l ) {
      for(int j=0;j<res[i].cols();++j )
        res[i][l][j] *= m_eval.m_weight[l]*m_eval.m_weight[j]*Lz/2;
      for( int k=0;k<res[i].rows();++k )
        res[i][l].scaleRow(k,m_eval.m_weight[k]);
    }
  }

  utilities::g_manager.unlock(temp);
}

SEM3DDivergenceEvaluator::SEM3DDivergenceEvaluator(const basics::geometryStack3D& geometry,
    const basics::Matrix& D,
    const basics::Matrix& GLL2G,
    const basics::Vector& weight,
    int rank, int size) :
  m_geometry(geometry), m_weight(weight), m_GLL2G(GLL2G),
  m_Dt("interpolate derivative matrix", D.rows()-2,D.cols()),
  m_rank(rank), m_size(size)
{
  basics::multTranspose(m_Dt,m_GLL2G,D,'N','N');
  m_G = geometry.getInterpolatedReferenceGeometryDerivatives(GLL2G);
  m_division = geometry.getDivisionInfo(size);
  m_W   = utilities::g_manager.aquireMatrixStack("temp W"  ,D.rows()-2,D.cols(),m_division[rank].elements.size());
  m_Wt  = utilities::g_manager.aquireMatricesStack("temp W1" ,D.rows()-2,D.cols()-2,D.cols(),m_division[rank].elements.size());
  m_buf = utilities::g_manager.aquireMatricesStack("temp buf",D.rows()-2,D.cols()-2,D.cols()-2,m_division[rank].elements.size());
  m_view = new basics::componentView<mnl::basics::Matrices,mnl::basics::matricesStack>(m_G);
  for( int i=0;i<m_division[rank].elements.size();++i ) {
    int eval = m_division[rank].elements[i];
    m_evals.push_back(new deformed3DDivergenceEvaluator((*m_view)[eval],m_Dt,m_GLL2G,
          m_weight,(*m_W)[i],(*m_Wt)[i],
          (*m_buf)[i]));
  }
}

SEM3DDivergenceEvaluator::~SEM3DDivergenceEvaluator()
{
  for( int i=0;i<m_G.size();++i )
    delete m_G[i];
  for( int i=0;i<m_evals.size();++i )
    delete m_evals[i];
  delete m_view;
  utilities::g_manager.unlock(m_W);
  utilities::g_manager.unlock(m_Wt);
  utilities::g_manager.unlock(m_buf);
}

void SEM3DDivergenceEvaluator::evaluate(basics::matricesStack& res,
    const basics::Field3<basics::matricesStack>& u) const
{
  int max=m_division[m_rank].elements.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    const basics::Field3<basics::Matrices> foo(const_cast<basics::Matrices*>(&u.X()[i]),
        const_cast<basics::Matrices*>(&u.Y()[i]),
        const_cast<basics::Matrices*>(&u.Z()[i]));
    m_evals[i]->evaluate(res[i],foo);
  }
}

