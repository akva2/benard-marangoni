#include "gradient.h"
#include "sem.h"

using namespace mnl;
using namespace std;

deformedGradientEvaluator::deformedGradientEvaluator(const basics::matrixStack& G,
    const basics::Matrix& Dt,
    const basics::Matrix& GLL2G,
    const basics::Vector& weight,
    basics::Matrix& W,
    basics::Matrix& t,
    basics::Matrix& s) :
  m_weight(weight), m_GLL2G(GLL2G), m_G(G), m_Dt(Dt),
  m_W(W), m_t(t), m_s(s)
{
}

void deformedGradientEvaluator::evaluate(basics::Field2<basics::Matrix>& res,
    const basics::Matrix& p) const
{
  /* D^T_x */
  m_W = p;
  massReference(m_W,m_weight);

  basics::multPointwise(m_t,m_W,m_G[3]); // xi x
  basics::multTranspose(m_s,m_Dt,m_t,'T','N');
  basics::multTranspose(res.X(),m_s,m_GLL2G,'N','N');
  basics::multPointwise(m_t,m_W,m_G[2]); // eta x
  basics::multTranspose(m_s,m_GLL2G,m_t,'T','N');
  basics::multTranspose(res.X(),m_s,m_Dt,'N','N',mnlRealMinusOne,mnlRealOne);

  /* D^T_y */
  basics::multPointwise(m_t,m_W,m_G[1]); // xi y
  basics::multTranspose(m_s,m_Dt,m_t,'T','N');
  basics::multTranspose(res.Y(),m_s,m_GLL2G,'N','N');
  basics::multPointwise(m_t,m_W,m_G[0]); // eta y
  basics::multTranspose(m_s,m_GLL2G,m_t,'T','N');
  basics::multTranspose(res.Y(),m_s,m_Dt,'N','N',mnlRealOne,mnlRealMinusOne);
}

deformed3DGradientEvaluator::deformed3DGradientEvaluator(const basics::matricesStack& G,
    const basics::Matrix& Dt,
    const basics::Matrix& GLL2G,
    const basics::Vector& weight,
    basics::Matrices& W,
    basics::Matrices& t,
    basics::Matrix& s,
    basics::Matrices& v) :
  m_weight(weight), m_GLL2G(GLL2G), m_G(G), m_Dt(Dt),
  m_W(W), m_t(t), m_s(s), m_v(v)
{
}

void deformed3DGradientEvaluator::evaluate(basics::Field3<basics::Matrices>& res,
    const basics::Matrices& p) const
{
  m_W = p;
  massReference(m_W,m_weight);
  evaluateComponent(res.X(),0);
  evaluateComponent(res.Y(),1);
  evaluateComponent(res.Z(),2);
}

void deformed3DGradientEvaluator::evaluateComponent(basics::Matrices& res, int index) const
{
  basics::multPointwise(m_t,m_W,m_G[index]);
  for( int l=0;l<m_t.matrices();++l ) {
    basics::multTranspose(m_s,m_Dt,m_t[l],'T','N');
    basics::multTranspose(m_v[l],m_s,m_GLL2G,'N','N');
  }
  basics::applyLocalGlobal(res,m_v,m_GLL2G,'N','N');

  basics::multPointwise(m_t,m_W,m_G[index+3]);
  for( int l=0;l<m_t.matrices();++l ) {
    basics::multTranspose(m_s,m_GLL2G,m_t[l],'T','N');
    basics::multTranspose(m_v[l],m_s,m_Dt,'N','N');
  }
  basics::applyLocalGlobal(res,m_v,m_GLL2G,'N','N',0,mnlRealOne,mnlRealOne);

  basics::multPointwise(m_t,m_W,m_G[index+6]);
  for( int l=0;l<m_t.matrices();++l ) {
    basics::multTranspose(m_s,m_GLL2G,m_t[l],'T','N');
    basics::multTranspose(m_v[l],m_s,m_GLL2G,'N','N');
  }
  basics::applyLocalGlobal(res,m_v,m_Dt,'N','N',0,mnlRealOne,mnlRealOne);
}

SEMGradientEvaluator::SEMGradientEvaluator(const basics::geometryStack& geometry,
    const basics::Matrix& D,
    const basics::Matrix& GLL2G,
    const basics::Vector& weight,
    const std::vector<basics::matrixStack*>* G,
    int rank, int size) :
  m_weight(weight), m_GLL2G(GLL2G), m_geometry(geometry), m_D(D),
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
  m_W   = utilities::g_manager.aquireMatrixStack("temp W",m_GLL2G.rows(), m_GLL2G.rows(),geometry.size());
  m_t  = utilities::g_manager.clone(*m_W);
  m_s  = utilities::g_manager.aquireMatrixStack("s",m_Dt.cols(),(*m_W)[0].cols(),geometry.size());
  m_division = m_geometry.getDivisionInfo(size);
  m_view = new basics::componentView<basics::Matrix,basics::matrixStack>(m_G);
  for( int i=0;i<geometry.size();++i )
    m_evals.push_back(new deformedGradientEvaluator((*m_view)[i],m_Dt,m_GLL2G,m_weight,
          (*m_W)[i],(*m_t)[i],(*m_s)[i]));
}

SEMGradientEvaluator::~SEMGradientEvaluator()
{
  if( m_mine ) {
    for( int i=0;i<m_G.size();++i )
      delete m_G[i];
  }
  for( int i=0;i<m_evals.size();++i )
    delete m_evals[i];
  delete m_view;
  utilities::g_manager.unlock(m_W);
  utilities::g_manager.unlock(m_t);
  utilities::g_manager.unlock(m_s);
}

void SEMGradientEvaluator::evaluate(basics::Field2<basics::matrixStack>& res, 
    const basics::matrixStack& u) const
{
  int max=m_division[m_rank].elements.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    basics::Field2<basics::Matrix> foo(&res.X()[i],&res.Y()[i]);
    m_evals[i]->evaluate(foo,u[i]);
  }
}

extrudedGradientEvaluator::extrudedGradientEvaluator(const SEMGradientEvaluator& eval, 
    int rank, int size) :
  m_eval(eval), m_rank(rank), m_size(size)
{
  m_division = m_eval.m_geometry.getDivisionInfo(size);
}

extrudedGradientEvaluator::~extrudedGradientEvaluator()
{
}

void extrudedGradientEvaluator::evaluate(basics::Field3<basics::matricesStack>& res, 
    const basics::matricesStack& p) const
{
  basics::Field3<basics::matricesStack>* temp = 
    utilities::g_manager.aquireMatricesStackField("temporary stack",res.X()[0].rows(),
        res.X()[0].cols(),
        p[0].matrices(),p.size());
  /* X and Y directions */
  int max=p.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    for( int l=0;l<p[0].matrices();++l ) {
      basics::Field2<basics::Matrix> dest(&temp->X().at(l)[i],&temp->Y().at(l)[i]);
      m_eval.m_evals[m_division[m_rank].elements[i]]->evaluate(dest,p.at(l)[i]);
      dest *= m_eval.m_weight[l]*m_eval.m_geometry.m_Lz/2;
    }
    applyLocalGlobal(res.X()[i],temp->X()[i],m_eval.m_GLL2G,'N','N');
    applyLocalGlobal(res.Y()[i],temp->Y()[i],m_eval.m_GLL2G,'N','N');
  }

  /* Z direction */
  const std::vector<basics::matrixStack*>& G = m_eval.m_G;
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    int elem = m_division[m_rank].elements[i];
    for( int l=0;l<p[i].matrices();++l ) {
      for( int j=0;j<p[i].cols();++j )
        for( int k=0;k<p[i].rows();++k )
          (*m_eval.m_W)[i][j][k] = p[i][l][j][k]*((*G[0])[elem][j][k]*(*G[3])[elem][j][k]-
              (*G[1])[elem][j][k]*(*G[2])[elem][j][k])
            *m_eval.m_weight[j]
            *m_eval.m_weight[k];

      basics::multTranspose((*m_eval.m_s)[i],m_eval.m_GLL2G,(*m_eval.m_W)[i],'T','N');
      basics::multTranspose(temp->Z().at(l)[i],(*m_eval.m_s)[i],m_eval.m_GLL2G,'N','N',
          m_eval.m_weight[l]);
    }
    applyLocalGlobal(res.Z()[i],temp->Z()[i],m_eval.m_Dt,'N','N',0);
  }

  utilities::g_manager.unlock(temp);
}

SEM3DGradientEvaluator::SEM3DGradientEvaluator(const basics::geometryStack3D& geometry,
    const basics::Matrix& D,
    const basics::Matrix& GLL2G,
    const basics::Vector& weight,
    const std::vector<basics::matricesStack*>* G,
    int rank, int size) :
  m_weight(weight), m_GLL2G(GLL2G), m_geometry(geometry), m_D(D),
  m_Dt("interpolate derivative matrix", D.rows()-2,D.cols()),
  m_rank(rank), m_size(size)
{
  basics::multTranspose(m_Dt,m_GLL2G,D,'N','N');

  if( G ) {
    m_G = *G;
    m_mine = false;
  }
  else  {
    m_G = geometry.getInterpolatedReferenceGeometryDerivatives(GLL2G);
    m_mine = true;
  }
  m_division = m_geometry.getDivisionInfo(size);
  m_W   = utilities::g_manager.aquireMatricesStack("temp W",m_GLL2G.rows(), m_GLL2G.rows(),m_GLL2G.rows(),m_division[m_rank].elements.size());
  m_t  = utilities::g_manager.clone(*m_W);
  m_s  = utilities::g_manager.aquireMatrixStack("s",m_Dt.cols(),(*m_W)[0].cols(),m_division[m_rank].elements.size());
  m_v  = utilities::g_manager.aquireMatricesStack("s",m_Dt.cols(),m_Dt.cols(),D.cols()-2,m_division[m_rank].elements.size());
  m_view = new basics::componentView<basics::Matrices,basics::matricesStack>(m_G);
  for( int i=0;i<m_division[rank].elements.size();++i ) {
    int eval=m_division[rank].elements[i];
    m_evals.push_back(new deformed3DGradientEvaluator((*m_view)[eval],m_Dt,m_GLL2G,m_weight,
          (*m_W)[i],(*m_t)[i],(*m_s)[i],(*m_v)[i]));
  }
}

SEM3DGradientEvaluator::~SEM3DGradientEvaluator()
{
  if( m_mine ) {
    for( int i=0;i<m_G.size();++i )
      delete m_G[i];
  }
  delete m_view;
  for( int i=0;i<m_evals.size();++i )
    delete m_evals[i];

  utilities::g_manager.unlock(m_W);
  utilities::g_manager.unlock(m_t);
  utilities::g_manager.unlock(m_s);
  utilities::g_manager.unlock(m_v);
}

void SEM3DGradientEvaluator::evaluate(basics::Field3<basics::matricesStack>& res, 
    const basics::matricesStack& u) const
{
  int max=m_division[m_rank].elements.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    basics::Field3<basics::Matrices> foo(&res.X()[i],&res.Y()[i],&res.Z()[i]);
    m_evals[i]->evaluate(foo,u[i]);
  }
}
