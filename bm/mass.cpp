#include "mass.h"
#include "sem.h"

using namespace mnl;
using namespace std;

SEMMassEvaluator::SEMMassEvaluator(const basics::geometryStack& geometry,
    const basics::Vector& weight,
    const basics::Matrix& GLL2G,
    int rank, int size) :
  m_weight(weight), m_geometry(geometry), m_rank(rank), m_size(size)
{
  /* setup jacobian on GL grid */
  basics::Matrix temp ("temp",weight.length(),weight.length()+2);
  basics::Matrix temp2("temp",weight.length(),weight.length());
  basics::matrixStack G("temp");
  G.add(temp2,4);
  m_division = geometry.getDivisionInfo(size);
  m_J = utilities::g_manager.aquireMatrixStack("interpolated jacobian",
      GLL2G.rows(),GLL2G.rows(),geometry.size());
  for( int i=0;i<geometry.size();++i ) {
    const basics::matrixStack& D = geometry[i].getGH().getGeometryDerivatives();
    basics::multTranspose(temp,GLL2G,D[0],'N','N');
    basics::multTranspose(G[0],temp,GLL2G,'N','T');
    basics::multTranspose(temp,GLL2G,D[1],'N','N');
    basics::multTranspose(G[1],temp,GLL2G,'N','T');
    basics::multTranspose(temp,GLL2G,D[2],'N','N');
    basics::multTranspose(G[2],temp,GLL2G,'N','T');
    basics::multTranspose(temp,GLL2G,D[3],'N','N');
    basics::multTranspose(G[3],temp,GLL2G,'N','T');
    for( int j=0;j<G[0].cols();++j )
      for( int k=0;k<G[0].rows();++k )
        (*m_J)[i][j][k] = m_weight[j]*m_weight[k]*(G[0][j][k]*G[3][j][k]-
            G[2][j][k]*G[1][j][k]);
  }
}

SEMMassEvaluator::~SEMMassEvaluator()
{
  utilities::g_manager.unlock(m_J);
}

void SEMMassEvaluator::evaluate(basics::matrixStack& res, const basics::matrixStack& u) const
{
  int max=res.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i )
    basics::multPointwise(res[i],u[i],(*m_J)[m_division[m_rank].elements[i]]);
}

SEMInverseMassEvaluator::SEMInverseMassEvaluator(const basics::geometryStack& geometry,
    const basics::Vector& weight,
    const basics::Matrix& GLL2G,
    int rank, int size) :
  SEMMassEvaluator(geometry,weight,GLL2G,rank,size)
{
  m_iJ = utilities::g_manager.clone(*m_J);
  for( int i=0;i<m_J->size();++i )
    for( int j=0;j<(*m_J)[i].cols();++j )
      for( int k=0;k<(*m_J)[i].cols();++k )
        (*m_iJ)[i][j][k] = Real(1)/(*m_J)[i][j][k];
}

SEMInverseMassEvaluator::~SEMInverseMassEvaluator()
{
  utilities::g_manager.unlock(m_iJ);
}

void SEMInverseMassEvaluator::evaluate(basics::matrixStack& res, const basics::matrixStack& u) const
{
  int max=res.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i )
    basics::multPointwise(res[i],u[i],(*m_iJ)[m_division[m_rank].elements[i]]);
}

void SEMInverseMassEvaluator::evaluateMass(basics::matrixStack& res, const basics::matrixStack& u) const
{
  SEMMassEvaluator::evaluate(res,u);
}

SEM3DMassEvaluator::SEM3DMassEvaluator(const basics::geometryStack3D& geometry,
    const basics::Vector& weight,
    const basics::Matrix& GLL2G,
    int rank, int size) :
  m_weight(weight), m_geometry(geometry), m_rank(rank), m_size(size)
{
  std::vector<basics::matricesStack*> D = m_geometry.getInterpolatedGeometryDerivatives(GLL2G);
  m_division = geometry.getDivisionInfo(size);
  m_J = utilities::g_manager.clone(*D[0]);
#pragma omp parallel for schedule(static)
  for( int i=0;i<geometry.size();++i ) {
    for( int gamma=0;gamma<(*D[0])[0].matrices();++gamma ) {
      for( int beta=0;beta<(*D[0])[0].cols();++beta) {
        for( int alpha=0;alpha<(*D[0])[0].rows();++alpha ) {
          Real xxi    = (*D[0])[i][gamma][beta][alpha];
          Real xeta   = (*D[1])[i][gamma][beta][alpha];
          Real xgamma = (*D[2])[i][gamma][beta][alpha];
          Real yxi    = (*D[3])[i][gamma][beta][alpha];
          Real yeta   = (*D[4])[i][gamma][beta][alpha];
          Real ygamma = (*D[5])[i][gamma][beta][alpha];
          Real zxi    = (*D[6])[i][gamma][beta][alpha];
          Real zeta   = (*D[7])[i][gamma][beta][alpha];
          Real zgamma = (*D[8])[i][gamma][beta][alpha];
          (*m_J)[i][gamma][beta][alpha] =   (xxi*(yeta*zgamma-zeta*ygamma)
              -xeta*(yxi*zgamma-zxi*ygamma)
              +xgamma*(yxi*zeta-zxi*yeta))
            *weight[gamma]*weight[beta]*weight[alpha];
        }
      }
    }
  }
  for( int i=0;i<D.size();++i )
    delete D[i];
}

SEM3DMassEvaluator::~SEM3DMassEvaluator()
{
  utilities::g_manager.unlock(m_J);
}

void SEM3DMassEvaluator::evaluate(basics::matricesStack& res, const basics::matricesStack& u) const
{
  int max=res.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i )
    basics::multPointwise(res[i],u[i],(*m_J)[m_division[m_rank].elements[i]]);
}

SEM3DInverseMassEvaluator::SEM3DInverseMassEvaluator(const basics::geometryStack3D& geometry,
    const basics::Vector& weight,
    const basics::Matrix& GLL2G,
    int rank, int size) :
  SEM3DMassEvaluator(geometry,weight,GLL2G,rank,size)
{
  m_iJ = utilities::g_manager.clone(*m_J);
  for( int i=0;i<m_J->size();++i )
    for( int l=0;l<(*m_J)[i].matrices();++l )
      for( int j=0;j<(*m_J)[i].cols();++j )
        for( int k=0;k<(*m_J)[i].cols();++k )
          (*m_iJ)[i][l][j][k] = Real(1)/(*m_J)[i][l][j][k];
}

SEM3DInverseMassEvaluator::~SEM3DInverseMassEvaluator()
{
  utilities::g_manager.unlock(m_iJ);
}

void SEM3DInverseMassEvaluator::evaluate(basics::matricesStack& res, const basics::matricesStack& u) const
{
  int max=res.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i )
    basics::multPointwise(res[i],u[i],(*m_iJ)[m_division[m_rank].elements[i]]);
}

void SEM3DInverseMassEvaluator::evaluateMass(basics::matricesStack& res, const basics::matricesStack& u) const
{
  SEM3DMassEvaluator::evaluate(res,u);
}

extrudedMassEvaluator::extrudedMassEvaluator(const SEMMassEvaluator& eval, int rank, int size) :
  m_eval(eval), m_rank(rank), m_size(size)
{
  m_division = m_eval.m_geometry.getDivisionInfo(size);
}

void extrudedMassEvaluator::evaluate(basics::matricesStack& res, const basics::matricesStack& p) const
{
  int max=res.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    for( int l=0;l<res[i].matrices();++l ) {
      basics::multPointwise(res[i][l],p[i][l],(*m_eval.m_J)[m_division[m_rank].elements[i]]);
      res[i][l] *= (m_eval.m_weight[l]*m_eval.m_geometry.m_Lz/2);
    }
  }
}

extrudedInverseMassEvaluator::extrudedInverseMassEvaluator(const SEMInverseMassEvaluator& eval, int rank, int size) :
  extrudedMassEvaluator(eval,rank,size), m_invEval(eval)
{
}

void extrudedInverseMassEvaluator::evaluate(basics::matricesStack& res, const basics::matricesStack& p) const
{
  int max=res.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    for( int l=0;l<res[i].matrices();++l ) {
      basics::multPointwise(res[i][l],p[i][l],(*m_invEval.m_iJ)[m_division[m_rank].elements[i]]);
      res[i][l] *= Real(1)/(m_eval.m_weight[l]*m_eval.m_geometry.m_Lz/2);
    }
  }
}

void extrudedInverseMassEvaluator::evaluateMass(basics::matricesStack& res, const basics::matricesStack& p) const
{
  extrudedMassEvaluator::evaluate(res,p);
}

