#include "convection.h"
#include "sem.h"

using namespace mnl;
using namespace std;

deformedConvectionEvaluator::deformedConvectionEvaluator(const basics::matrixStack& G,
    const basics::Matrix& D,
    const basics::Vector& weight,
    const basics::Field2<basics::Matrix>& u,
    basics::Matrix& buffer,
    basics::Matrix& buffer2) :
  m_weight(weight), m_D(D), m_buffer(buffer), m_buffer2(buffer2), m_G(G), m_u(u)
{
}

void deformedConvectionEvaluator::evaluate(basics::Matrix& result, const basics::Matrix& u) const
{
  /* ux */
  basics::multTranspose(m_buffer, m_D,u,'N','N');
  basics::multPointwise(m_buffer2,m_buffer,m_G[3]);
  basics::multTranspose(m_buffer,u,m_D,'N','T');
  basics::multPointwise(m_buffer2,m_buffer,m_G[2],1,-1);
  basics::multPointwise(result,m_buffer2,m_u.X());

  /* uy */
  basics::multTranspose(m_buffer, u,m_D,'N','T');
  basics::multPointwise(m_buffer2,m_buffer,m_G[0]);
  basics::multTranspose(m_buffer,m_D,u,'N','N');
  basics::multPointwise(m_buffer2,m_buffer,m_G[1],1,-1);
  basics::multPointwise(result,m_buffer2,m_u.Y(),1,1);

  result *= -1;
  massReference(result,m_weight);
}

deformed3DConvectionEvaluator::
deformed3DConvectionEvaluator(const basics::matricesStack& G,
    const basics::Matrix& D,
    const basics::Vector& weight,
    const basics::Field3<basics::Matrices>& u,
    basics::Matrices& buffer,
    basics::Matrices& buffer2) :
  m_weight(weight), m_D(D), m_buffer(buffer), m_buffer2(buffer2), m_G(G), m_u(u)
{
}

void deformed3DConvectionEvaluator::evaluate(basics::Matrices& res,
    const basics::Matrices& u) const
{
  evaluateComponent(res,u,0);
  basics::multPointwise(res,m_u.X());
  evaluateComponent(m_buffer2,u,1);
  basics::multPointwise(m_buffer2,m_u.Y());
  res += m_buffer2;
  evaluateComponent(m_buffer2,u,2);
  basics::multPointwise(m_buffer2,m_u.Z());
  res += m_buffer2;
  res *= -1;
  massReference(res,m_weight);
}

void deformed3DConvectionEvaluator::evaluateComponent(basics::Matrices& res,
    const basics::Matrices& field,
    int index) const
{
  for( int l=0;l<field.matrices();++l )
    basics::multTranspose(m_buffer[l],m_D,field[l],'N','N');
  basics::multPointwise(res,m_buffer,m_G[index]);
  for( int l=0;l<field.matrices();++l )
    basics::multTranspose(m_buffer[l],field[l],m_D,'N','T');
  basics::multPointwise(m_buffer,m_G[index+3]);
  res += m_buffer;
  basics::applyLocalGlobal(m_buffer,field,m_D,'N','T');
  basics::multPointwise(m_buffer,m_G[index+6]);
  res += m_buffer;
}

convectionEvaluator::convectionEvaluator(const basics::geometryStack& geometry,
    const basics::Matrix& D,
    const basics::Vector& weight,
    int rank, int size,
    const basics::geometryStack3D* geometry3D) :
  m_geometry(geometry), m_weight(weight),
  m_D(D), m_rank(rank), m_size(size), m_view(NULL),
  m_geometry3D(NULL), m_deformed(false), m_order(1)
{
  if( geometry3D ) {
    m_grid = geometry3D->getDivisionInfo(size);
    m_geometry3D = geometry3D;
    m_deformed = true;
    m_view = new basics::componentView<basics::Matrices,basics::matricesStack>(m_geometry3D->getReferenceGeometryDerivatives());
  }
  else
    m_grid = geometry.getDivisionInfo(size);

  m_buffer  = utilities::g_manager.aquireMatrixStack("convection buffer",D.rows(),D.cols(),m_grid[m_rank].elements.size());
  m_buffer2 = utilities::g_manager.clone(*m_buffer);
  m_buffer3 = utilities::g_manager.aquireMatricesStack("convection buffer3D",D.rows(),D.cols(),D.rows(),m_grid[m_rank].elements.size());
  m_buffer4 = utilities::g_manager.clone(*m_buffer3);
}

convectionEvaluator::~convectionEvaluator()
{
  delete m_view;

  utilities::g_manager.unlock(m_buffer);
  utilities::g_manager.unlock(m_buffer2);
  utilities::g_manager.unlock(m_buffer3);
  utilities::g_manager.unlock(m_buffer4);
}

void convectionEvaluator::evaluate(basics::matrixStack& res,
    const basics::matrixStack& u) const
{
  basics::Field2<basics::matrixStack>* field = (basics::Field2<basics::matrixStack>*)m_interp;
  int max=u.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    basics::Field2<basics::Matrix> foo(&field->X()[i],&field->Y()[i]);
    deformedConvectionEvaluator eval(m_geometry[i].getGH().getGeometryDerivatives(),
        m_D,m_weight,foo,(*m_buffer)[i],(*m_buffer2)[i]);
    eval.evaluate(res[i],u[i]);
  }
  if( m_bc == SEMLaplacianEvaluator::PERIODIC ) {
    m_geometry.periodicDssum(res);
    m_geometry.invMassP(res);
  } else {
    m_geometry.dssum(res);
    m_geometry.invMass(res);
    m_geometry.mask(res);
  }
}

void convectionEvaluator::evaluate(basics::matricesStack& res,
    const basics::matricesStack& u)
{
  if( m_deformed ) {
    evaluateDeformed(res,u);
    return;
  }

  /* extruded */
  basics::Field3<basics::matricesStack>* field = (basics::Field3<basics::matricesStack>*)m_interp;

  /* z */
  int max=u.size();
#pragma omp parallel for schedule(static)
  for( int n=0;n<max;++n ) {
    basics::applyLocalGlobal((*m_buffer3)[n],u[n],m_D,'N','T',0,-Real(2)/m_geometry.m_Lz);
    basics::multPointwise(res[n],(*m_buffer3)[n],field->Z()[n]);
  }
  m_geometry.mass(res,m_grid[m_rank].elements);

  /* x and y */
  for( int l=0;l<u[0].matrices();++l ) {
#pragma omp parallel for schedule(static)
    for( int i=0;i<max;++i ) {
      basics::Field2<basics::Matrix> foo(&field->X().at(l)[i],&field->Y().at(l)[i]);
      deformedConvectionEvaluator
        eval(m_geometry[m_grid[m_rank].elements[i]].getGH().getGeometryDerivatives(),
            m_D,m_weight,foo,(*m_buffer)[i],(*m_buffer2)[i]);
      eval.evaluate(m_buffer3->at(l)[i],u.at(l)[i]);
    }
    res.at(l).axpy(m_weight[l]*m_geometry.m_Lz/2,m_buffer3->at(l));
  }
  if( m_bc == SEMLaplacianEvaluator::PERIODIC ) {
    m_geometry.periodicDssum(res);
    m_geometry.invMassP(res);
  } else {
    dssum(res,m_geometry,m_rank,m_size);
    m_geometry.invMass(res,m_grid[m_rank].elements);
    mask(res,m_geometry,m_rank,m_size,m_bc);
  }
}

void convectionEvaluator::evaluateDeformed(basics::matricesStack& res,
    const basics::matricesStack& u)
{
  basics::Field3<basics::matricesStack>* field = (basics::Field3<basics::matricesStack>*)m_interp;
  int max=m_grid[m_rank].elements.size();
#pragma omp parallel for schedule(static)
  for( int i=0;i<max;++i ) {
    int elem=m_grid[m_rank].elements[i];
    basics::Field3<basics::Matrices> foo(&field->X()[i],&field->Y()[i],&field->Z()[i]);
    deformed3DConvectionEvaluator
      eval((*m_view)[elem],m_D,m_weight,foo,(*m_buffer3)[i],(*m_buffer4)[i]);
    eval.evaluate(res[i],u[i]);
  }
  dssum(res,*m_geometry3D,m_rank,m_size);
  m_geometry3D->invMass(res,m_grid[m_rank].elements);
  mask(res,*m_geometry3D,m_rank,m_size,m_bc);
}

