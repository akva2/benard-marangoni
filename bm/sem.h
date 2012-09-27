#ifndef SEM_H_
#define SEM_H_

#include <sstream>
#include <numeric>

#ifdef HAS_MPI
#include <mpi.h>
#endif

#include "mass.h"
#include "laplacian.h"

void removeHydrostaticMode(mnl::basics::Matrix& res,
    const SEMInverseMassEvaluator& eval, int i);
void removeHydrostaticMode(mnl::basics::matrixStack& res,
    const SEMInverseMassEvaluator& eval);
std::string errorReport(mnl::basics::matrixStack& u, 
    const mnl::basics::geometryStack& geometry,
    const mnl::basics::function2D& eu,
    const mnl::basics::Vector& grid, Real t, 
    const SEMLaplacianEvaluator* eval=NULL);
std::string errorReport(mnl::basics::matricesStack& u, 
    const mnl::basics::geometryStack& geometry,
    const mnl::basics::function3D& eu,
    const mnl::basics::Vector& grid, Real t, 
    const extrudedLaplacianEvaluator* eval=NULL);
std::string errorReport(mnl::basics::matricesStack& u, 
    const mnl::basics::geometryStack3D& geometry,
    const mnl::basics::function3D& eu,
    const mnl::basics::Vector& grid, Real t, 
    const SEM3DLaplacianEvaluator* eval=NULL);
std::string errorReport(mnl::basics::Field2<mnl::basics::matrixStack>& u, 
    const mnl::basics::geometryStack& geometry,
    const mnl::basics::function2D& eu, 
    const mnl::basics::function2D& ey,
    const mnl::basics::Vector& grid, Real t, 
    const SEMLaplacianEvaluator* eval2=NULL);
std::string errorReport(mnl::basics::Field2<mnl::basics::matrixStack>& u, 
    mnl::basics::matrixStack& p,
    const mnl::basics::geometryStack& geometry,
    const mnl::basics::function2D& eu, 
    const mnl::basics::function2D& ey,
    const mnl::basics::function2D& ep,
    const mnl::basics::Vector& grid, 
    const mnl::basics::Vector& gridGL, Real t, 
    const SEMInverseMassEvaluator& eval,
    const SEMLaplacianEvaluator* eval2=NULL);
std::string errorReport(mnl::basics::Field3<mnl::basics::matricesStack>& u, 
    mnl::basics::matricesStack& p,
    const mnl::basics::geometryStack& geometry,
    const mnl::basics::function3D& eu, 
    const mnl::basics::function3D& ey,
    const mnl::basics::function3D& ez,
    const mnl::basics::function3D& ep,
    const mnl::basics::Vector& grid, 
    const mnl::basics::Vector& gridGL, Real t, 
    const extrudedInverseMassEvaluator& eval,
    const extrudedLaplacianEvaluator* eval2=NULL);
std::string errorReport(mnl::basics::Field3<mnl::basics::matricesStack>& u, 
    mnl::basics::matricesStack& p,
    const mnl::basics::geometryStack3D& geometry,
    const mnl::basics::function3D& eu, 
    const mnl::basics::function3D& ey,
    const mnl::basics::function3D& ez,
    const mnl::basics::function3D& ep,
    const mnl::basics::Vector& grid, 
    const mnl::basics::Vector& gridGL, Real t, 
    const SEM3DInverseMassEvaluator& eval,
    const SEM3DLaplacianEvaluator* eval2=NULL);
std::string errorReport(mnl::basics::Field3<mnl::basics::matricesStack>& u, 
    const mnl::basics::geometryStack& geometry,
    const mnl::basics::function3D& eu, 
    const mnl::basics::function3D& ey,
    const mnl::basics::function3D& ez,
    const mnl::basics::Vector& grid, Real t, 
    const extrudedLaplacianEvaluator* eval2=NULL);

void G2GLL(mnl::basics::matricesStack& u, const mnl::basics::matricesStack& p,
    const mnl::basics::Vector& gridGLL, const mnl::basics::Vector& gridGL);

Real doReduce(Real input, int rank, int size, int tag);
void vectorSum(mnl::basics::Vector& result, int rank, int size, int tag);
void sendStack(mnl::basics::matricesStack& result, const mnl::basics::matricesStack& input, 
    const int* scount, const int* sdispl, 
    const int* rcount, const int* rdispl,
    SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET,
    std::vector<mnl::basics::matrixStack*>* stack=NULL,
    std::vector< std::vector<int> >* planes=NULL);
std::vector<mnl::basics::matrixStack*> sendAndSetupStack(mnl::basics::matricesStack& result, 
    const mnl::basics::matricesStack& input,
    SEMLaplacianEvaluator::BC bc,
    int rank, int size,
    const int* scount,
    const int* sdispl,
    const int* rcount,
    const int* rdispl,
    std::vector< std::vector<int> >* planes=NULL);
void sendStack(mnl::basics::Field3<mnl::basics::matricesStack>& result,
    mnl::basics::matricesStack& buffer,
    const mnl::basics::matricesStack& input, 
    const int* scount, const int* sdispl, 
    const int* rcount, const int* rdispl,
    SEMLaplacianEvaluator::BC bc1=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET,
    SEMLaplacianEvaluator::BC bc2=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET,
    SEMLaplacianEvaluator::BC bc3=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET,
    std::vector<mnl::basics::matrixStack*>* stack=NULL,
    std::vector< std::vector<int> >* planes=NULL);
std::vector<mnl::basics::matrixStack*> sendAndSetupStack(mnl::basics::matricesStack& result, 
    mnl::basics::matricesStack& buffer,
    const mnl::basics::Field3<mnl::basics::matricesStack>& input,
    SEMLaplacianEvaluator::BC bc1,
    SEMLaplacianEvaluator::BC bc2,
    SEMLaplacianEvaluator::BC bc3,
    int rank, int size,
    const int* scount,
    const int* sdispl,
    const int* rcount,
    const int* rdispl,
    std::vector< std::vector<int> >* planes=NULL);

void dssum(mnl::basics::matricesStack& result, const mnl::basics::geometryStack& geometry,
    int rank, int size, SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
void dssum(mnl::basics::matricesStack& result, const mnl::basics::geometryStack3D& geometry,
    int rank, int size, SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
void dssumCoarse(mnl::basics::matricesStack& result, int N,
    const mnl::basics::geometryStack& geometry,
    int rank, int size, bool neumann=false);

void mask(mnl::basics::matricesStack& result,
    const mnl::basics::geometryStack& geometry, int rank=0, int size=1,
    SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);

void mask(mnl::basics::matricesStack& result,
    const mnl::basics::geometryStack3D& geometry, int rank=0, int size=1,
    SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);

void mask(mnl::basics::Field3<mnl::basics::matricesStack>& result, 
    const mnl::basics::geometryStack& geometry, int rank=0, int size=1,
    SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
void mask(mnl::basics::Field3<mnl::basics::matricesStack>& result, 
    const mnl::basics::geometryStack3D& geometry, int rank=0, int size=1,
    SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);

  template<class T>
void dssum(mnl::basics::Field3<mnl::basics::matricesStack>& result,
    const T& geometry, int rank, int size,
    SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET)
{
  dssum(result.X(),geometry,rank,size,bc);
  dssum(result.Y(),geometry,rank,size,bc);
  dssum(result.Z(),geometry,rank,size,bc);
}

void massReference(mnl::basics::Matrix& res, const mnl::basics::Vector& weight);
void massReference(mnl::basics::Matrices& res, const mnl::basics::Vector& weight);

class MPIDotter {
  public:
    MPIDotter(int rank, int size, int tag=0) :
      m_rank(rank), m_size(size), m_tag(tag)
  {
  }

    template<class T>
      Real operator()(const T& u, const T& u2) const
      {
        Real result=0;
        for( int i=0;i<u.size();++i )
          result += u[i].getDotter()(u[i],u2[i]);
        Real resultTotal;
#ifdef HAS_MPI
        MPI_Allreduce(&result,&resultTotal,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
        resultTotal = doReduce(result,m_rank,m_size,m_tag);
#endif
        return( resultTotal );
      }
    int m_tag;
    int m_rank;
    int m_size;
};

template<class G>
class MPIGeometryDotter {
  public:
    MPIGeometryDotter(const G& geometry, const std::vector<int>& r,
        int rank, int size, int tag=0) :
      m_geometry(geometry), m_range(r), m_size(size), m_rank(rank), m_tag(tag)
  {
  }

    template<class T>
      Real operator()(const T& u, const T& u2) const
      {
        Real result = m_geometry.getDotter()(u,u2,m_range);
        Real resultTotal=0;
#ifdef HAS_MPI
        MPI_Allreduce(&result,&resultTotal,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
        resultTotal = result;
#endif

        return( resultTotal );
      }

    int m_tag;
    int m_rank;
    int m_size;
  protected:
    const G& 				m_geometry;
    const std::vector<int>&	m_range;
};

  template<class T>
void removeHydrostaticMode(mnl::basics::matricesStack& res, 
    const T& eval)
{
  mnl::basics::matricesStack* temp  = mnl::utilities::g_manager.clone(res);
  mnl::basics::matricesStack* temp2 = mnl::utilities::g_manager.clone(res);

  *temp = 1;
  eval.evaluateMass(*temp2,*temp);
  MPIDotter dotter(eval.m_rank,eval.m_size);
  Real area = dotter(*temp,*temp2);

  eval.evaluateMass(*temp2,res);
  Real mean = dotter(*temp,*temp2);
  mean /= area;

  res -= mean;

  mnl::utilities::g_manager.unlock(temp);
  mnl::utilities::g_manager.unlock(temp2);

}

#endif

