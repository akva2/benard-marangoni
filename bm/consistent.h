#ifndef CONSISTENT_H_
#define CONSISTENT_H_

#include "mnl/geometry.h"
#include "mnl/util.h"

#include "gradient.h"
#include "divergence.h"
#include "mass.h"
#include "sem.h"

class deformedConsistentPressureEvaluator : public 
                                            mnl::utilities::Evaluator<mnl::basics::Matrix> {
                                              public:
                                                deformedConsistentPressureEvaluator(const deformedDivergenceEvaluator& divergence,
                                                    const deformedGradientEvaluator& gradient,
                                                    const mnl::basics::spectralElement2D& elem,
                                                    const mnl::basics::Vector& weight);
                                                virtual ~deformedConsistentPressureEvaluator();

                                                virtual void evaluate(mnl::basics::Matrix& res, const mnl::basics::Matrix& p) const;
                                              protected:
                                                mnl::basics::Field2D* m_W;
                                                const deformedDivergenceEvaluator&    m_divergence;
                                                const deformedGradientEvaluator&   	  m_gradient;
                                                const mnl::basics::spectralElement2D& m_element;
                                                const mnl::basics::Vector&			  m_weight;
                                            };

class SEMConsistentPressureEvaluator : public 
                                       mnl::utilities::Evaluator<mnl::basics::matrixStack> {
                                         public:
                                           SEMConsistentPressureEvaluator(const mnl::basics::geometryStack& geometry,
                                               const mnl::basics::Matrix& D,
                                               const mnl::basics::Matrix& GLL2G,
                                               const mnl::basics::Vector& weightGL,
                                               const mnl::basics::Vector& weight,
                                               mnl::basics::matrixStack* J=NULL,
                                               const Real nu=0,
                                               SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET,
                                               int rank=0, int size=1);
                                           virtual ~SEMConsistentPressureEvaluator();

                                           virtual void evaluate(mnl::basics::matrixStack& res, 
                                               const mnl::basics::matrixStack& p) const;

                                           void evaluateSimple(mnl::basics::matrixStack& res, 
                                               const mnl::basics::matrixStack& p) const;

                                           virtual bool isL2() const
                                           {
                                             return m_nu < 1.e-10;
                                           }

                                           virtual void filter(mnl::basics::matrixStack& res) const
                                           {
                                             res -= res.sum()/res.length();
                                           }

                                           int solve(mnl::basics::matrixStack& res, 
                                               mnl::utilities::Evaluator<mnl::basics::matrixStack>* pre=NULL);

                                           template<class T>
                                             static void getRestriction(mnl::basics::Vector& result,
                                                 const T& p,
                                                 const mnl::basics::coarseGrid& desc,
                                                 int rank=0, int size=1, int tag=0)
                                             {
                                               result.clear();

                                               /* restrict */
                                               for( int i=0;i<desc.size();++i ) { 
                                                 for( int n=0;n<desc[i].elements.size();++n )
                                                   result[desc[i].size3] += p[desc[i].elements[n]].sum();
                                               }
                                               vectorSum(result,rank,size,tag);
                                             }

                                           template<class T>
                                             static void getProlongiation(T& res, 
                                                 const mnl::basics::Vector& temp, 
                                                 const mnl::basics::coarseGrid& desc,
                                                 bool bMasked=true)
                                             {
                                               /* prolong */
                                               for( int i=0;i<desc.size();++i ) {
                                                 for( int n=0;n<desc[i].elements.size();++n ) {
                                                   if( bMasked && desc[i].size3 == temp.length()-1 )
                                                     res[desc[i].elements[n]].clear();
                                                   else
                                                     res[desc[i].elements[n]] = temp[desc[i].size3];
                                                 }
                                               }
                                             }

                                           template<class T>
                                             static void orthogonalize(T& res,
                                                 const mnl::basics::coarseGrid& desc)
                                             {
                                               for( int i=0;i<desc.size();++i ) {
                                                 Real sum=0;
                                                 int length=0;
                                                 for( int n=0;n<desc[i].elements.size();++n ) {
                                                   sum    += res[desc[i].elements[n]].sum();
                                                   length += res[desc[i].elements[n]].length();
                                                 }
                                                 for( int n=0;n<desc[i].elements.size();++n )
                                                   res[desc[i].elements[n]] -= sum/length;
                                               }
                                             }

                                           SEMDivergenceEvaluator 		m_divergence;
                                           const mnl::basics::geometryStack&	m_geometry;
                                           SEMGradientEvaluator   		m_gradient;
                                           const mnl::basics::Vector&	m_weight;
                                           Real						m_nu;
                                           const mnl::basics::Matrix&	m_GLL2G;
                                           mnl::basics::Matrix* 		m_E;
                                           mnl::basics::matrixStack*	m_J;
                                           mnl::basics::matrixStack*	m_buf;
                                           mnl::basics::matrixStack*	m_buf2;
                                           mnl::basics::matrixStack*	m_buf3;
                                           mnl::basics::matrixStack*	m_buf4;
                                           mnl::basics::matrixStack*	m_buf5;
                                           mnl::basics::Vector*		m_restricted;
                                           mnl::basics::Vector*		m_restricted1;
                                           bool 						m_mine;
                                           int							m_rank;
                                           int							m_size;
                                           mnl::basics::coarseGrid 	m_division;
                                           bool 						m_deflated;
                                           SEMLaplacianEvaluator::BC 	m_bc;

                                           std::vector<mnl::basics::coarseGrid> 		   m_desc;
                                           mnl::basics::Field2<mnl::basics::matrixStack>* m_W;
                                       };

class SEMTensoredConsistentPressurePreconditioner : 
  public mnl::utilities::Evaluator<mnl::basics::matrixStack> {
    public:
      SEMTensoredConsistentPressurePreconditioner(const mnl::basics::geometryStack& geometry,
          const mnl::basics::Vector& weight,
          const mnl::basics::Vector& weightGL,
          const mnl::basics::Vector& grid,
          const mnl::basics::Vector& gridGL,
          const mnl::basics::Matrix& GLL2G,
          const mnl::basics::Matrix& D,
          const Real nu, int rank=0, int size=1);
      virtual ~SEMTensoredConsistentPressurePreconditioner();

      virtual void evaluate(mnl::basics::matrixStack& res, const mnl::basics::matrixStack& u) const;

      static mnl::basics::Matrix constructInvMass(const mnl::basics::Vector& L,
          const mnl::basics::Vector& weight,
          SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
      static mnl::basics::Matrix constructIntMass(const mnl::basics::Vector& L,
          const mnl::basics::Vector& weightGL,
          const mnl::basics::Matrix& GLL2G,
          SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
      static mnl::basics::Matrix constructDivergence(const mnl::basics::Vector& weightGL,
          const mnl::basics::Matrix& GLL2G,
          const mnl::basics::Matrix& D, int elem,
          SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
      static mnl::basics::Matrix constructRestriction(int N, int elem);
      static mnl::basics::Matrix constructE0(const mnl::basics::Matrix& E, const mnl::basics::Matrix& R);
      static void constructEn(mnl::basics::Matrix& E, int N, int elem);

      std::vector<legendreLegendreW::poissonSolver*> m_SP;
      const mnl::basics::Vector& m_grid;
      const mnl::basics::geometryStack& m_geometry;
      mnl::basics::matrixStack* m_coarse;
      mnl::basics::coarseGrid m_desc;
      Real m_nu;
      int m_rank;
      int m_size;
  };

class SEMFEMConsistentPressurePreconditioner :
  public mnl::utilities::Evaluator<mnl::basics::matrixStack> {
    public:
      SEMFEMConsistentPressurePreconditioner(mnl::basics::geometryStack& geometry,
          const mnl::basics::Vector& weight,
          const mnl::basics::Vector& weightGL,
          const mnl::basics::Vector& grid,
          const mnl::basics::Vector& gridGL,
          const Real nu, int rank=0, int size=1);

      virtual ~SEMFEMConsistentPressurePreconditioner();

      virtual void evaluate(mnl::basics::matrixStack& res,
          const mnl::basics::matrixStack& p) const;

      int m_rank, m_size;
      Real m_nu;
      mnl::basics::geometryStack&    m_geometry;
      mnl::basics::matrixStack*	   m_coarse;
      mnl::basics::matrixStack*	   m_coarse2;
      mnl::basics::matrixStack*	   m_coarse3;
      linearElements::poissonSolver* m_SP;
      linearElements::poissonSolver* m_SP2;
      mnl::basics::Matrix*		   m_A;
      mnl::basics::Matrix* 		   m_restrict;
      mnl::basics::matrixStack*	   m_LG;
      mnl::basics::matrixStack*	   m_work;
      mnl::basics::Vector*		   m_restricted;
      mnl::basics::Vector*		   m_restricted1;
      std::vector<mnl::basics::coarseGrid> m_group;
      std::vector<Real>			   m_Lx;
      std::vector<Real>			   m_Ly;
      mnl::basics::Matrix			   m_GLL2G;
      mnl::basics::matrixStack*	   m_velTemp;
      mnl::basics::matrixStack*	   m_velBuf;
      mnl::basics::matrixStack*	   m_velBuf2;
      std::vector<mnl::basics::Matrix*> m_RT;
  };

class extrudedTensoredConsistentPressurePreconditioner : 
  public mnl::utilities::Evaluator<mnl::basics::matricesStack> {
    public:
      extrudedTensoredConsistentPressurePreconditioner(const mnl::basics::geometryStack& geometry,
          const mnl::basics::Vector& weight,
          const mnl::basics::Vector& weightGL,
          const mnl::basics::Vector& grid,
          const mnl::basics::Vector& gridGL,
          const mnl::basics::Matrix& GLL2G,
          const mnl::basics::Matrix& D,
          int rank=0, int size=1);
      virtual ~extrudedTensoredConsistentPressurePreconditioner();

      virtual void evaluate(mnl::basics::matricesStack& res, const mnl::basics::matricesStack& u) const;
    protected:
      std::vector<Real> m_LyLx;
      std::vector<legendreLegendreW::poissonSolver*> m_SP;
      const mnl::basics::Vector& m_grid;
      const mnl::basics::geometryStack& m_geometry;
      mnl::basics::Matrix* m_E;
      mnl::basics::matricesStack* m_coarse;
      mnl::basics::coarseGrid m_desc;
      int m_rank;
      int m_size;
  };


class extrudedConsistentPressureEvaluator : public 
                                            mnl::utilities::Evaluator<mnl::basics::matricesStack> {
                                              public:
                                                extrudedConsistentPressureEvaluator(mnl::basics::geometryStack& geometry,
                                                    const mnl::basics::Matrix& D,
                                                    const mnl::basics::Matrix& GLL2G,
                                                    const mnl::basics::Vector& weightGL,
                                                    const mnl::basics::Vector& weight,
                                                    const mnl::basics::Vector& grid,
                                                    const mnl::basics::Vector& gridGL,
                                                    bool preconditioned=true,
                                                    int rank=0, int size=1,
                                                    SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET,
                                                    int elems=1);
                                                virtual ~extrudedConsistentPressureEvaluator();

                                                virtual void evaluate(mnl::basics::matricesStack& res, 
                                                    const mnl::basics::matricesStack& p) const;

                                                virtual bool isL2() const
                                                {
                                                  return true;
                                                }

                                                void filter(mnl::basics::matricesStack& res) const
                                                {
                                                  applyFilter(res);
                                                }

                                                static void applyFilter(mnl::basics::matricesStack& res)
                                                {
                                                  Real sum = res.sum();
                                                  int length=res.length();
#ifdef HAS_MPI
                                                  Real mysum = sum;
                                                  int mylength=length;
                                                  MPI_Allreduce(&mysum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                                                  MPI_Allreduce(&mylength,&length,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
                                                  res -= sum/length;
                                                }

                                                void solve(mnl::basics::matricesStack& res);
                                                int solve3(mnl::basics::matricesStack& res, 
                                                    mnl::basics::matricesStack& u,
                                                    extrudedTensoredConsistentPressurePreconditioner* pre=NULL);

                                                legendreLegendreW::poissonSolver 			m_SP;
                                                SEMConsistentPressureEvaluator				m_eval;
                                                std::vector<mnl::utilities::Evaluator<mnl::basics::matrixStack>*> m_pre;
                                                std::vector<SEMConsistentPressureEvaluator*> m_evals;
                                                extrudedGradientEvaluator		 			m_gradient;
                                                extrudedDivergenceEvaluator	 	 			m_divergence;

                                                mnl::basics::Vector*						m_restricted;
                                                mnl::basics::Vector*						m_restricted1;
                                                std::vector<mnl::basics::coarseGrid> 		m_desc;
                                                mnl::basics::coarseGrid 					m_division;
                                                mnl::basics::Matrix* 						m_E;
                                                mnl::basics::matricesStack*					m_buf1;
                                                mnl::basics::matricesStack*					m_buf2;
                                                mnl::basics::Field3<mnl::basics::matricesStack>* m_W;
                                                extrudedTensoredConsistentPressurePreconditioner* m_pre3;
                                                bool m_deflated;
                                                int m_rank;
                                                int m_size;
                                                std::vector<int*> m_scount;
                                                std::vector<int*> m_sdispl;
                                                std::vector<int*> m_rcount;
                                                std::vector<int*> m_rdispl;
                                                std::pair< std::vector<int>, std::vector< std::vector<int> > > m_counts;
                                                SEMLaplacianEvaluator::BC m_bc;
                                                mnl::utilities::iterationStatistics m_stat;
                                              private:
                                                void evaluateSimple(mnl::basics::matricesStack& res, const mnl::basics::matricesStack& p) const;
                                            };

class extrudedFEMConsistentPressurePreconditioner :
  public mnl::utilities::Evaluator<mnl::basics::matricesStack> {
    public:
      extrudedFEMConsistentPressurePreconditioner(mnl::basics::geometryStack& geometry,
          const mnl::basics::Vector& weight,
          const mnl::basics::Vector& weightGL,
          const mnl::basics::Vector& grid,
          const mnl::basics::Vector& gridGL,
          int rank=0, int size=1);

      virtual ~extrudedFEMConsistentPressurePreconditioner();

      virtual void evaluate(mnl::basics::matricesStack& res,
          const mnl::basics::matricesStack& p) const;
      int m_rank, m_size;
      mnl::basics::geometryStack&    m_geometry;
      mnl::basics::matricesStack*	   m_coarse;
      mnl::basics::matricesStack*	   m_coarse2;
      mnl::basics::matricesStack*	   m_coarse3;
      linearElements::poissonSolver* m_SP;
      linearElements::poissonSolver* m_SP2;
      mnl::basics::Matrix*		   m_A;
      mnl::basics::Matrix* 		   m_restrict;
      mnl::basics::matricesStack*	   m_LG;
      mnl::basics::matricesStack*	   m_work;
      mnl::basics::matricesStack*	   m_work2;
      mnl::basics::matricesStack*	   m_work3;
      mnl::basics::Vector*		   m_restricted;
      mnl::basics::Vector*		   m_restricted1;
      std::vector<mnl::basics::coarseGrid> m_group;
      std::vector<Real>			   m_Lx;
      std::vector<Real>			   m_Ly;
      mnl::basics::Matrix			   m_GLL2G;
      mnl::basics::Matrix*		   m_velTemp;
      mnl::basics::Matrices*		   m_velTemp2;
      mnl::basics::matricesStack*	   m_velBuf;
      mnl::basics::matricesStack*	   m_velBuf2;
      std::vector<mnl::basics::Matrix*> m_RT;
  };

class extrudedConsistentPressurePreconditioner :
  public mnl::utilities::Evaluator<mnl::basics::matricesStack> {
    public:
      extrudedConsistentPressurePreconditioner(extrudedConsistentPressureEvaluator& eval,
          const mnl::basics::geometryStack3D& geometry) :
        m_eval(eval), m_geometry(geometry)
    {
    }

      virtual void evaluate(mnl::basics::matricesStack& res,
          const mnl::basics::matricesStack& p) const
      {
        mnl::basics::matricesStack* temp 
          = m_geometry.localToGlobalZ(p,m_eval.m_rank,m_eval.m_size,true);
        m_eval.solve(*temp);
        m_geometry.globalToLocalZ(res,temp,m_eval.m_rank,m_eval.m_size,true);
        Real sum = res.sum();
        int length=res.length();
#ifdef HAS_MPI
        Real mysum = sum;
        MPI_Allreduce(&mysum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        length *= m_eval.m_size;
#endif
        res -= sum/length;
      }
    protected:
      extrudedConsistentPressureEvaluator& m_eval;
      const mnl::basics::geometryStack3D& m_geometry;
  };

class SEM3DConsistentPressureEvaluator: public 
                                        mnl::utilities::Evaluator<mnl::basics::matricesStack> {
                                          public:
                                            SEM3DConsistentPressureEvaluator(mnl::basics::geometryStack3D& geometry,
                                                const mnl::basics::Matrix& D,
                                                const mnl::basics::Matrix& GLL2G,
                                                const mnl::basics::Vector& weightGL,
                                                const mnl::basics::Vector& weight,
                                                int rank=0, int size=1,
                                                SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
                                            virtual ~SEM3DConsistentPressureEvaluator();

                                            virtual void evaluate(mnl::basics::matricesStack& res, 
                                                const mnl::basics::matricesStack& p) const;

                                            virtual bool isL2() const
                                            {
                                              return true;
                                            }

                                            virtual void filter(mnl::basics::matricesStack& res) const
                                            {
                                              Real sum = res.sum();
                                              int length=res.length();
#ifdef HAS_MPI
                                              Real mysum = sum;
                                              MPI_Allreduce(&mysum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                                              length *= m_size;
#endif
                                              res -= sum/length;
                                            }

                                            SEM3DGradientEvaluator		 			m_gradient;
                                            SEM3DDivergenceEvaluator	 	 		m_divergence;

                                            mnl::basics::Field3<mnl::basics::matricesStack>* m_W;
                                            int m_rank;
                                            int m_size;
                                            SEMLaplacianEvaluator::BC m_bc;
                                            mnl::basics::coarseGrid 	m_division;
                                        };

class SEM3DFEMConsistentPressurePreconditioner :
  public mnl::utilities::Evaluator<mnl::basics::matricesStack> {
    public:
      SEM3DFEMConsistentPressurePreconditioner(mnl::basics::geometryStack3D& geometry,
          const mnl::basics::Vector& weight,
          const mnl::basics::Vector& weightGL,
          const mnl::basics::Vector& grid,
          const mnl::basics::Vector& gridGL,
          int rank=0, int size=1);

      virtual ~SEM3DFEMConsistentPressurePreconditioner();

      virtual void evaluate(mnl::basics::matricesStack& res,
          const mnl::basics::matricesStack& p) const;

      int m_rank, m_size;
      mnl::basics::geometryStack3D&  m_geometry;
      mnl::basics::matricesStack*	   m_coarse;
      mnl::basics::matricesStack*	   m_coarse2;
      mnl::basics::matricesStack*	   m_coarse3;
      std::vector<linearElements::poissonSolver*> m_SP;
      mnl::basics::Matrix*		   m_A;
      mnl::basics::Matrix* 		   m_restrict;
      mnl::basics::matricesStack*	   m_LG;
      mnl::basics::matricesStack*	   m_work;
      mnl::basics::matricesStack*	   m_work2;
      mnl::basics::Vector*		   m_restricted;
      mnl::basics::Vector*		   m_restricted1;
      std::vector<mnl::basics::coarseGrid> m_group;
      std::vector<Real>			   m_Lx;
      std::vector<Real>			   m_Ly;
      std::vector<Real>			   m_Lz;
      mnl::basics::Matrix			   m_GLL2G;
      mnl::basics::Matrix*		   m_velTemp;
      mnl::basics::Matrices*		   m_velTemp2;
      mnl::basics::matricesStack*	   m_velBuf;
      mnl::basics::matricesStack*	   m_velBuf2;
  };

#endif

