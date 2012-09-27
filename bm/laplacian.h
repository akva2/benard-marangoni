#ifndef LAPLACIAN_H_
#define LAPLACIAN_H_

#include "mnl/cgsolver.h"
#include "mnl/geometry.h"
#include "mnl/util.h"

#include "poissonsolver-lgw.h"
#include "legendrelegendrew.h"
#include "poissonsolver-fe.h"

class deformedLaplacianEvaluator : public mnl::utilities::Evaluator<mnl::basics::Matrix> {
  public:
    deformedLaplacianEvaluator(const mnl::basics::matrixStack& G,
        const mnl::basics::Matrix& D,
        const mnl::basics::Vector& weight,
        mnl::basics::Matrix& Uxi,
        mnl::basics::Matrix& Ueta,
        mnl::basics::Matrix& t1,
        mnl::basics::Matrix& t2,
        mnl::basics::Matrix& t3);

    void evaluate(mnl::basics::Matrix& res,
        const mnl::basics::Matrix& u, bool bMask) const;

    virtual void evaluate(mnl::basics::Matrix& res, const mnl::basics::Matrix& u) const
    {
      evaluate(res,u,true);
    }
  protected:
    const mnl::basics::Matrix& m_D;
    const mnl::basics::Vector& m_weight;
    const mnl::basics::matrixStack& m_G;
    mnl::basics::Matrix& m_Uxi;
    mnl::basics::Matrix& m_Ueta;
    mnl::basics::Matrix& m_t1;
    mnl::basics::Matrix& m_t2;
    mnl::basics::Matrix& m_t3;
};

class deformed3DLaplacianEvaluator : public mnl::utilities::Evaluator<mnl::basics::Matrices> {
  public:
    deformed3DLaplacianEvaluator(const mnl::basics::matricesStack& G,
        const mnl::basics::Matrix& D,
        const mnl::basics::Vector& weight,
        mnl::basics::Matrices& Uxi,
        mnl::basics::Matrices& Ueta,
        mnl::basics::Matrices& Ugamma,
        mnl::basics::Matrices& t1,
        mnl::basics::Matrices& t2,
        mnl::basics::Matrices& t3);

    void evaluate(mnl::basics::Matrices& res,
        const mnl::basics::Matrices& u, bool bMask) const;

    virtual void evaluate(mnl::basics::Matrices& res, const mnl::basics::Matrices& u) const
    {
      evaluate(res,u,true);
    }
  protected:
    const mnl::basics::matricesStack& m_G;
    const mnl::basics::Matrix& m_D;
    const mnl::basics::Vector& m_weight;
    mnl::basics::Matrices& m_Uxi;
    mnl::basics::Matrices& m_Ueta;
    mnl::basics::Matrices& m_Ugamma;
    mnl::basics::Matrices& m_t1;
    mnl::basics::Matrices& m_t2;
    mnl::basics::Matrices& m_t3;
};

class deformedTensoredLaplacianPreconditioner : 
  public mnl::utilities::Evaluator<mnl::basics::Matrix> {
    public:
      deformedTensoredLaplacianPreconditioner(const mnl::basics::spectralElement2D& G,
          const legendreLegendreW::poissonSolver& SP,
          const Real nu);

      virtual void evaluate(mnl::basics::Matrix& res, const mnl::basics::Matrix& u) const;
    protected:
      Real 									m_LyLx;
      Real									m_nu;
      const legendreLegendreW::poissonSolver& m_SP;
  };

class SEMLaplacianEvaluator : public 
                              mnl::utilities::Evaluator<mnl::basics::matrixStack> {
                                public:
                                  typedef enum BC {
                                    HOMOGENOUS_DIRICHLET,
                                    HOMOGENOUS_NEUMANN,
                                    DIRICHLET_DIRICHLET_BOTTOM,
                                    NEUMANN_DIRICHLET_BOTTOM,
                                    PERIODIC,
                                    PERIODIC_DIRICHLET,
                                    PERIODIC_DIRICHLET_BOTTOM
                                  } BC;

                                  SEMLaplacianEvaluator(mnl::basics::geometryStack& geometry,
                                      const mnl::basics::Matrix& D,	
                                      const mnl::basics::Vector& weight,
                                      const std::vector<mnl::basics::matrixStack*>* G=NULL,
                                      int rank=0, int size=1, BC bc=HOMOGENOUS_DIRICHLET);
                                  virtual ~SEMLaplacianEvaluator();

                                  void evaluate(mnl::basics::matrixStack& res, 
                                      const mnl::basics::matrixStack& u, bool bMask, bool doSum) const;

                                  virtual void evaluate(mnl::basics::matrixStack& res, 
                                      const mnl::basics::matrixStack& u) const
                                  {
                                    evaluate(res,u,m_bc==HOMOGENOUS_DIRICHLET?true:false,true);
                                  }

                                  inline void evaluate(mnl::basics::Field2<mnl::basics::matrixStack>& res,
                                      mnl::basics::Field2<mnl::basics::matrixStack>& u,
                                      bool bMask, bool doSum)
                                  {
                                    evaluate(res.X(),u.X(),bMask,doSum);
                                    evaluate(res.Y(),u.Y(),bMask,doSum);
                                  }

                                  Real m_nu;
                                  mnl::basics::geometryStack& m_geometry;
                                  const mnl::basics::Matrix& m_D;
                                  const mnl::basics::Vector& m_weight;
                                  std::vector<mnl::basics::matrixStack*> m_G;
                                  std::vector<deformedLaplacianEvaluator*> m_evals;
                                  mnl::basics::componentView<mnl::basics::Matrix,mnl::basics::matrixStack>* m_view;
                                  mnl::basics::coarseGrid m_division;
                                  mnl::basics::matrixStack* m_buffer;
                                  mnl::basics::matrixStack* m_nosum;
                                  bool m_mine;
                                  int m_rank; 
                                  int m_tag;
                                  int m_size;
                                  BC m_bc;
                                  mnl::basics::matrixStack*   m_Uxi;
                                  mnl::basics::matrixStack*   m_Ueta;
                                  mnl::basics::matrixStack*   m_t1;
                                  mnl::basics::matrixStack*   m_t2;
                                  mnl::basics::matrixStack*   m_t3;
                              };

class SEMFEMLaplacianPreconditioner : public
                                      mnl::utilities::Evaluator<mnl::basics::matrixStack> {
                                        public:
                                          SEMFEMLaplacianPreconditioner(mnl::basics::geometryStack& geometry,
                                              const mnl::basics::Vector& grid,
                                              const mnl::basics::Vector& weight,
                                              const Real nu,
                                              int rank=0, int size=1,
                                              SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
                                          virtual ~SEMFEMLaplacianPreconditioner();

                                          virtual void evaluate(mnl::basics::matrixStack& res,
                                              const mnl::basics::matrixStack& u) const;
                                          virtual mnl::basics::matrixStack* buffer()
                                          {
                                            return( m_nosum );
                                          }

                                          static void setupRestrictionOperators(std::vector<mnl::basics::Matrix*>& result,
                                              const mnl::basics::geometryStack& geometry,
                                              const mnl::basics::Vector& grid,
                                              int size1, int size2, 
                                              int size3, int size4,
                                              SEMLaplacianEvaluator::BC bc);

                                          Real m_nu;
                                          mnl::basics::matrixStack*	  m_nosum;
                                          int							  m_tag;
                                          SEMLaplacianEvaluator::BC	  m_bc;
                                          mnl::basics::geometryStack&   m_geometry;
                                          mnl::basics::matrixStack*	  m_coarse;
                                          mnl::basics::matrixStack*	  m_coarse2;
                                          mnl::basics::matrixStack*	  m_LG;
                                          mnl::basics::matrixStack*	  m_work;
                                          mnl::basics::Matrix* 		  m_restrict;
                                          mnl::basics::Vector*		  m_restricted;
                                          mnl::basics::Vector*		  m_restricted1;
                                          std::vector<mnl::basics::coarseGrid> m_group;
                                          std::vector<Real>			   m_Lx;
                                          std::vector<Real>			   m_Ly;
                                          linearElements::poissonSolver* m_SP;
                                          linearElements::poissonSolver* m_SP2;
                                          mnl::basics::Matrix*		   m_A;
                                          int							   m_rank;
                                          int							   m_size;
                                          std::vector<mnl::basics::Matrix*> m_RT;
                                      };

class SEMLaplacianDiagonalPreconditioner : public mnl::utilities::Evaluator<mnl::basics::matrixStack> {
  public:
    SEMLaplacianDiagonalPreconditioner(const mnl::basics::geometryStack& geometry, 
        const Real nu, const legendreLegendreW::poissonSolver& SP);

    void evaluate(mnl::basics::matrixStack& res, const mnl::basics::matrixStack& u) const;
  protected:
    mnl::basics::matrixStack m_iHD;
    const mnl::basics::geometryStack& m_geometry;
    Real m_nu;
};

class SEMTensoredLaplacianPreconditioner : 
  public mnl::utilities::Evaluator<mnl::basics::matrixStack> {
    public:
      SEMTensoredLaplacianPreconditioner(const mnl::basics::geometryStack& geometry,
          const Real nu, const legendreLegendreW::poissonSolver& SP);

      virtual ~SEMTensoredLaplacianPreconditioner();

      virtual void evaluate(mnl::basics::matrixStack& res, const mnl::basics::matrixStack& u) const;
    protected:
      std::vector<deformedTensoredLaplacianPreconditioner*> m_pre;
      SEMLaplacianDiagonalPreconditioner 					  m_eval;
  };



class extrudedLaplacianEvaluator : public 
                                   mnl::utilities::Evaluator<mnl::basics::matricesStack> {
                                     public:
                                       extrudedLaplacianEvaluator(SEMLaplacianEvaluator& eval, 
                                           const legendreLegendreW::poissonSolver& SP,
                                           const Real nu=0,
                                           const mnl::basics::Matrix* Az=NULL,
                                           bool preconditioned=true,
                                           int  rank=0, int size=1,
                                           SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
                                       virtual ~extrudedLaplacianEvaluator();

                                       virtual void evaluate(mnl::basics::matricesStack& res, 
                                           const mnl::basics::matricesStack& u) const
                                       {
                                         evaluate(res,u,true,true);
                                       }
                                       void evaluate(mnl::basics::matricesStack& res, const mnl::basics::matricesStack& u, bool bMask, bool doSum) const;
                                       inline void evaluate(mnl::basics::Field3<mnl::basics::matricesStack>& res,
                                           mnl::basics::Field3<mnl::basics::matricesStack>& u,
                                           bool bMask, bool doSum)
                                       {
                                         evaluate(res.X(),u.X(),bMask,doSum);
                                         evaluate(res.Y(),u.Y(),bMask,doSum);
                                         evaluate(res.Z(),u.Z(),bMask,doSum);
                                       }

                                       void solve(mnl::basics::matricesStack& res, mnl::basics::matricesStack& u, int comp=0);
                                       inline void solve(mnl::basics::Field3<mnl::basics::matricesStack>& res,
                                           mnl::basics::Field3<mnl::basics::matricesStack>& u)
                                       {
                                         std::string name = m_name;
                                         m_name = name+" X";
                                         solve(res.X(),u.X(),0);
                                         m_name = name+" Y";
                                         solve(res.Y(),u.Y(),1);
                                         m_name = name+" Z";
                                         solve(res.Z(),u.Z(),2);
                                         m_name = name;
                                       }

                                       void solve(mnl::basics::Field3<mnl::basics::matricesStack>& res,
                                           mnl::basics::Field3<mnl::basics::matricesStack>& u,
                                           extrudedLaplacianEvaluator& X,
                                           extrudedLaplacianEvaluator& Y,
                                           extrudedLaplacianEvaluator& Z);

                                       std::string m_name;
                                       SEMLaplacianEvaluator& m_eval;
                                       const legendreLegendreW::poissonSolver& m_SP;
                                       std::vector<mnl::utilities::Evaluator<mnl::basics::matrixStack>*> m_pre;
                                       std::vector<SEMLaplacianEvaluator*> m_evals;
                                       mnl::basics::coarseGrid m_division;
                                       Real m_nu;
                                       const mnl::basics::Matrix* m_Az;
                                       mnl::basics::matricesStack* m_nosum;
                                       int m_rank;
                                       int m_size;
                                       std::vector<int*> m_scount;
                                       std::vector<int*> m_sdispl;
                                       std::vector<int*> m_rcount;
                                       std::vector<int*> m_rdispl;
                                       SEMLaplacianEvaluator::BC m_bc;
                                       std::vector<mnl::utilities::iterationStatistics*> m_stat;
                                       std::pair< std::vector<int>, std::vector< std::vector<int> > > m_counts;
                                       mnl::utilities::iterationStatistics* m_stat3;
                                       std::vector<int*> m_scount3;
                                       std::vector<int*> m_sdispl3;
                                       std::vector<int*> m_rcount3;
                                       std::vector<int*> m_rdispl3;
                                       std::pair< std::vector<int>, std::vector< std::vector<int> > > m_counts3;
                                   };

class extrudedFEMLaplacianPreconditioner :
  public mnl::utilities::Evaluator<mnl::basics::matricesStack> {
    public:
      extrudedFEMLaplacianPreconditioner(mnl::basics::geometryStack& geometry,
          const mnl::basics::Vector& weight,
          const mnl::basics::Vector& grid,
          Real nu=0, int rank=0, int size=1,
          SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);

      virtual ~extrudedFEMLaplacianPreconditioner();

      virtual void evaluate(mnl::basics::matricesStack& res,
          const mnl::basics::matricesStack& u) const;

      virtual mnl::basics::matricesStack* buffer()
      {
        return( m_nosum );
      }

      Real m_nu;
      int m_rank, m_size;
      mnl::basics::geometryStack&    m_geometry;
      mnl::basics::matricesStack*	   m_coarse;
      mnl::basics::matricesStack*	   m_coarse2;
      linearElements::poissonSolver* m_SP;
      linearElements::poissonSolver* m_SP2;
      mnl::basics::Matrix*		   m_A;
      mnl::basics::Matrix* 		   m_restrict;
      mnl::basics::matricesStack*	   m_LG;
      mnl::basics::matricesStack*	   m_work;
      mnl::basics::matricesStack*	   m_work2;
      mnl::basics::matricesStack*	   m_work3;
      mnl::basics::Vector* 		   m_restricted;
      mnl::basics::Vector*		  m_restricted1;
      std::vector<mnl::basics::coarseGrid> m_group;
      std::vector<Real>			   m_Lx;
      std::vector<Real>			   m_Ly;
      mnl::basics::matricesStack*	   m_nosum;
      SEMLaplacianEvaluator::BC	   m_bc;
      std::vector<mnl::basics::Matrix*> m_RT;
  };

class extrudedLaplacianPreconditioner :
  public mnl::utilities::Evaluator<mnl::basics::matricesStack> {
    public:
      extrudedLaplacianPreconditioner(extrudedLaplacianEvaluator& eval,
          const mnl::basics::geometryStack3D& geometry);

      virtual void evaluate(mnl::basics::matricesStack& res,
          const mnl::basics::matricesStack& u) const;

      virtual mnl::basics::matricesStack* buffer()
      {
        return( m_eval.m_nosum );
      }
      extrudedLaplacianEvaluator& m_eval;
      const mnl::basics::geometryStack3D& m_geometry;
  };

class SEM3DLaplacianEvaluator :
  public mnl::utilities::Evaluator<mnl::basics::matricesStack> {
    public:
      SEM3DLaplacianEvaluator(const mnl::basics::geometryStack3D& geometry,
          const mnl::basics::Vector& weight,
          const mnl::basics::Matrix& D,
          Real nu=0, int rank=0, int size=1,
          SEMLaplacianEvaluator::BC bc=SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET);
      virtual ~SEM3DLaplacianEvaluator();

      void evaluate(mnl::basics::Field3<mnl::basics::matricesStack>& res,
          const mnl::basics::Field3<mnl::basics::matricesStack>& u,
          bool mask, bool doSum) const
      {
        evaluate(res.X(),u.X(),mask,doSum);
        evaluate(res.Y(),u.Y(),mask,doSum);
        evaluate(res.Z(),u.Z(),mask,doSum);
      }
      void evaluate(mnl::basics::matricesStack& res,
          const mnl::basics::matricesStack& u,
          bool mask, bool doSum) const;

      virtual void evaluate(mnl::basics::matricesStack& res,
          const mnl::basics::matricesStack& u) const
      {
        evaluate(res,u,true,true);
      }

      mnl::basics::matricesStack* m_nosum;
      Real m_nu;
      SEMLaplacianEvaluator::BC m_bc;

      int m_rank;
      int m_size;
      std::vector<mnl::basics::matricesStack*> m_G;
      std::vector<deformed3DLaplacianEvaluator*> m_evals;
      mnl::basics::componentView<mnl::basics::Matrices,mnl::basics::matricesStack>* m_view;
      const mnl::basics::geometryStack3D& m_geometry;
      const mnl::basics::Matrix&	  m_D;
      const mnl::basics::Vector&	  m_weight;
      mnl::basics::matricesStack*   m_Uxi;
      mnl::basics::matricesStack*   m_Ueta;
      mnl::basics::matricesStack*   m_Ugamma;
      mnl::basics::matricesStack*   m_t1;
      mnl::basics::matricesStack*   m_t2;
      mnl::basics::matricesStack*   m_t3;
      mnl::basics::coarseGrid m_division;
  };

class SEM3DFEMLaplacianPreconditioner :
  public mnl::utilities::Evaluator<mnl::basics::matricesStack> {
    public:
      SEM3DFEMLaplacianPreconditioner(mnl::basics::geometryStack3D& geometry,
          const mnl::basics::Vector& weight,
          const mnl::basics::Vector& grid,
          Real nu=0, int rank=0, int size=1);

      virtual ~SEM3DFEMLaplacianPreconditioner();

      virtual void evaluate(mnl::basics::matricesStack& res,
          const mnl::basics::matricesStack& u) const;

      virtual mnl::basics::matricesStack* buffer()
      {
        return( m_nosum );
      }

      Real m_nu;
      int m_rank, m_size;
      mnl::basics::geometryStack3D&  m_geometry;
      mnl::basics::matricesStack*	   m_coarse;
      mnl::basics::matricesStack*	   m_coarse2;
      std::vector<linearElements::poissonSolver*> m_SP;
      mnl::basics::Matrix*		   m_A;
      mnl::basics::Matrix* 		   m_restrict;
      mnl::basics::matricesStack*	   m_LG;
      mnl::basics::matricesStack*	   m_work;
      mnl::basics::matricesStack*	   m_work2;
      mnl::basics::Vector*		   m_restricted;
      std::vector<mnl::basics::coarseGrid> m_group;
      std::vector<Real>			   m_Lx;
      std::vector<Real>			   m_Ly;
      std::vector<Real>			   m_Lz;
      mnl::basics::matricesStack*	   m_nosum;
      SEMLaplacianEvaluator::BC	   m_bc;
  };

#endif

