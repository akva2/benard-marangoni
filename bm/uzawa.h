#ifndef UZAWA_H_
#define UZAWA_H_

#include "divergence.h"
#include "laplacian.h"
#include "gradient.h"
#include "consistent.h"

#include <limits>

class SEMCahouetChabardPreconditioner : public mnl::utilities::Evaluator<mnl::basics::matrixStack> {
  public:
    SEMCahouetChabardPreconditioner(mnl::basics::geometryStack& geometry,
        const mnl::basics::Matrix& D,
        const mnl::basics::Matrix& GLL2G,
        const mnl::basics::Vector& weightGL,
        const mnl::basics::Vector& weight,
        const mnl::basics::Vector& grid,
        const mnl::basics::Vector& gridGL);

    virtual void evaluate(mnl::basics::matrixStack& res, const mnl::basics::matrixStack&p) const;

    Real m_Dt;
    SEMInverseMassEvaluator		   m_invMass;
  protected:
    SEMConsistentPressureEvaluator m_eval;
    SEMTensoredConsistentPressurePreconditioner m_pre;
};

class SEMUzawaEvaluator : public mnl::utilities::Evaluator<mnl::basics::matrixStack> {
  public:
    SEMUzawaEvaluator(mnl::basics::geometryStack& geometry,
        const mnl::basics::Matrix& D,
        const mnl::basics::Matrix& GLL2G,
        const mnl::basics::Vector& grid_gll,
        const mnl::basics::Vector& weight_gll,
        const mnl::basics::Vector& weight_gl,
        SEMLaplacianEvaluator& laplacian,
        bool preconditioned=true,
        int rank=0, int size=1);
    virtual ~SEMUzawaEvaluator();

    virtual void evaluate(mnl::basics::matrixStack& res, const mnl::basics::matrixStack& p) const;

    void source(mnl::basics::Field2<mnl::basics::matrixStack>& u, 
        const mnl::basics::Vector& grid, const mnl::basics::function2D& ux, 
        const mnl::basics::function2D& uy, const mnl::basics::function2D& p, Real t) const;

    mnl::basics::coarseGrid m_division;
    mnl::basics::Field2<mnl::basics::matrixStack>* 			  m_W;

    std::vector<mnl::basics::matrixStack*>	m_G;
    SEMDivergenceEvaluator 				   	m_divergence;
    SEMLaplacianEvaluator& 			  	   	m_laplacian;
    SEMGradientEvaluator   			  	   	m_gradient;
    SEMFEMLaplacianPreconditioner*		   	m_pre;
    const mnl::basics::geometryStack& 	   	m_geometry;
    int									   	m_rank;
    int									   	m_size;
};

class extrudedCahouetChabardPreconditioner : public mnl::utilities::Evaluator<mnl::basics::matricesStack> {
  public:
    extrudedCahouetChabardPreconditioner(const extrudedInverseMassEvaluator& invMass,
        extrudedConsistentPressureEvaluator& consistent,
        Real nu);

    virtual void evaluate(mnl::basics::matricesStack& res,
        const mnl::basics::matricesStack& p) const;

    Real m_nu;
    const extrudedInverseMassEvaluator& m_invMass;
    extrudedConsistentPressureEvaluator& m_consistent;
};

class extrudedUzawaEvaluator : public mnl::utilities::Evaluator<mnl::basics::matricesStack> {
  public:
    extrudedUzawaEvaluator(mnl::basics::geometryStack& geometry,
        const mnl::basics::Matrix& GLL2G,
        const mnl::basics::Vector& weight_gl,
        const legendreLegendreW::poissonSolver& SP,
        extrudedLaplacianEvaluator& laplacian,
        int rank=0, int size=1,
        extrudedLaplacianEvaluator* laplacian2=NULL);

    virtual void evaluate(mnl::basics::matricesStack& res, const mnl::basics::matricesStack& p) const;

    void source(mnl::basics::Field3<mnl::basics::matricesStack>& u, 
        const mnl::basics::Vector&     grid, const mnl::basics::function3D& ux, 
        const mnl::basics::function3D& uy,   const mnl::basics::function3D& uz, 
        const mnl::basics::function3D& p,    Real t, int elem);
    void source(mnl::basics::Field3<mnl::basics::matricesStack>& u, 
        const mnl::basics::Vector&     grid, const mnl::basics::function3D& ux, 
        const mnl::basics::function3D& uy,   const mnl::basics::function3D& uz, 
        const mnl::basics::function3D& p,    Real t, int i, int elem);

    virtual bool isL2() const
    {
      return true;
    }

    void filter(mnl::basics::matricesStack& res) const
    {
      res -= res.sum()/res.length();
    }

    const mnl::basics::geometryStack&  	m_geometry;
    SEMUzawaEvaluator		   			m_uzawa2D;
    extrudedGradientEvaluator			m_gradient;
    extrudedLaplacianEvaluator& 		m_laplacian;
    extrudedLaplacianEvaluator* 		m_laplacian2;
    extrudedDivergenceEvaluator 		m_divergence;
    int									m_rank;
    int									m_size;
};

#endif

