#ifndef MASS_H_
#define MASS_H_

#include "mnl/geometry.h"

class SEMMassEvaluator : public mnl::utilities::Evaluator<mnl::basics::matrixStack> {
  public:
    friend class extrudedMassEvaluator;
    friend class extrudedInverseMassEvaluator;

    SEMMassEvaluator(const mnl::basics::geometryStack& geometry,
        const mnl::basics::Vector& weight,
        const mnl::basics::Matrix& GLL2G,
        int rank=0, int size=1);
    virtual ~SEMMassEvaluator();

    virtual void evaluate(mnl::basics::matrixStack& res, const mnl::basics::matrixStack& u) const;
    const mnl::basics::matrixStack& getJacobian() const
    {
      return( *m_J );
    }
    int m_rank;
    int m_size;
    const mnl::basics::geometryStack& m_geometry;
  protected:
    mnl::basics::matrixStack* 		  m_J;
    mnl::basics::coarseGrid 		  m_division;
    const mnl::basics::Vector&        m_weight;
};

class SEMInverseMassEvaluator : public SEMMassEvaluator {
  public:
    friend class extrudedInverseMassEvaluator;

    SEMInverseMassEvaluator(const mnl::basics::geometryStack& geometry,
        const mnl::basics::Vector& weight,
        const mnl::basics::Matrix& GLL2G,
        int rank=0, int size=1);
    virtual ~SEMInverseMassEvaluator();

    virtual void evaluate(mnl::basics::matrixStack& res, const mnl::basics::matrixStack& p) const;
    virtual void evaluateMass(mnl::basics::matrixStack& res, const mnl::basics::matrixStack& p) const;

    const mnl::basics::matrixStack& getInverseJacobian() const
    {
      return( *m_iJ );
    }
  protected:
    mnl::basics::matrixStack* m_iJ;
};

class SEM3DMassEvaluator : public mnl::utilities::Evaluator<mnl::basics::matricesStack> {
  public:
    SEM3DMassEvaluator(const mnl::basics::geometryStack3D& geometry,
        const mnl::basics::Vector& weight,
        const mnl::basics::Matrix& GLL2G,
        int rank=0, int size=1);
    virtual ~SEM3DMassEvaluator();

    virtual void evaluate(mnl::basics::matricesStack& res, const mnl::basics::matricesStack& u) const;
    const mnl::basics::matricesStack& getJacobian() const
    {
      return( *m_J );
    }
    int m_rank;
    int m_size;
    const mnl::basics::geometryStack3D& m_geometry;
  protected:
    mnl::basics::matricesStack* 	  m_J;
    mnl::basics::coarseGrid 		  m_division;
    const mnl::basics::Vector&        m_weight;
};

class SEM3DInverseMassEvaluator : public SEM3DMassEvaluator {
  public:
    SEM3DInverseMassEvaluator(const mnl::basics::geometryStack3D& geometry,
        const mnl::basics::Vector& weight,
        const mnl::basics::Matrix& GLL2G,
        int rank=0, int size=1);
    virtual ~SEM3DInverseMassEvaluator();

    virtual void evaluate(mnl::basics::matricesStack& res, const mnl::basics::matricesStack& p) const;
    virtual void evaluateMass(mnl::basics::matricesStack& res, const mnl::basics::matricesStack& p) const;

    const mnl::basics::matricesStack& getInverseJacobian() const
    {
      return( *m_iJ );
    }
  protected:
    mnl::basics::matricesStack* m_iJ;
};

class extrudedMassEvaluator : public
                              mnl::utilities::Evaluator<mnl::basics::matricesStack> {
                                public:
                                  extrudedMassEvaluator(const SEMMassEvaluator& eval, int rank=0, int size=1);

                                  virtual void evaluate(mnl::basics::matricesStack& res, const mnl::basics::matricesStack& p) const;

                                  const SEMMassEvaluator& m_eval;
                                  int m_rank;
                                  int m_size;
                                  mnl::basics::coarseGrid m_division;
                              };

class extrudedInverseMassEvaluator : public 
                                     extrudedMassEvaluator {
                                       public:
                                         extrudedInverseMassEvaluator(const SEMInverseMassEvaluator& eval, int rank=0, int size=1);

                                         virtual void evaluate(mnl::basics::matricesStack& res, const mnl::basics::matricesStack& p) const;
                                         virtual void evaluateMass(mnl::basics::matricesStack& res, const mnl::basics::matricesStack& p) const;

                                         const SEMInverseMassEvaluator& m_invEval;
                                     };

#endif

