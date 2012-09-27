#ifndef GRADIENT_H_
#define GRADIENT_H_

#include "laplacian.h"

#include <sstream>

class deformedGradientEvaluator : public 
                                  mnl::utilities::scalarToFieldEvaluator<mnl::basics::Matrix, mnl::basics::Field2> {
                                    public:
                                      deformedGradientEvaluator(const mnl::basics::matrixStack& G,
                                          const mnl::basics::Matrix& Dt,
                                          const mnl::basics::Matrix& GLL2G,
                                          const mnl::basics::Vector& weight,
                                          mnl::basics::Matrix& W,
                                          mnl::basics::Matrix& t,
                                          mnl::basics::Matrix& s);

                                      virtual void evaluate(mnl::basics::Field2<mnl::basics::Matrix>& res, const mnl::basics::Matrix& p) const;
                                    protected:
                                      mnl::basics::Matrix& m_W;
                                      mnl::basics::Matrix& m_t;
                                      mnl::basics::Matrix& m_s;

                                      const mnl::basics::matrixStack& m_G;
                                      const mnl::basics::Matrix& m_GLL2G;
                                      const mnl::basics::Matrix& m_Dt;
                                      const mnl::basics::Vector& m_weight;
                                  };

class deformed3DGradientEvaluator : public 
                                    mnl::utilities::scalarToFieldEvaluator<mnl::basics::Matrices, mnl::basics::Field3> {
                                      public:
                                        deformed3DGradientEvaluator(const mnl::basics::matricesStack& G,
                                            const mnl::basics::Matrix& Dt,
                                            const mnl::basics::Matrix& GLL2G,
                                            const mnl::basics::Vector& weight,
                                            mnl::basics::Matrices& W,
                                            mnl::basics::Matrices& t,
                                            mnl::basics::Matrix& s,
                                            mnl::basics::Matrices& v);

                                        virtual void evaluate(mnl::basics::Field3<mnl::basics::Matrices>& res, const mnl::basics::Matrices& p) const;
                                      protected:
                                        mnl::basics::Matrices& m_W;
                                        mnl::basics::Matrices& m_t;
                                        mnl::basics::Matrix& m_s;
                                        mnl::basics::Matrices& m_v;

                                        const mnl::basics::matricesStack& m_G;
                                        const mnl::basics::Matrix& m_GLL2G;
                                        const mnl::basics::Matrix& m_Dt;
                                        const mnl::basics::Vector& m_weight;
                                      private:
                                        void evaluateComponent(mnl::basics::Matrices& res, int index) const;
                                    };

class SEMGradientEvaluator : public
                             mnl::utilities::scalarToFieldEvaluator<mnl::basics::matrixStack,mnl::basics::Field2> {
                               public:
                                 friend class extrudedGradientEvaluator;

                                 SEMGradientEvaluator(const mnl::basics::geometryStack& geometry,
                                     const mnl::basics::Matrix& D,
                                     const mnl::basics::Matrix& GLL2G,
                                     const mnl::basics::Vector& weightGL,
                                     const std::vector<mnl::basics::matrixStack*>* G = NULL,
                                     int rank=0, int size=1);
                                 virtual ~SEMGradientEvaluator();

                                 virtual void evaluate(mnl::basics::Field2<mnl::basics::matrixStack>& res,
                                     const mnl::basics::matrixStack& u) const;

                                 const mnl::basics::Matrix& m_GLL2G;
                                 std::vector<mnl::basics::matrixStack*> m_G;
                               protected:
                                 bool m_mine;
                                 std::vector<deformedGradientEvaluator*> m_evals;
                                 mnl::basics::componentView<mnl::basics::Matrix,mnl::basics::matrixStack>* m_view;
                                 mnl::basics::Matrix m_Dt;
                                 const mnl::basics::Matrix& m_D;
                                 const mnl::basics::Vector& m_weight;
                                 const mnl::basics::geometryStack& m_geometry;
                                 mnl::basics::coarseGrid m_division;
                                 int m_rank;
                                 int m_size;
                                 mnl::basics::matrixStack* m_W;
                                 mnl::basics::matrixStack* m_t;
                                 mnl::basics::matrixStack* m_s;
                             };

class extrudedGradientEvaluator : public
                                  mnl::utilities::scalarToFieldEvaluator<mnl::basics::matricesStack,mnl::basics::Field3> {
                                    public:
                                      friend class SEMGradientEvaluator;

                                      extrudedGradientEvaluator(const SEMGradientEvaluator& eval, int rank=0, int size=1);
                                      virtual ~extrudedGradientEvaluator();

                                      virtual void evaluate(mnl::basics::Field3<mnl::basics::matricesStack>& res, 
                                          const mnl::basics::matricesStack& p) const;
                                    protected:
                                      const SEMGradientEvaluator& m_eval;
                                      int m_rank;
                                      int m_size;
                                      mnl::basics::coarseGrid m_division;
                                  };

class SEM3DGradientEvaluator : public
                               mnl::utilities::scalarToFieldEvaluator<mnl::basics::matricesStack,mnl::basics::Field3> {
                                 public:
                                   SEM3DGradientEvaluator(const mnl::basics::geometryStack3D& geometry,
                                       const mnl::basics::Matrix& D,
                                       const mnl::basics::Matrix& GLL2G,
                                       const mnl::basics::Vector& weightGL,
                                       const std::vector<mnl::basics::matricesStack*>* G = NULL,
                                       int rank=0, int size=1);
                                   virtual ~SEM3DGradientEvaluator();

                                   virtual void evaluate(mnl::basics::Field3<mnl::basics::matricesStack>& res,
                                       const mnl::basics::matricesStack& p) const;

                                   const mnl::basics::Matrix& m_GLL2G;
                                   std::vector<mnl::basics::matricesStack*> m_G;
                                 protected:
                                   bool m_mine;
                                   std::vector<deformed3DGradientEvaluator*> m_evals;
                                   mnl::basics::componentView<mnl::basics::Matrices,mnl::basics::matricesStack>* m_view;
                                   mnl::basics::Matrix m_Dt;
                                   const mnl::basics::Matrix& m_D;
                                   const mnl::basics::Vector& m_weight;
                                   const mnl::basics::geometryStack3D& m_geometry;
                                   mnl::basics::coarseGrid m_division;
                                   int m_rank;
                                   int m_size;
                                   mnl::basics::matricesStack* m_W;
                                   mnl::basics::matricesStack* m_t;
                                   mnl::basics::matrixStack* m_s;
                                   mnl::basics::matricesStack* m_v;
                               };

#endif

