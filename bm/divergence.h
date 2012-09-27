#ifndef DIVERGENCE_H_
#define DIVERGENCE_H_

#include "mnl/geometry.h"
#include "mnl/buffers.h"

class deformedDivergenceEvaluator : public 
                                    mnl::utilities::fieldToScalarEvaluator<mnl::basics::Matrix,mnl::basics::Field2> {
                                      public:
                                        deformedDivergenceEvaluator(const mnl::basics::matrixStack& G,
                                            const mnl::basics::Matrix& Dt,
                                            const mnl::basics::Matrix& GLL2G,
                                            const mnl::basics::Vector& weightGL,
                                            mnl::basics::Matrix& W,
                                            mnl::basics::Matrix& W1,
                                            mnl::basics::Matrix& buf);

                                        void evaluate(mnl::basics::Matrix& res,
                                            const mnl::basics::Field2<mnl::basics::Matrix>& u, bool mass) const;
                                        virtual void evaluate(mnl::basics::Matrix& res,
                                            const mnl::basics::Field2<mnl::basics::Matrix>& u) const
                                        {
                                          evaluate(res,u,true);
                                        }
                                      protected:
                                        const mnl::basics::matrixStack& m_G;
                                        const mnl::basics::Matrix& m_GLL2G;
                                        const mnl::basics::Matrix& m_Dt;
                                        const mnl::basics::Vector& m_weight;
                                        mnl::basics::Matrix& m_W;
                                        mnl::basics::Matrix& m_Wt;
                                        mnl::basics::Matrix& m_buf;
                                    };

class deformed3DDivergenceEvaluator : public
                                      mnl::utilities::fieldToScalarEvaluator<mnl::basics::Matrices,mnl::basics::Field3> {
                                        public:
                                          deformed3DDivergenceEvaluator(const mnl::basics::matricesStack& G,
                                              const mnl::basics::Matrix& Dt,
                                              const mnl::basics::Matrix& GLL2G,
                                              const mnl::basics::Vector& weightGL,
                                              mnl::basics::Matrix& W,
                                              mnl::basics::Matrices& W1,
                                              mnl::basics::Matrices& buf);

                                          virtual ~deformed3DDivergenceEvaluator();

                                          virtual void evaluate(mnl::basics::Matrices& res, const mnl::basics::Field3<mnl::basics::Matrices>& u) const;
                                        protected:
                                          const mnl::basics::matricesStack& m_G;
                                          const mnl::basics::Matrix& m_GLL2G;
                                          const mnl::basics::Matrix& m_Dt;
                                          const mnl::basics::Vector& m_weight;
                                          mnl::basics::Matrix& m_W;
                                          mnl::basics::Matrices& m_Wt;
                                          mnl::basics::Matrices& m_buf;
                                        private:
                                          void evaluateComponent(mnl::basics::Matrices& res,
                                              const mnl::basics::Matrices& field,
                                              int index) const;
                                      };

class SEMDivergenceEvaluator : public
                               mnl::utilities::fieldToScalarEvaluator<mnl::basics::matrixStack,mnl::basics::Field2> {
                                 public:
                                   friend class extrudedDivergenceEvaluator;

                                   SEMDivergenceEvaluator(const mnl::basics::geometryStack& geometry,
                                       const mnl::basics::Matrix& D,
                                       const mnl::basics::Matrix& GLL2G,
                                       const mnl::basics::Vector& weight,
                                       const std::vector<mnl::basics::matrixStack*>* G = NULL,
                                       int rank=0, int size=1);
                                   virtual ~SEMDivergenceEvaluator();

                                   virtual void evaluate(mnl::basics::matrixStack& res, const mnl::basics::Field2<mnl::basics::matrixStack>& u) const;
                                 protected:
                                   std::vector<mnl::basics::matrixStack*> m_G;
                                   std::vector<deformedDivergenceEvaluator*> m_evals;
                                   mnl::basics::componentView<mnl::basics::Matrix,mnl::basics::matrixStack>* m_view;
                                   const mnl::basics::geometryStack& m_geometry;
                                   const mnl::basics::Matrix& m_GLL2G;
                                   mnl::basics::Matrix m_Dt;
                                   const mnl::basics::Vector& m_weight;
                                   bool m_mine;
                                   int m_rank;
                                   int m_size;
                                   mnl::basics::coarseGrid m_division;
                                   mnl::basics::matrixStack* m_W;
                                   mnl::basics::matrixStack* m_Wt;
                                   mnl::basics::matrixStack* m_buf;
                               };

class extrudedDivergenceEvaluator : public 
                                    mnl::utilities::fieldToScalarEvaluator<mnl::basics::matricesStack,mnl::basics::Field3> {
                                      public:
                                        friend class SEMDivergenceEvaluator;

                                        extrudedDivergenceEvaluator(const SEMDivergenceEvaluator& eval, int rank=0, int size=1);
                                        virtual ~extrudedDivergenceEvaluator();

                                        virtual void evaluate(mnl::basics::matricesStack& res, 
                                            const mnl::basics::Field3<mnl::basics::matricesStack>& u) const;
                                      protected:
                                        const SEMDivergenceEvaluator& m_eval;
                                        int m_rank;
                                        int m_size;
                                        mnl::basics::coarseGrid m_division;
                                    };

class SEM3DDivergenceEvaluator : public
                                 mnl::utilities::fieldToScalarEvaluator<mnl::basics::matricesStack,mnl::basics::Field3> {
                                   public:
                                     SEM3DDivergenceEvaluator(const mnl::basics::geometryStack3D& geometry,
                                         const mnl::basics::Matrix& D,
                                         const mnl::basics::Matrix& GLL2G,
                                         const mnl::basics::Vector& weight,
                                         int rank=0, int size=1);
                                     virtual ~SEM3DDivergenceEvaluator();

                                     virtual void evaluate(mnl::basics::matricesStack& res,
                                         const mnl::basics::Field3<mnl::basics::matricesStack>& u) const;
                                     std::vector<mnl::basics::matricesStack*> m_G;
                                     const mnl::basics::geometryStack3D& m_geometry;
                                   protected:
                                     const mnl::basics::Matrix& m_GLL2G;
                                     mnl::basics::Matrix m_Dt;
                                     const mnl::basics::Vector& m_weight;
                                     bool m_mine;
                                     int m_rank;
                                     int m_size;
                                     mnl::basics::componentView<mnl::basics::Matrices,mnl::basics::matricesStack>* m_view;
                                     std::vector<deformed3DDivergenceEvaluator*> m_evals;
                                     mnl::basics::coarseGrid m_division;
                                     mnl::basics::matrixStack* m_W;
                                     mnl::basics::matricesStack* m_Wt;
                                     mnl::basics::matricesStack* m_buf;
                                 };

#endif

