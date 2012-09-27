#ifndef CONVECTION_H_
#define CONVECTION_H_

#include "mnl/geometry.h"
#include "laplacian.h"

class deformedConvectionEvaluator : public 
                                    mnl::utilities::Evaluator<mnl::basics::Matrix> {
                                      public:
                                        deformedConvectionEvaluator(const mnl::basics::matrixStack& G,
                                            const mnl::basics::Matrix& D,
                                            const mnl::basics::Vector& weight,
                                            const mnl::basics::Field2<mnl::basics::Matrix>& u,
                                            mnl::basics::Matrix& buffer,
                                            mnl::basics::Matrix& buffer2);

                                        virtual void evaluate(mnl::basics::Matrix& res, const mnl::basics::Matrix& u) const;
                                      protected:
                                        const mnl::basics::matrixStack& m_G;
                                        const mnl::basics::Vector& m_weight;
                                        const mnl::basics::Matrix& m_D;
                                        const mnl::basics::Field2<mnl::basics::Matrix>& m_u;
                                        mnl::basics::Matrix& m_buffer;
                                        mnl::basics::Matrix& m_buffer2;
                                    };

class deformed3DConvectionEvaluator : public
                                      mnl::utilities::Evaluator<mnl::basics::Matrices> {
                                        public:
                                          deformed3DConvectionEvaluator(const mnl::basics::matricesStack& G,
                                              const mnl::basics::Matrix& D,
                                              const mnl::basics::Vector& weight,
                                              const mnl::basics::Field3<mnl::basics::Matrices>& u,
                                              mnl::basics::Matrices& buffer,
                                              mnl::basics::Matrices& buffer2);

                                          virtual void evaluate(mnl::basics::Matrices& res, const mnl::basics::Matrices& u) const;
                                        protected:
                                          const mnl::basics::matricesStack& m_G;
                                          const mnl::basics::Vector& m_weight;
                                          const mnl::basics::Matrix& m_D;
                                          const mnl::basics::Field3<mnl::basics::Matrices>& m_u;
                                          mnl::basics::Matrices& m_buffer;
                                          mnl::basics::Matrices& m_buffer2;
                                        private:
                                          void evaluateComponent(mnl::basics::Matrices& res,
                                              const mnl::basics::Matrices& field,
                                              int index) const;
                                      };

class convectionEvaluator {
  public:
    convectionEvaluator(const mnl::basics::geometryStack& geometry,
        const mnl::basics::Matrix& D,
        const mnl::basics::Vector& weight,
        int rank=0, int size=1,
        const mnl::basics::geometryStack3D* geometry3D=NULL);
    ~convectionEvaluator();

    void evaluate(mnl::basics::matrixStack& res, const mnl::basics::matrixStack& u) const;
    void evaluate(mnl::basics::matricesStack& res, const mnl::basics::matricesStack& u);
    void evaluateDeformed(mnl::basics::matricesStack& res,
        const mnl::basics::matricesStack& u);
    template <class M, class M2>
      void solve(M& res, const M& phi, const mnl::utilities::ringBuffer<M2>& u, int n,
          Real Dt, Real T, Real t0, Real Dtau=0)
      {
        if( Dtau == 0 )
          Dtau = Dt;
        int steps=(T-t0)/Dtau;
        M* buf[2];
        buf[0] = mnl::utilities::g_manager.clone(phi);
        buf[1] = mnl::utilities::g_manager.clone(phi);
        M* temp   = mnl::utilities::g_manager.clone(phi);
        M2* interp = mnl::utilities::g_manager.clone(u.get(n));
        m_interp = (void*)interp;
        *buf[0] = phi;
        int i;
        for( i=0;i<steps;++i ) {
          extrapolate(*interp,u,n,Dt,t0+i*Dtau);
          /* first step RK4 */
          evaluate(*temp,*buf[i%2]);
          *buf[(i+1)%2] = *temp;
          *temp *= Dtau/Real(2);
          *temp += *buf[i%2];

          extrapolate(*interp,u,n,Dt,t0+(i+.5f)*Dtau);
          /* second step RK4 */
          evaluate(res,*temp);
          buf[(i+1)%2]->axpy(Real(2),res);
          res *= Dtau/Real(2);
          res += *buf[i%2];

          /* third step RK4 */
          evaluate(*temp,res);
          buf[(i+1)%2]->axpy(Real(2),*temp);
          *temp *= Dtau;
          *temp += *buf[i%2];

          extrapolate(*interp,u,n,Dt,t0+(i+1)*Dtau);
          /* fourth step RK4 */
          evaluate(res,*temp);
          *buf[(i+1)%2] += res;
          *buf[(i+1)%2] *= Dtau/Real(6);
          *buf[(i+1)%2] += *buf[i%2];
        }
        if( m_rank == 0 )
          std::cout << i << std::endl;
        res = *buf[i%2];

        mnl::utilities::g_manager.unlock(temp);
        mnl::utilities::g_manager.unlock(buf[0]);
        mnl::utilities::g_manager.unlock(buf[1]);
        mnl::utilities::g_manager.unlock(interp);
        m_interp = NULL;
      }

    template<class M, class M2>
      void solve(mnl::basics::Field3<M>& res,
          const mnl::basics::Field3<M>& phi,
          const mnl::utilities::ringBuffer<M2>& u,
          int n, Real Dt, Real T, Real t0, Real Dtau)
      {
        solve(res.X(),phi.X(),u,n,Dt,T,t0,Dtau);
        solve(res.Y(),phi.Y(),u,n,Dt,T,t0,Dtau);
        solve(res.Z(),phi.Z(),u,n,Dt,T,t0,Dtau);
      }

    SEMLaplacianEvaluator::BC m_bc;
    bool m_deformed;
    const mnl::basics::geometryStack3D* m_geometry3D;
    mnl::basics::componentView<mnl::basics::Matrices,mnl::basics::matricesStack>* m_view;
    int m_order;
  protected:
    template <class T>
      void extrapolate(T& res, const mnl::utilities::ringBuffer<T>& input,
          int n, Real Dt, Real t)
      {
        if( m_order == 0 )
          res = input.get(n);
        else {
          Real scale = t/Dt;
          res = input.get(n,0);
          res -= input.get(n,-1);
          res *= scale;
          res += input.get(n,0);
        }
      }

    const mnl::basics::geometryStack& m_geometry;
    const mnl::basics::Matrix& m_D;
    const mnl::basics::Vector& m_weight;
    mnl::basics::matrixStack* m_buffer;
    mnl::basics::matrixStack* m_buffer2;
    mnl::basics::matricesStack* m_buffer3;
    mnl::basics::matricesStack* m_buffer4;
    void* m_interp;
    mnl::basics::coarseGrid m_grid;
    int m_rank;
    int m_size;
};

#endif

