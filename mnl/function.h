/***************************************************************************
 *   Copyright (C) 2005 by Arne Morten Kvarving                            *
 *   spiff@mspiggy                                                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef MNL_FUNCTION_H_
#define MNL_FUNCTION_H_

#include "config.h"
#include "field.h"

#include <pthread.h>
#include <sys/time.h>
#include <errno.h>

namespace mnl {
  namespace basics {
    class function1D {
      public:
        virtual Real val(Real x, Real t) const = 0;
        virtual Real diffx(Real x, Real t) 
        {
          return 0;
        }
        virtual Real difft(Real x, Real t) const
        {
          return 0;
        }
        virtual Real diff2x(Real x, Real t) const
        {
          return 0;
        }

        inline Real operator ()(Real x, Real t) const
        {
          return val(x,t);
        }
    };
    class function2D {
      public:
        virtual Real val(Real x, Real y, Real t) const = 0;
        virtual Real diffx(Real x, Real y, Real t) const
        {
          return 0;
        }
        virtual Real diffy(Real x, Real y, Real t) const
        {
          return 0;
        }
        virtual Real difft(Real x, Real y, Real t) const
        {
          return 0;
        }
        virtual Real diff2x(Real x, Real y, Real t) const
        {
          return 0;
        }
        virtual Real diff2y(Real x, Real y, Real t) const
        {
          return 0;
        }

        inline Real operator ()(Real x, Real y, Real t) const
        {
          return val(x,y,t);
        }

        int id() const
        {
          return m_id;
        }

        int m_id;
    };
    class function1Dto2D {
      public:
        virtual std::pair<Real,Real> val(Real x, Real t) const = 0;
        virtual std::pair<Real,Real> diffx(Real x, Real t) const {return std::make_pair<Real,Real>(0.f,0.f);};
        virtual std::pair<Real,Real> difft(Real x, Real t) const {return std::make_pair<Real,Real>(0.f,0.f);};
        virtual std::pair<Real,Real> diff2x(Real x, Real t) const {return std::make_pair<Real,Real>(0.f,0.f);};

        inline std::pair<Real,Real> operator ()(Real x, Real t) const
        {
          return val(x,t);
        }
    };
    class function3D {
      public:
        virtual Real val(Real x, Real y, Real z, Real t) const = 0;
        virtual Real diffx(Real x, Real y, Real z, Real t) const
        {
          return 0;
        }
        virtual Real diffy(Real x, Real y, Real z, Real t) const
        {
          return 0;
        }
        virtual Real diffz(Real x, Real y, Real z, Real t) const
        {
          return 0;
        }
        virtual Real difft(Real x, Real y, Real z, Real t) const
        {
          return 0;
        }
        virtual Real diff2x(Real x, Real y, Real z, Real t) const
        {
          return 0;
        }
        virtual Real diff2y(Real x, Real y, Real z, Real t) const
        {
          return 0;
        }
        virtual Real diff2z(Real x, Real y, Real z, Real t) const
        {
          return 0;
        }

        inline Real operator()(Real x, Real y, Real z, Real t) const
        {
          return val(x,y,z,t);
        }

        int id() const
        {
          return m_id;
        }

        int m_id;
    };

    template<class T>
      class Functor1 {
        public:
          virtual void init(const T& data) {}
          virtual void operator ()(T& input)=0;
          inline void operator()(Field2<T>& input)
          {
            (*this)(input.X());
            (*this)(input.Y());
          }
          inline void operator()(Field3<T>& input)
          {
            (*this)(input.X());
            (*this)(input.Y());
            (*this)(input.Z());
          }
      };
    template<class T>
      class DoNothing1 : public Functor1<T> {
        public:
          virtual void operator()(T& input) {}
      };
    class Vorticity2D : public basics::function2D {
      public:
        Vorticity2D(const basics::function2D& ux, const basics::function2D& uy) :
          m_ux(ux), m_uy(uy)
      {
        m_id = 1;
      }

        Real val(Real x, Real y, Real t) const
        {
          return m_uy.diffx(x,y,t)-m_ux.diffy(x,y,t);
        }
      protected:
        const basics::function2D& m_ux;
        const basics::function2D& m_uy;
    };
  }

  namespace utilities {
    template<class T>
      class Evaluator {
        public:
          Evaluator() {};
          virtual ~Evaluator() {};

          virtual void evaluate(T&, const T&) const=0;
          inline void evaluateField(basics::Field2<T>& res, const basics::Field2<T>& u) const
          {
            evaluate(res.X(),u.X());
            evaluate(res.Y(),u.Y());
          }
          inline void evaluateField(basics::Field3<T>& res, const basics::Field3<T>& u) const
          {
            evaluate(res.X(),u.X());
            evaluate(res.Y(),u.Y());
            evaluate(res.Z(),u.Z());
          }

          virtual bool isL2() const
          {
            return( false );
          }

          virtual void filter(T& foo) const
          {
          }

          virtual void post(T& foo) const
          {
          }

          virtual T* buffer()
          {
            return( NULL );
          }
      };

    class Thread {
      public:
        virtual void run() = 0;
    };

    extern "C" void* runthread(void* ctx);

    template<class T>
      class EvaluatorThread : public Thread {
        public:
          virtual void run()
          {
            m_die=false;
            struct timeval now;
            struct timespec timeout;

            while( !m_die ) {
              pthread_mutex_lock(&m_mutex);
              int j=0;
              do {
                pthread_cond_wait(&m_wake,&m_mutex);
              } while( !m_started );
              pthread_mutex_unlock(&m_mutex);
              if( m_out && m_in && m_eval )
                m_eval->evaluate(*m_out,*m_in);

              pthread_mutex_lock(&m_mutexdone);

              pthread_cond_broadcast(&m_done);
              pthread_mutex_unlock(&m_mutexdone);
            }
          }

          pthread_t m_thread_id;
          pthread_attr_t m_thread_attr;
          pthread_mutex_t m_mutex;
          pthread_mutex_t m_mutexdone;
          pthread_mutexattr_t m_mutex_attr;
          pthread_cond_t m_done;
          pthread_cond_t* m_end;
          pthread_cond_t m_wake;
          int m_id;
          T* m_out;
          const T* m_in;
          Evaluator<T>* m_eval;
          bool m_started;
          bool m_die;
      };

    template<class T>
      class ThreadedEvaluator : public Evaluator<T> {
        public:
          ThreadedEvaluator(int size, std::vector< Evaluator<T>* >& evals) :
            m_evals(evals)
        {
          for( int i=0;i<size;++i ) {
            EvaluatorThread<T>* thread = new EvaluatorThread<T>;
            pthread_cond_init(&thread->m_done,NULL);
            pthread_cond_init(&thread->m_wake,NULL);
            pthread_mutexattr_init(&thread->m_mutex_attr);
            pthread_mutex_init(&thread->m_mutex,&thread->m_mutex_attr);
            pthread_mutex_init(&thread->m_mutexdone,&thread->m_mutex_attr);
            pthread_attr_init(&thread->m_thread_attr);
            pthread_attr_setscope(&thread->m_thread_attr,PTHREAD_SCOPE_SYSTEM);
            thread->m_eval = m_evals[i];
            thread->m_out =  NULL;
            thread->m_in  =  NULL;
            thread->m_id  = i;
            thread->m_started = false;
            pthread_create(&thread->m_thread_id,&thread->m_thread_attr,
                runthread,thread);
            m_thread.push_back(thread);
          }
        }

          virtual ~ThreadedEvaluator()
          {
            for( unsigned int i=0;i<m_thread.size();++i ) {
              m_thread[i]->m_die = true;
              m_thread[i]->m_in = NULL;
              m_thread[i]->m_out = NULL;
              m_thread[i]->m_eval = NULL;
            }
            WakeThreads();
            WaitForThreads();
            for( unsigned int i=0;i<m_thread.size();++i ) {
              pthread_mutex_destroy(&m_thread[i]->m_mutexdone);
              pthread_mutex_destroy(&m_thread[i]->m_mutex);
              pthread_mutexattr_destroy(&m_thread[i]->m_mutex_attr);
              pthread_cond_destroy(&m_thread[i]->m_wake);
              pthread_cond_destroy(&m_thread[i]->m_done);
              pthread_attr_destroy(&m_thread[i]->m_thread_attr);
            }
          }

          virtual void evaluate(T& out, const T& in) const
          {
            (const_cast<ThreadedEvaluator*>(this))->evaluate(out,in);
          }

          void evaluate(T& out, const T& in)
          {
            for( int i=0;i<m_thread.size();++i ) {
              m_thread[i]->m_in = &in;
              m_thread[i]->m_out = &out;
              m_thread[i]->m_eval = m_evals[i];
            }
            WakeThreads();
            WaitForThreads();
            m_evals[0]->post(out);
          }

          std::vector< Evaluator<T>* > m_evals;
        protected:
          std::vector<EvaluatorThread<T>*> m_thread;

          void WakeThreads()
          {
            for( int i=0;i<m_thread.size();++i ) {
              pthread_mutex_lock(&m_thread[i]->m_mutexdone);
              m_thread[i]->m_started = true;
              pthread_cond_broadcast(&m_thread[i]->m_wake);
            }
          }

          void WaitForThreads()
          {
            for( unsigned int i=0;i<m_thread.size();++i ) {
              pthread_cond_wait(&m_thread[i]->m_done,&m_thread[i]->m_mutexdone);
              pthread_mutex_unlock(&m_thread[i]->m_mutexdone);
            }
          }
      };

    template<class T, template<class T3> class T2 >
      class fieldToScalarEvaluator {
        public:
          fieldToScalarEvaluator() {};

          virtual void evaluate(T&, const T2<T>&) const=0;
      };

    template<class T, template<class T3> class T2 >
      class scalarToFieldEvaluator {
        public:
          scalarToFieldEvaluator() {};

          virtual void evaluate(T2<T>&, const T&) const=0;
      };
  }
}

#endif

