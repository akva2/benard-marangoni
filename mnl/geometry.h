/***************************************************************************
 *   Copyright (C) 2005-2008 by Arne Morten Kvarving                       *
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
#ifndef MNL_GEOMETRY_H_
#define MNL_GEOMETRY_H_

#include "function.h"
#include "gordonhall.h"
#include "buffers.h"
#include "memtracker.h"

#include <string>
#include <map>

#define AXPY_ROW_FROM_ROW_NO_EDGES(y,x,boffs) \
  BLASRPFX(axpy,BLASINT(N2),BLASREAL(mnlRealOne),BLASPREAL(x)+N*boffs,BLASINT(N),BLASPREAL(y)+N*boffs,BLASINT(N));

#define AXPY_ROW_FROM_ROW(y,x) \
  BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealOne),BLASPREAL(x),BLASINT(N),BLASPREAL(y),BLASINT(N));

#define COPY_ROW_FROM_ROW_NO_EDGES(y,x,boffs) \
  BLASRPFX(copy,BLASINT(N2),BLASPREAL(x)+N*boffs,BLASINT(N),BLASPREAL(y)+N*boffs,BLASINT(N));

#define COPY_ROW_FROM_ROW(y,x) \
  BLASRPFX(copy,BLASINT(N),BLASPREAL(x),BLASINT(N),BLASPREAL(y),BLASINT(N));

#define COPY_ROW_FROM_COL(y,x) \
  BLASRPFX(copy,BLASINT(N),BLASPREAL(x),BLASINT(mnlIntOne),BLASPREAL(y),BLASINT(N));

#define AXPY_ROW_FROM_COL(y,x) \
  BLASRPFX(axpy,BLASINT(N2),BLASREAL(mnlRealOne),BLASPREAL(x),BLASINT(mnlIntOne),BLASPREAL(y),BLASINT(N));

#define COPY_ROW_FROM_COL_BACKWARDS(y,x) \
  BLASRPFX(copy,BLASINT(N),BLASPREAL(x),BLASINT(mnlIntMinusOne),BLASPREAL(y),BLASINT(N));

#define AXPY_ROW_FROM_COL_BACKWARDS(y,x) \
  BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealOne),BLASPREAL(x),BLASINT(mnlIntMinusOne),BLASPREAL(y),BLASINT(N));

#define COPY_COL_FROM_COL(y,x) \
  BLASRPFX(copy,BLASINT(N),BLASPREAL(x),BLASINT(mnlIntOne),BLASPREAL(y),BLASINT(mnlIntOne));
#define AXPY_COL_FROM_COL(y,x) \
  BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealOne),BLASPREAL(x),BLASINT(mnlIntOne),BLASPREAL(y),BLASINT(mnlIntOne));

#define COPY_COL_FROM_ROW(y,x) \
  BLASRPFX(copy,BLASINT(N2),BLASPREAL(x),BLASINT(N),BLASPREAL(y),BLASINT(mnlIntOne));

#define COPY_COL_FROM_ROW_BACKWARDS(y,x) \
  BLASRPFX(copy,BLASINT(N),BLASPREAL(x),BLASINT(N),BLASPREAL(y),BLASINT(mnlIntMinusOne));

#define AXPY_NO_EDGES(y,x,loffs) \
  BLASRPFX(axpy,BLASINT(N2),BLASREAL(mnlRealOne),BLASPREAL(x)+loffs,BLASINT(mnlIntOne),BLASPREAL(y)+loffs,BLASINT(mnlIntOne));

#define COPY_NO_EDGES(y,x,loffs) \
  BLASRPFX(copy,BLASINT(N2),BLASPREAL(x)+loffs,BLASINT(mnlIntOne),BLASPREAL(y)+loffs,BLASINT(mnlIntOne));

#define AXPY_BACKWARDS(y,x,loffs) \
  BLASRPFX(axpy,BLASINT(N2),BLASREAL(mnlRealOne),BLASPREAL(x)+loffs,BLASINT(mnlIntMinusOne),BLASPREAL(y),BLASINT(mnlIntOne));

#define COPY_BACKWARDS(y,x,loffs) \
  BLASRPFX(copy,BLASINT(N2),BLASPREAL(x)+loffs,BLASINT(mnlIntMinusOne),BLASPREAL(y)+loffs,BLASINT(mnlIntOne));

#define AXPY_ROW_TO_COLUMN(y,x) \
  BLASRPFX(axpy,BLASINT(N),BLASREAL(mnlRealOne),BLASPREAL(x),BLASINT(N),BLASPREAL(y),BLASINT(mnlIntOne));

namespace mnl {
  namespace basics {
    class spectralElement2D {
      public:
        spectralElement2D(int N, int M) :
          m_gh(N,M)
      {
      }

        spectralElement2D(int N, int M, function1Dto2D& edge1,function1Dto2D& edge2,
            function1Dto2D& edge3,function1Dto2D& edge4) :
          m_gh(N,M,edge1,edge2,edge3,edge4)
      {
      }

        utilities::GordonHall& getGH()
        {
          return m_gh;
        }

        const utilities::GordonHall& getGH() const
        {
          return m_gh;
        }

        void mass(Matrix& res, const Vector& weight) const;
        void invMass(Matrix& res, const Vector& weight) const;

        void evaluate(Matrix& res, const Vector& grid,
            const Vector& ggrid, const function2D& source,
            Real t, Real Lz=2) const;

        void evaluate(Matrix& res, const Vector& grid, const Vector& ggrid,
            const function3D& source, Real z, Real t) const;

        void evaluate(Matrices& res, const Vector& grid, const Vector& ggrid,
            const function3D& source, Real t, Real Lz=2) const;

        Real getVolume(const Vector& weight) const;
        std::pair<Real,Real> getSize(const Vector& weight, const Vector& grid) const;
        std::pair<Real,Real> getAdjustedSize(const Vector& weight, const Vector& grid) const;
      protected:
        utilities::GordonHall m_gh;
    };

    class spectralElement3D {
      public:
        spectralElement3D(int N, int M, int P) :
          m_gh(N,M,P)
      {
      }

        utilities::GordonHall3D& getGH()
        {
          return m_gh;
        }

        const utilities::GordonHall3D& getGH() const
        {
          return m_gh;
        }

        void evaluate(Matrices& res, const function3D& source, Real t) const;

        void mass(Matrices& res, const Vector& weight) const;

        Real getVolume(const Vector& weight) const;
        std::vector<Real> getSize(const Vector& weight, const Vector& grid) const;
        std::vector<Real> getAdjustedSize(const Vector& weight, const Vector& grid) const;
      protected:
        utilities::GordonHall3D m_gh;
    };

    typedef struct {
      enum OverlapType {
        FULL_OVERLAP      = 0,
        NO_BORDER_OVERLAP = 1
      } type;
      enum FlipType {
        NO_FLIP = 0,
        FLIP_XY = 1
      } type2;
      enum CopyType {
        NORMAL					= 0,
        BACKWARDS_ROW			= 1,
        BACKWARDS_COL_BACKWARDS = 2,
        ROW_BACKWARDS 	        = 3

      } type3;
      std::vector<int> elements;
      int size1; /* coarse # points */
      int size2; /* something */
      int size3; /* coarse domains */
      int size4; /* neumann coarse operator size */
    } coarseDescriptor;

    typedef std::vector<coarseDescriptor> coarseGrid;

    class geometryStack {
      public:
        class geometryDotter {
          public:
            geometryDotter(const geometryStack& stack, bool periodic=false) :
              m_stack(stack), m_periodic(periodic)
          {
          }

            Real operator()(const matrixStack& A, const matrixStack& B) const
            {
              if( m_periodic )
                return m_stack.dot(A,B,m_stack.m_multP);
              else
                return m_stack.dot(A,B,m_stack.m_mult);
            }

            Real operator()(const matricesStack& A, const matricesStack& B) const
            {
              return m_stack.dot(A,B);
            }

            Real operator()(const matrixStack& A, const matrixStack& B, const std::vector<int>& r) const
            {
              Real result=0;
              matrixStack* mult = m_periodic?m_stack.m_multP:m_stack.m_mult;
              for( int i=0;i<A.size();++i )
                result += m_stack.dot(A[i],B[i],r[i],&((*mult)[i]));

              return( result );
            }

            Real operator()(const matricesStack& A, const matricesStack& B, const std::vector<int>& r) const
            {
              Real result=0;
              for( int i=0;i<A.size();++i )
                result += m_stack.dot(A[i],B[i],r[i]);

              return( result );
            }
          protected:
            const geometryStack& m_stack;
            bool m_periodic;
        };


        geometryStack(const Vector& grid, const Vector& weight);

        virtual ~geometryStack();

        virtual std::string getOperatorFile() const
        {
          return("nofile");
        }

        void setLz(Real Lz)
        {
          m_Lz = Lz;
        }

        std::vector<spectralElement2D*>& getGridVector()
        {
          return m_grid;
        }

        const geometryDotter& getDotter(bool periodic=false) const
        {
          if( periodic )
            return m_dotterP;

          return m_dotter;
        }

        std::vector<matrixStack*> getInterpolatedGeometryDerivatives(const Matrix& GLL2G) const;

        const matrixStack& getMultiplicities() const
        {
          return( *m_mult );
        }

        virtual void dssum(matrixStack& op, int rank=0, int size=1, int tag=0, void* group=NULL) const = 0;
        virtual void periodicDssum(matrixStack& op) const
        {
        }
        void dssum(matricesStack& element, int rank=0, int size=1, int tag=0, void* group=NULL) const
        {
          for( int l=0;l<element[0].matrices();++l )
            dssum(element.at(l),rank,size,tag,group);
        }

        void periodicDssum(matricesStack& element) const
        {
          for( int l=0;l<element[0].matrices();++l )
            periodicDssum(element.at(l));
        }

        void dssum(Field2<matrixStack>& op, int rank=0, int size=1, int tag=0, void* group=NULL) const
        {
          dssum(op.X(),rank,size,tag,group);
          dssum(op.Y(),rank,size,tag,group);
        }

        void periodicDssum(Field2<matrixStack>& op) const
        {
          periodicDssum(op.X());
          periodicDssum(op.Y());
        }

        void dssum(Field3<matricesStack>& op, int rank=0, int size=1, int tag=0, void* group=NULL) const
        {
          dssum(op.X(),rank,size,tag,group);
          dssum(op.Y(),rank,size,tag,group);
          dssum(op.Z(),rank,size,tag,group);
        }
        virtual void mask(matrixStack& op, int rank=0, int size=1) const = 0;
        virtual void mask(matricesStack& op, int rank=0, int size=1) const = 0;

        virtual std::vector<coarseGrid> getCoarseGroups(int size) const
        {
          std::vector<coarseGrid> result;

          return( result );
        }

        virtual matrixStack* getCoarseBuffer(int rows, int cols,
            const coarseGrid& desc,
            const std::string& name,
            bool L2=false, bool neumann=false) const
        {
          return( NULL );
        }

        virtual matricesStack* getCoarseBuffer(int rows, int cols, int matrices,
            const coarseGrid& desc,
            const std::string& name, bool L2=false,
            bool neumann=false) const
        {
          return( NULL );
        }

        virtual void fineToCoarse(matrixStack& result,
            const matrixStack& u,
            const coarseGrid& desc,
            bool neumann) const
        {
        }

        virtual void fineToCoarse(matricesStack& result,
            const matricesStack& u,
            const coarseGrid& desc,
            bool neumann) const
        {
        }

        virtual void fineToCoarseL2(matrixStack& result,
            const matrixStack& p,
            const coarseGrid& desc) const
        {
        }
        virtual void fineToCoarseL2(matricesStack& result,
            const matricesStack& p,
            const coarseGrid& desc) const
        {
        }
        virtual void fineToCoarseRestriction(matrixStack& result,
            matrixStack& buffer,
            const matrixStack& u,
            const coarseGrid& desc,
            bool neumann,
            int rank=0, int size=1) const
        {
        }

        virtual void fineToCoarseRestriction(matricesStack& result,
            matricesStack& buffer,
            const matricesStack& u,
            const coarseGrid& desc,
            bool neumann,
            int rank=0, int size=1) const
        {
        }


        virtual void coarseToFine(matrixStack& u,
            const matrixStack& foos,
            const coarseGrid& desc,
            bool neumann) const
        {
        }
        virtual void coarseToFine(matricesStack& u,
            const matricesStack& foos,
            const coarseGrid& desc,
            bool neumann) const
        {
        }
        virtual void coarseToFineL2(matrixStack& p,
            const mnl::basics::matrixStack& foos,
            const coarseGrid& desc) const
        {
        }
        virtual void coarseToFineL2(matricesStack& p,
            const mnl::basics::matricesStack& foos,
            const coarseGrid& desc) const
        {
        }

        virtual void dssumCoarse(matricesStack& result, const coarseGrid& groups,
            int rank=0, int size=1, int tag=0) const
        {
        }

        virtual void dssumCoarse(matrixStack& result,
            const coarseGrid& groups,
            int rank=0, int size=1, int tag=0) const
        {
        }

        virtual bool hasCoarseSolver() const
        {
          return( false );
        }

        virtual coarseDescriptor getRestrictionGridInfo() const
        {
          coarseDescriptor result;
          result.size1 = 0;
          result.size2 = 0;
          result.size3 = 0;

          return( result );
        }

        virtual coarseGrid getDivisionInfo(int siz) const
        {
          coarseGrid result;
          coarseDescriptor desc;
          desc.size1 = size();
          for( int i=0;i<size();++i )
            desc.elements.push_back(i);
          result.push_back(desc);

          return( result );
        }

        virtual void setupRestrictionOperators(std::vector<Matrix*>& result,
            const Vector& coarseGrid,
            const std::vector<const Vector*> grids) const
        {
        }

        virtual void getRestriction(basics::Vector& result,
            basics::matrixStack& work,
            const basics::matrixStack& LG,
            const basics::matrixStack& input,
            basics::Matrix& temp,
            const coarseGrid& group,
            const std::vector<Matrix*>& RT)
        {
        }

        virtual void getRestriction(basics::Vector& result,
            basics::matricesStack& work,
            basics::matricesStack& work2,
            const basics::matricesStack& LG,
            const basics::matricesStack& input,
            basics::Matrix& temp,
            const coarseGrid& group,
            const std::vector<Matrix*>& RT)
        {
        }

        virtual void gatherAndRestrict(basics::Vector& result,
            basics::matricesStack& work,
            basics::matricesStack& work2,
            basics::matricesStack& work3,
            const basics::matricesStack& LG,
            const basics::matricesStack& input,
            basics::Matrix& temp,
            const std::vector<Matrix*>& RT,
            int rank, int size)
        {
        }

        virtual void getProlongiation(basics::matrixStack& result,
            basics::matrixStack& work,
            const basics::matrixStack& LG,
            const basics::Vector& input,
            basics::Matrix& temp,
            const std::vector<Matrix*>& RT)
        {
        }

        virtual void getProlongiation(basics::matricesStack& result,
            basics::matricesStack& work,
            basics::matricesStack& work2,
            const basics::matricesStack& LG,
            const basics::Vector& input,
            basics::Matrix& temp,
            const std::vector<Matrix*>& RT)
        {
        }

        virtual void prolongAndScatter(basics::matricesStack& result,
            basics::matricesStack& work,
            basics::matricesStack& work2,
            basics::matricesStack& work3,
            const basics::matricesStack& LG,
            const basics::Vector& input,
            basics::Matrix& temp,
            const std::vector<Matrix*>& RT,
            int rank, int size)
        {
        }

        inline void maskField(Field2<matrixStack>& op, int rank=0, int size=1) const
        {
          mask(op.X(),rank,size);
          mask(op.Y(),rank,size);
        }

        inline void maskField(Field3<matricesStack>& op, int rank=0, int size=1) const
        {
          mask(op.X(),rank,size);
          mask(op.Y(),rank,size);
          mask(op.Z(),rank,size);
        }

        virtual void mass(matricesStack& u, bool dosum=true) const;
        virtual void mass(matrixStack& res, bool dosum=true) const;
        virtual void mass(matrixStack& res, const std::vector<int>& r) const;
        virtual void mass(matricesStack& res, const std::vector<int>& r, int start=0) const;
        inline void mass(Field2<matrixStack>& res, bool dosum=true) const
        {
          mass(res.X(),dosum);
          mass(res.Y(),dosum);
        }
        inline void mass(Field2<matrixStack>& res, const std::vector<int>& r, int size=1, int plane=1) const
        {
          mass(res.X(),r);
          mass(res.Y(),r);
        }
        inline void mass(Field3<matricesStack>& res, bool dosum=true) const
        {
          mass(res.X(),dosum);
          mass(res.Y(),dosum);
          mass(res.Z(),dosum);
        }
        inline void mass(Field3<matricesStack>& res, const std::vector<int>& r) const
        {
          mass(res.X(),r);
          mass(res.Y(),r);
          mass(res.Z(),r);
        }
        virtual void invMass(matrixStack& res) const;
        virtual void invMassP(matrixStack& res) const;
        virtual void invMass(matrixStack& res, const std::vector<int>& r) const;
        virtual void invMassP(matrixStack& res, const std::vector<int>& r) const;
        virtual void invMass(matricesStack& res) const;
        virtual void invMassP(matricesStack& res) const;
        virtual void invMass(matricesStack& res, const std::vector<int>& r) const;

        inline void invMass(Field2<matrixStack>& res) const
        {
          invMass(res.X());
          invMass(res.Y());
        }
        inline void invMass(Field2<matrixStack>& res, const std::vector<int>& r) const
        {
          invMass(res.X(),r);
          invMass(res.Y(),r);
        }
        inline void invMass(Field3<matricesStack>& res) const
        {
          invMass(res.X());
          invMass(res.Y());
          invMass(res.Z());
        }
        inline void invMass(Field3<matricesStack>& res, const std::vector<int>& r) const
        {
          invMass(res.X(),r);
          invMass(res.Y(),r);
          invMass(res.Z(),r);
        }
        inline void invMassP(Field3<matricesStack>& res) const
        {
          invMassP(res.X());
          invMassP(res.Y());
          invMassP(res.Z());
        }

        template<class T, class TFUNC>
          void evaluate(T& U, const Vector& grid, const TFUNC& source, Real t) const
          {
            int max=m_grid.size();
#pragma omp parallel for schedule(static)
            for( int i=0;i<max;++i )
              m_grid[i]->evaluate(U[i],grid,m_ggrid,source,t,m_Lz);
          }

        inline void evaluate(basics::Field2<basics::matrixStack>& U,
            const Vector& grid,
            const basics::function2D& sourcex,
            const basics::function2D& sourcey, Real t) const
        {
          evaluate(U.X(),grid,sourcex,t);	
          evaluate(U.Y(),grid,sourcey,t);	
        }

        inline void evaluate(basics::Field3<basics::matricesStack>& U,
            const Vector& grid,
            const basics::function3D& sourcex,
            const basics::function3D& sourcey,
            const basics::function3D& sourcez, Real t) const
        {
          evaluate(U.X(),grid,sourcex,t);	
          evaluate(U.Y(),grid,sourcey,t);	
          evaluate(U.Z(),grid,sourcez,t);	
        }
        template<class T, class TFUNC>
          void evaluate(T& U, const Vector& grid,
              const TFUNC& source, Real t, const std::vector<int>& r) const
          {
            int max=r.size();
#pragma omp parallel for schedule(static)
            for( int i=0;i<max;++i )
              m_grid[r[i]]->evaluate(U[i],grid,m_ggrid,source,t,m_Lz);
          }
        template<class T, class TFUNC>
          void evaluate(T& U, const Vector& grid, const TFUNC& source,
              Real z, Real t, const std::vector<int>& r) const
          {
            int max=r.size();
#pragma omp parallel for schedule(static)
            for( int i=0;i<max;++i )
              m_grid[r[i]]->evaluate(U[i],grid,m_ggrid,source,z,t,m_Lz);
          }
        inline void evaluate(basics::Field2<basics::matrixStack>& U,
            const Vector& grid,
            const basics::function2D& sourcex,
            const basics::function2D& sourcey, Real t,
            const std::vector<int>& r) const
        {
          evaluate(U.X(),grid,sourcex,t,r);	
          evaluate(U.Y(),grid,sourcey,t,r);	
        }
        inline void evaluate(basics::Field3<basics::matricesStack>& U,
            const Vector& grid,
            const basics::function3D& sourcex,
            const basics::function3D& sourcey,
            const basics::function3D& sourcez,
            Real t, const std::vector<int>& r) const
        {
          evaluate(U.X(),grid,sourcex,t,r);	
          evaluate(U.Y(),grid,sourcey,t,r);	
          evaluate(U.Z(),grid,sourcez,t,r);	
        }

        void localToGlobal(Vector& result, const basics::matrixStack& input,
            const basics::matrixStack& LG) const;
        void globalToLocal(basics::matrixStack& result, const Vector& input,
            const basics::matrixStack& LG) const;
        void localToGlobal(Vector& result, const basics::matricesStack& input,
            const basics::matricesStack& LG) const;
        void globalToLocal(basics::matricesStack& result, const Vector& input,
            const basics::matricesStack& LG) const;

        int size() const
        {
          return m_grid.size();
        }

        int rows() const
        {
          return m_grid[0]->getGH().getMapping().rows();
        }

        int cols() const
        {
          return m_grid[0]->getGH().getMapping().cols();
        }

        inline spectralElement2D& operator [](int index)
        {
          assert( index < m_grid.size() );

          return( *m_grid[index] );
        }

        inline const spectralElement2D& operator [](int index) const
        {
          assert( index < m_grid.size() );

          return( *m_grid[index] );
        }

        void compute(Real t);

        void computeMultiplicities();
        void computeMass();

        Real dot(const matrixStack& in1, const matrixStack& in2, const matrixStack* mult=NULL) const;
        Real dot(const matricesStack& in1, const matricesStack& in2) const;
        Real dot(const Matrix& in1, const Matrix& in2, int i, const Matrix* mult=NULL) const;
        Real dot(const Matrices& in1, const Matrices& in2, int i) const;
        matrixStack* m_mass;
        matrixStack* m_imass;
        matrixStack* m_imassP;
        matrixStack* m_mult;
        matrixStack* m_multP;
        Real m_Lz;

        static void interpolate(mnl::basics::Matrix& result,
            const mnl::basics::Matrix& G,
            const mnl::basics::Matrix& GLL2G);
      protected:
        std::vector<spectralElement2D*> m_grid;
        geometryDotter m_dotter;
        geometryDotter m_dotterP;
        const Vector& m_ggrid;
        const Vector& m_weight;
        bool m_mine;
    };

    class geometryStack3D {
      public:
        class geometryDotter {
          public:
            geometryDotter(const geometryStack3D& stack) :
              m_stack(stack)
          {
          }

            Real operator()(const matricesStack& A, const matricesStack& B) const
            {
              return m_stack.dot(A,B);
            }

            Real operator()(const matricesStack& A, const matricesStack& B, const std::vector<int>& r) const
            {
              Real result=0;
              for( int i=0;i<A.size();++i )
                result += m_stack.dot(A[i],B[i],r[i]);

              return( result );
            }
          protected:
            const geometryStack3D& m_stack;
        };

        geometryStack3D(const mnl::basics::Vector& m_weight);
        virtual ~geometryStack3D();

        std::vector<spectralElement3D*>& getGridVector()
        {
          return m_grid;
        }

        void computeMultiplicities();
        void computeMass();
        void computeReferenceGeometryDerivatives();

        int size() const
        {
          return m_grid.size();
        }

        Real max() const
        {
          Real result=0;
          for( int i=0;i<size();++i )
            if( m_grid[i]->getGH().getMapping().X().max() > result )
              result =  m_grid[i]->getGH().getMapping().X().max();

          return( result );
        }

        inline spectralElement3D& operator [](int index)
        {
          assert( index < m_grid.size() );

          return( *m_grid[index] );
        }

        inline const spectralElement3D& operator [](int index) const
        {
          assert( index < m_grid.size() );

          return( *m_grid[index] );
        }

        Real dot(const matricesStack& in1, const matricesStack& in2) const;
        Real dot(const Matrices& in1, const Matrices& in2, int i) const;

        virtual void mass(matricesStack& res) const;
        inline void mass(Field3<matricesStack>& res) const
        {
          mass(res.X());
          mass(res.Y());
          mass(res.Z());
        }
        virtual void mass(matricesStack& res, const std::vector<int>& r) const;
        inline void mass(Field3<matricesStack>& res, const std::vector<int>& r) const
        {
          mass(res.X(),r);
          mass(res.Y(),r);
          mass(res.Z(),r);
        }

        virtual void invMass(matricesStack& res) const;
        virtual void invMass(matricesStack& res, const std::vector<int>& r) const;
        inline void invMass(Field3<matricesStack>& res) const
        {
          invMass(res.X());
          invMass(res.Y());
          invMass(res.Z());
        }

        inline void invMass(Field3<matricesStack>& res, const std::vector<int>& r) const
        {
          invMass(res.X(),r);
          invMass(res.Y(),r);
          invMass(res.Z(),r);
        }

        virtual void dssum(matricesStack& op, int rank=0, int size=1) const = 0;
        virtual void dssum(matrixStack& op, int rank=0, int size=1) const = 0;
        virtual void mask(matricesStack& op, int rank=0, int size=1) const = 0;
        virtual void mask(matrixStack& op, int rank=0, int size=1) const
        {
        }

        virtual coarseGrid getDivisionInfo(int siz) const
        {
          coarseGrid result;
          coarseDescriptor desc;
          desc.size1 = size();
          for( int i=0;i<size();++i )
            desc.elements.push_back(i);
          result.push_back(desc);

          return( result );
        }

        virtual std::vector<coarseGrid> getCoarseGroups(int size) const
        {
          std::vector<coarseGrid> result;

          return( result );
        }

        virtual matricesStack* getCoarseBuffer(int rows, int cols, int matrices,
            const coarseGrid& desc,
            const std::string& name, bool L2=false) const
        {
          return( NULL );
        }

        virtual void fineToCoarse(matricesStack& result,
            const matricesStack& u,
            const coarseGrid& desc) const
        {
        }

        virtual void fineToCoarseL2(matricesStack& result,
            const matricesStack& p,
            const coarseGrid& desc) const
        {
        }

        virtual void fineToCoarseRestriction(matricesStack& result,
            matricesStack& buffer,
            const matricesStack& u,
            const coarseGrid& desc,
            int rank=0, int size=1) const
        {
        }

        virtual void coarseToFine(matricesStack& u,
            const matricesStack& foos,
            const coarseGrid& desc) const
        {
        }

        virtual void coarseToFineL2(matricesStack& p,
            const mnl::basics::matricesStack& foos,
            const coarseGrid& desc) const
        {
        }

        virtual void dssumCoarse(matricesStack& result, const coarseGrid& groups,
            int rank=0, int size=1, int tag=0) const
        {
        }

        virtual bool hasCoarseSolver() const
        {
          return( false );
        }

        virtual std::string getOperatorFile() const
        {
          return("nofile");
        }

        virtual coarseDescriptor getRestrictionGridInfo() const
        {
          coarseDescriptor result;
          result.size1 = 0;
          result.size2 = 0;
          result.size3 = 0;

          return( result );
        }

        virtual void setupRestrictionOperators(const Vector& coarseGrid,
            const std::vector<const Vector*> grids)
        {
        }

        virtual void getRestriction(basics::Vector& result,
            basics::matricesStack& work,
            basics::matricesStack& work2,
            const basics::matricesStack& LG,
            const basics::matricesStack& input,
            basics::Matrix& temp,
            const coarseGrid& group)
        {
        }

        virtual void getProlongiation(basics::matricesStack& result,
            basics::matricesStack& work,
            basics::matricesStack& work2,
            const basics::matricesStack& LG,
            const basics::Vector& input,
            basics::Matrix& temp)
        {
        }

        const geometryDotter& getDotter() const
        {
          return m_dotter;
        }

        void evaluate(matricesStack& U, const function3D& source, Real t) const
        {
          int max=m_grid.size();
#pragma omp parallel for schedule(static)
          for( int i=0;i<max;++i )
            m_grid[i]->evaluate(U[i],source,t);
        }

        void evaluate(matricesStack& U, const function3D& source, Real t,
            const std::vector<int>& r) const
        {
          int max=r.size();
#pragma omp parallel for schedule(static)
          for( int i=0;i<max;++i )
            m_grid[r[i]]->evaluate(U[i],source,t);
        }

        inline void evaluate(basics::Field3<basics::matricesStack>& U,
            const basics::function3D& sourcex,
            const basics::function3D& sourcey,
            const basics::function3D& sourcez, Real t) const
        {
          evaluate(U.X(),sourcex,t);
          evaluate(U.Y(),sourcey,t);
          evaluate(U.Z(),sourcez,t);
        }

        inline void evaluate(basics::Field3<basics::matricesStack>& U,
            const basics::function3D& sourcex,
            const basics::function3D& sourcey,
            const basics::function3D& sourcez,
            Real t, const std::vector<int>& r) const
        {
          evaluate(U.X(),sourcex,t,r);
          evaluate(U.Y(),sourcey,t,r);
          evaluate(U.Z(),sourcez,t,r);
        }

        std::vector<matricesStack*> getInterpolatedGeometryDerivatives(const Matrix& GLL2G) const;
        const std::vector<matricesStack*>& getReferenceGeometryDerivatives() const;
        std::vector<matricesStack*>& getReferenceGeometryDerivatives();
        std::vector<matricesStack*> getInterpolatedReferenceGeometryDerivatives(const Matrix& GLL2G) const;

        std::vector<mnl::basics::matricesStack*> setupFakeStacks(mnl::basics::matricesStack& op,
            int levels=-1, int perlevel=-1) const;
        void killFakeStacks(std::vector<mnl::basics::matricesStack*> op) const;
        mnl::basics::Field3<mnl::basics::matricesStack>*
          localToGlobalZ(const mnl::basics::Field3<mnl::basics::matricesStack>& input,
              int rank=0, int size=1, bool L2=false, bool dosum=false) const;
        mnl::basics::matricesStack* localToGlobalZ(const mnl::basics::matricesStack& input,
            int rank=0, int size=1, bool L2=false, bool dosum=false,
            mnl::basics::matricesStack* result=NULL) const;
        void globalToLocalZ(mnl::basics::Field3<mnl::basics::matricesStack>& result,
            mnl::basics::Field3<mnl::basics::matricesStack>* input,
            int rank=0, int size=1,bool L2=false) const;
        void globalToLocalZ(mnl::basics::matricesStack& result,
            mnl::basics::matricesStack* input,
            int rank=0, int size=1,bool L2=false) const;

        int m_level;
        int m_perlevel;
        matricesStack* m_mass;
        matricesStack* m_imass;
        matricesStack* m_mult;

        static void interpolate(mnl::basics::Matrices& result,
            const mnl::basics::Matrices& G,
            const mnl::basics::Matrix& GLL2G);
      protected:
        std::vector<spectralElement3D*> m_grid;
        std::vector<basics::matricesStack*> m_referenceDerivatives;
        geometryDotter m_dotter;
        const mnl::basics::Vector& m_weight;
      private:
    };
  }
}

#endif

