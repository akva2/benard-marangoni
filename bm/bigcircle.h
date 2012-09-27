#ifndef BIGCIRCLE_H_
#define BIGCIRCLE_H_

#include "mnl/geometry.h"
#include "mnl/hdf5.h"

class bigCircleGeometry : public mnl::basics::geometryStack {
  public:
    using geometryStack::dssum;
    using geometryStack::periodicDssum;

    bigCircleGeometry(int N, int M, const mnl::basics::Vector& ggrid, const mnl::basics::Vector& weight);
    virtual ~bigCircleGeometry();

    void innerDssum(mnl::basics::matrixStack& element, int group) const;
    virtual void dssum(mnl::basics::matrixStack& element, int rank=0, int size=1, int tag=0, void* group=NULL) const;
    virtual void periodicDssum(mnl::basics::matrixStack& element) const;
    inline void  dssum_row_to_row(mnl::basics::matrixStack& element, int group1, 
        int group2, int edge1, int edge2, int row1, int row2, int skip=0) const;
    inline void  dssum_col_to_col(mnl::basics::matrixStack& element, int group1, int group2,
        int edge1, int edge2, int col1, int col2, int skip=0) const;
    inline void  dssum_col_to_row(mnl::basics::matrixStack& element, int group1, 
        int group2, int edge1, int edge2, int col1, int row2) const;
    inline void dssum_row_to_colb(mnl::basics::matrixStack& element, int group1, 
        int group2, int edge1, int edge2, int row1, int col2) const;

    virtual void mask(mnl::basics::matrixStack& op, int rank=0, int size=1) const;
    virtual void mask(mnl::basics::matricesStack& op, int rank=0, int size=1) const;

    virtual std::vector<mnl::basics::coarseGrid> getCoarseGroups(int size) const;
    virtual mnl::basics::matrixStack* getCoarseBuffer(int rows, int cols,
        const mnl::basics::coarseGrid& desc,
        const std::string& name, bool L2=false, bool neumann=false) const;
    virtual mnl::basics::matricesStack* getCoarseBuffer(int rows, int cols, int matrices,
        const mnl::basics::coarseGrid& desc,
        const std::string& name, bool L2=false, bool neumann=false) const;
    virtual void fineToCoarse(mnl::basics::matrixStack& result,
        const mnl::basics::matrixStack& u,
        const mnl::basics::coarseGrid& desc,
        bool neumann) const;
    virtual void fineToCoarse(mnl::basics::matricesStack& result,
        const mnl::basics::matricesStack& u,
        const mnl::basics::coarseGrid& desc,
        bool neumann) const;
    virtual void fineToCoarseL2(mnl::basics::matrixStack& result,
        const mnl::basics::matrixStack& p,
        const mnl::basics::coarseGrid& desc) const;
    virtual void fineToCoarseL2(mnl::basics::matricesStack& result,
        const mnl::basics::matricesStack& p,
        const mnl::basics::coarseGrid& desc) const;
    virtual void fineToCoarseRestriction(mnl::basics::matrixStack& result,
        mnl::basics::matrixStack& buffer,
        const mnl::basics::matrixStack& u,
        const mnl::basics::coarseGrid& desc,
        bool neumann,
        int rank=0, int size=1) const;
    virtual void fineToCoarseRestriction(mnl::basics::matricesStack& result,
        mnl::basics::matricesStack& buffer,
        const mnl::basics::matricesStack& u,
        const mnl::basics::coarseGrid& desc,
        bool neumann,
        int rank=0, int size=1) const;
    virtual void coarseToFine(mnl::basics::matrixStack& u, 
        const mnl::basics::matrixStack& coarse,
        const mnl::basics::coarseGrid& desc,
        bool neumann) const;
    virtual void coarseToFine(mnl::basics::matricesStack& u, 
        const mnl::basics::matricesStack& coarse,
        const mnl::basics::coarseGrid& desc,
        bool neumann) const;
    virtual void coarseToFineL2(mnl::basics::matrixStack& p,
        const mnl::basics::matrixStack& coarse,
        const mnl::basics::coarseGrid& desc) const;
    virtual void coarseToFineL2(mnl::basics::matricesStack& p,
        const mnl::basics::matricesStack& coarse,
        const mnl::basics::coarseGrid& desc) const;

    virtual void dssumCoarse(mnl::basics::matricesStack& result,
        const mnl::basics::coarseGrid& groups,
        int rank=0, int size=1, int tag=0) const
    {
      for( int l=0;l<result[0].matrices();++l )
        dssumCoarse(result.at(l),groups,rank,size,tag);
    }

    virtual void dssumCoarse(mnl::basics::matrixStack& result,
        const mnl::basics::coarseGrid& groups,
        int rank=0, int size=1, int tag=0) const;

    virtual bool hasCoarseSolver() const
    {
      return( true );
    }
    std::string getOperatorFile() const
    {
      return "bigcircle.hdf5";
    }

    virtual mnl::basics::coarseDescriptor getRestrictionGridInfo() const;
    virtual mnl::basics::coarseGrid getDivisionInfo(int size) const;

    void setupRestrictionOperators(std::vector<mnl::basics::Matrix*>& result,
        const mnl::basics::Vector& coarseGrid,
        const std::vector<const mnl::basics::Vector*> grids) const;

    virtual void getRestriction(mnl::basics::Vector& result,
        mnl::basics::matrixStack& work,
        const mnl::basics::matrixStack& LG,
        const mnl::basics::matrixStack& input,
        mnl::basics::Matrix& temp,
        const mnl::basics::coarseGrid& group,
        const std::vector<mnl::basics::Matrix*>& RT);
    virtual void getRestriction(mnl::basics::Vector& result,
        mnl::basics::matricesStack& work,
        mnl::basics::matricesStack& work2,
        const mnl::basics::matricesStack& LG,
        const mnl::basics::matricesStack& input,
        mnl::basics::Matrix& temp,
        const mnl::basics::coarseGrid& group,
        const std::vector<mnl::basics::Matrix*>& RT);

    virtual void gatherAndRestrict(mnl::basics::Vector& result,
        mnl::basics::matricesStack& work,
        mnl::basics::matricesStack& work2,
        mnl::basics::matricesStack& work3,
        const mnl::basics::matricesStack& LG,
        const mnl::basics::matricesStack& input,
        mnl::basics::Matrix& temp,
        const std::vector<mnl::basics::Matrix*>& RT,
        int rank, int size);

    virtual void getProlongiation(mnl::basics::matrixStack& result,
        mnl::basics::matrixStack& work,
        const mnl::basics::matrixStack& LG,
        const mnl::basics::Vector& input,
        mnl::basics::Matrix& temp,
        const std::vector<mnl::basics::Matrix*>& RT);
    virtual void getProlongiation(mnl::basics::matricesStack& result,
        mnl::basics::matricesStack& work,
        mnl::basics::matricesStack& work2,
        const mnl::basics::matricesStack& LG,
        const mnl::basics::Vector& input,
        mnl::basics::Matrix& temp,
        const std::vector<mnl::basics::Matrix*>& RT);

    void doProlongiation(mnl::basics::matrixStack& result,
        mnl::basics::matrixStack& work,
        mnl::basics::Matrix& temp,
        const std::vector<mnl::basics::Matrix*>& RT);

    void prolongAndScatter(mnl::basics::matricesStack& result,
        mnl::basics::matricesStack& work,
        mnl::basics::matricesStack& work2,
        mnl::basics::matricesStack& work3,
        const mnl::basics::matricesStack& LG,
        const mnl::basics::Vector& input,
        mnl::basics::Matrix& temp,
        const std::vector<mnl::basics::Matrix*>& RT,
        int rank, int size);


    void dssum_row_to_rowMPI(mnl::basics::matrixStack& element, 
        Real* temp, char* temp2,
        int group, int edge, int row, int& pos, bool dosum=false) const;
    void dssum_col_to_colMPI(mnl::basics::matrixStack& element,
        Real* temp, char* temp2,
        int group, int edge, int col, int& pos, bool dosum=false) const;

    int m_division;
    int m_size1;
    int m_size2;
    std::vector<mnl::utilities::Range> m_range;

    void processFineToCoarse(mnl::basics::Matrix& foo, 
        const mnl::basics::matrixStack& element, 
        const mnl::basics::coarseDescriptor& desc, bool neumann) const;
    void processCoarseToFine(mnl::basics::matrixStack& element,
        const mnl::basics::Matrix& foo,
        const mnl::basics::coarseDescriptor& desc, bool neumann) const;
    void processFineToCoarseL2(mnl::basics::Matrix& foo,
        const mnl::basics::matrixStack& element, 
        const mnl::basics::coarseDescriptor& desc) const;
    void processCoarseToFineL2(mnl::basics::matrixStack& p,
        const mnl::basics::Matrix& foo, 
        const mnl::basics::coarseDescriptor& desc) const;
    void processFineToCoarseRestrictionInner(mnl::basics::matrixStack& buffer, int n) const;
    void processFineToCoarseRestrictionOuter1(mnl::basics::matrixStack& buffer, int n) const;
    void processFineToCoarseRestrictionOuter2(mnl::basics::matrixStack& buffer, int n) const;
    void dssumRowToRowMPI(mnl::basics::Matrix& result, char* temp, Real* temp2,
        int row, int& pos, bool dosum=false) const;
    void dssumColToColMPI(mnl::basics::Matrix& result, char* temp, Real* temp2,
        int col, int& pos, bool dosum=false) const;
  private:
    void addInnerGroup(mnl::basics::coarseGrid& result, int n, int offs=0) const;
    void addOuterGroup1(mnl::basics::coarseGrid& result, int n, int offs=0) const;
    void addOuterGroup2(mnl::basics::coarseGrid& result, int n, int offs=0) const;
    void addOuterGroup3(mnl::basics::coarseGrid& result, int n, int offs=0) const;
    void addOuterGroup4(mnl::basics::coarseGrid& result, int n, int offs=0) const;
};

class bigCircleGeometry3D : public mnl::basics::geometryStack3D {
  public:
    bigCircleGeometry3D(const mnl::basics::Vector& weight, bigCircleGeometry& geom2D);
    virtual ~bigCircleGeometry3D();

    virtual void dssum(mnl::basics::matricesStack& op, int rank=0, int size=1) const;
    virtual void dssum(mnl::basics::matrixStack& op, int rank=0, int size=1) const;
    virtual void mask(mnl::basics::matricesStack& op, int rank=0, int size=1) const;
    virtual void mask(mnl::basics::matrixStack& op, int rank=0, int size=1) const;
    virtual mnl::basics::coarseGrid getDivisionInfo(int size) const;
    virtual std::vector<mnl::basics::coarseGrid> getCoarseGroups(int size) const;
    virtual mnl::basics::matricesStack* getCoarseBuffer(int rows, int cols, int matrices,
        const mnl::basics::coarseGrid& desc,
        const std::string& name, bool L2=false) const;
    virtual void fineToCoarse(mnl::basics::matricesStack& result,
        const mnl::basics::matricesStack& u,
        const mnl::basics::coarseGrid& desc, bool neumann) const;
    virtual void fineToCoarseL2(mnl::basics::matricesStack& result,
        const mnl::basics::matricesStack& p,
        const mnl::basics::coarseGrid& desc) const;
    virtual void fineToCoarseRestriction(mnl::basics::matricesStack& result,
        mnl::basics::matricesStack& buffer,
        const mnl::basics::matricesStack& u,
        const mnl::basics::coarseGrid& desc,
        int rank=0, int size=1) const;
    virtual void coarseToFine(mnl::basics::matricesStack& u, 
        const mnl::basics::matricesStack& coarse,
        const mnl::basics::coarseGrid& desc, bool neumann) const;
    virtual void coarseToFineL2(mnl::basics::matricesStack& p,
        const mnl::basics::matricesStack& coarse,
        const mnl::basics::coarseGrid& desc) const;
    virtual void dssumCoarse(mnl::basics::matricesStack& result,
        const mnl::basics::coarseGrid& groups,
        int rank=0, int size=1, int tag=0) const;

    virtual bool hasCoarseSolver() const
    {
      return( true );
    }

    std::string getOperatorFile() const
    {
      return "bigcircle.hdf5";
    }

    virtual mnl::basics::coarseDescriptor getRestrictionGridInfo() const;

    virtual void getRestriction(mnl::basics::Vector& result,
        mnl::basics::matricesStack& work,
        mnl::basics::matricesStack& work2,
        const mnl::basics::matricesStack& LG,
        const mnl::basics::matricesStack& input,
        mnl::basics::Matrix& temp,
        const mnl::basics::coarseGrid& group,
        const std::vector<mnl::basics::Matrix*>& RT);
    virtual void getProlongiation(mnl::basics::matricesStack& result,
        mnl::basics::matricesStack& work,
        mnl::basics::matricesStack& work2,
        const mnl::basics::matricesStack& LG,
        const mnl::basics::Vector& input,
        mnl::basics::Matrix& temp,
        const std::vector<mnl::basics::Matrix*>& RT);
  protected:
    bigCircleGeometry& m_geometry2D;
    int m_division;
    int m_size1;
    int m_size2;
};

#endif

