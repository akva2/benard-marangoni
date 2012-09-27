#pragma once

#include "mnl/geometry.h"

class quadraticGeometry : public mnl::basics::geometryStack {
  public:
    using geometryStack::dssum;
    using geometryStack::periodicDssum;

    quadraticGeometry(int N, int M, const mnl::basics::Vector& ggrid, const mnl::basics::Vector& weight);
    virtual ~quadraticGeometry();

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

    virtual bool hasCoarseSolver() const
    {
      return( false );
    }
    std::string getOperatorFile() const
    {
      return "quadratic.hdf5";
    }

    int m_size1;
    int m_size2;
};

