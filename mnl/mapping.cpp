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

#include "mapping.h"

using namespace std;
namespace mnl {
  namespace utilities {
    Mapping::Mapping()
    {
    }

    Mapping::~Mapping()
    {
    }

    void Mapping::addElement(int index, const Range& nodes, const Element& element)
    {
      assert( m_mapping.find(index) == m_mapping.end() );

      m_mapping.insert(make_pair(index,make_pair<Range,Element>(nodes,element)));
    }

    void Mapping::print()
    {
      for( map<int,std::pair<Range,Element> >::iterator iter=m_mapping.begin(); iter != m_mapping.end(); ++iter ) {
        std::cout << "Element " << iter->first << ":" << std::endl;
        iter->second.first.print();
        iter->second.second.getOperator().print();
      }
    }

    void Mapping::assemble(basics::Matrix& op)
    {
      for( map<int,pair<Range,Element> >::iterator iter=m_mapping.begin(); iter != m_mapping.end(); ++iter) {
        Element& element = iter->second.second;
        Range& nodes = iter->second.first;
        for( int i=0;i<element.getOperator().cols();++i )
          for( int j=0;j<element.getOperator().rows();++j )
            op[nodes[i]][nodes[j]] += 2.f/element.size()*element.getOperator()[i][j];
      }
    }

    void Mapping::assembleRHS(basics::Vector& vec)
    {
      for( map<int,std::pair<Range,Element> >::iterator iter=m_mapping.begin(); iter != m_mapping.end(); ++iter) {
        Element& element = iter->second.second;
        Range& nodes = iter->second.first;
        basics::Vector rhs = element.size()/2.f*basics::Matrix::diag(element.getWeight())*element.getRHS();
        for( int i=0;i<element.getRHS().length();++i )
          vec[nodes[i]] += rhs[i];
      }
    }

    void Mapping::assembleGrid(basics::Vector& vec)
    {
      int pos = 0;
      for( map<int,pair<Range,Element> >::iterator iter=m_mapping.begin(); iter != m_mapping.end(); ++iter ) {
        Element& element = iter->second.second;	
        for( int i=0;i<element.getGrid().length();++i )
          vec[i+pos] = element.x0()+(element.getGrid()[i]+1)*element.size()/2.f;
        pos += element.getGrid().length()-1;
      }
    }

    Mapping2D::Mapping2D()
    {
    }

    Mapping2D::~Mapping2D()
    {
    }

    void Mapping2D::addElement(const int index, const Range2D& nodes, const Element2D& element)
    {
      assert( m_mapping.find(index) == m_mapping.end() );

      m_mapping.insert(make_pair<int,pair<Range2D,Element2D> >(index,make_pair<Range2D,Element2D>(nodes,element)));
    }

    void Mapping2D::assemble(basics::Matrix& op)
    {
      for( map<int,pair<Range2D,Element2D> >::iterator iter=m_mapping.begin();iter!=m_mapping.end();++iter ) {
        Element2D& element = iter->second.second;
        Range2D& nodes = iter->second.first;
        for( int i=0;i<element.getOperator().cols();++i )
          for( int j=0;j<element.getOperator().rows();++j ) {
            cout << "nodes " << nodes[0][i] << " " << nodes[0][j] << endl;
            op[nodes[0][i]][nodes[0][j]] += 2.f/element.width()*2.f/element.height()*element.getOperator()[i][j];
          }
      }
    }

    void Mapping2D::assembleRHS(basics::Vector& vec)
    {
      for( map<int,pair<Range2D,Element2D> >::iterator iter=m_mapping.begin(); iter != m_mapping.end(); ++iter) {
        Element2D& element = iter->second.second;
        Range2D& nodes = iter->second.first;
        basics::Vector rhs = element.width()/2.f*element.height()/2.f*kron(basics::Matrix::diag(element.getWeightX()),basics::Matrix::diag(element.getWeightY()))*element.getRHS();
        for( int i=0;i<element.getRHS().length();++i )
          vec[nodes[0][i]] += rhs[i];
      }
    }

    void Mapping2D::assembleGrid(basics::Vector& vec)
    {
    }

    void Mapping2D::print()
    {
      for( map<int,pair<Range2D,Element2D> >::iterator iter=m_mapping.begin(); iter != m_mapping.end(); ++iter ) {
        cout << "Element " << iter->first << ":" << endl;
        iter->second.first.print();
      }
    }
  }
}

