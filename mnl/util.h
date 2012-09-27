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
#ifndef MNL_UTIL_H_
#define MNL_UTIL_H_

#include <vector>
#include <string>
#include <numeric>

namespace mnl {
  namespace utilities {
    int sendmail(const std::string& mailaddy, const std::string& mailsubject, const std::string& mailbody);

    class iterationStatistics {
      public:
        iterationStatistics(int entries, int size,
            std::vector<int*>& scount,
            std::vector<int*>& sdispl,
            std::vector<int*>& rcount,
            std::vector<int*>& rdispl);
        ~iterationStatistics();

        inline int* get()
        {
          return m_iters;
        }

        void reset()
        {
          for( int i=0;i<m_entries;++i )
            m_iters[i] = 1;
        }

        std::vector<int> getStarts(int size) const;
        std::vector<int> getCounts(int size) const;
        std::pair< std::vector<int>,std::vector< std::vector<int> > > getCounts2(int size) const;
        void getDisplacements(std::vector<int*>& sdispl,
            std::vector<int*>& scount,
            std::vector<int*>& rdispl,
            std::vector<int*>& rcount,
            int rank, int size, int mult);
        void getDisplacements2(std::vector<int*>& sdispl,
            std::vector<int*>& scount,
            std::vector<int*>& rdispl,
            std::vector<int*>& rcount,
            int rank, int size, int mult);

        void getDisplacementsInt(std::vector<int*>& sdispl,
            std::vector<int*>& scount,
            std::vector<int*>& rdispl,
            std::vector<int*>& rcount,
            int rank, int size, int mult,
            std::vector<int>& c);
        void cleanDisplacements(std::vector<int*>& sdispl,
            std::vector<int*>& scount,
            std::vector<int*>& rdispl,
            std::vector<int*>& rcount);
        void exchange();
        int  m_entries;
      protected:
        int* m_iters;
        int* m_sdispl;
        int* m_scount;
        int* m_rdispl;
        int* m_rcount;
    };
  };
};

#endif // MNL_UTIL_H_

