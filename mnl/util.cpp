#include "util.h"

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <algorithm>

#ifdef HAS_MPI
#include <mpi.h>
#endif

using namespace std;

namespace mnl {
  namespace utilities {
    int sendmail(const string& mailaddy, const string& mailsubject,
        const string& mailbody)
    {
      return system(("echo \"Simulation done @ `hostname`\n"+mailbody+"\" | mail -s"+mailsubject+" "+mailaddy).c_str());
    }

    iterationStatistics::iterationStatistics(int entries, int size,
        vector<int*>& scount,
        vector<int*>& sdispl,
        vector<int*>& rcount,
        vector<int*>& rdispl)
      : m_entries(entries)
    {
      m_iters = new int[entries];
      reset();

      m_scount = new int[size];
      m_sdispl = new int[size];
      m_rcount = new int[size];
      m_rdispl = new int[size];
      for( int i=0;i<2;++i ) {
        scount.push_back(new int[size]);
        sdispl.push_back(new int[size]);
        rcount.push_back(new int[size]);
        rdispl.push_back(new int[size]);
      }
    }

    iterationStatistics::~iterationStatistics()
    {
      delete[] m_iters;
      delete[] m_sdispl;
      delete[] m_scount;
      delete[] m_rdispl;
      delete[] m_rcount;
    }

    bool sort_planes(const pair<int,int>& left, const pair<int,int>& right)
    {
      return left.first > right.first;
    }

    vector<int> iterationStatistics::getCounts(int size) const
    {
      //            vector<int> result;
      //            int meanplanes = m_entries/size;
      //            for( int i=0;i<size;++i )
      //                result.push_back(meanplanes);
      //            for( int i=0;i<m_entries%size;++i )
      //                result[i]++;
      vector<int> result;
      int sum = accumulate(m_iters,m_iters+m_entries,0);
      double mean = double(sum)/size;
      int iCurr = 0;
      int j=0;
      int oldj=0;
      for( int i=0;i<size;++i) {
        while( j < m_entries && iCurr+m_iters[j] < (i+1)*mean )
          iCurr += m_iters[j++];
        //                if( j < m_entries && iCurr+m_iters[j]-(i+1)*mean < (i+1)*mean-iCurr )
        //                    ++j;
        if( j < m_entries-1 && j == oldj )
          iCurr += m_iters[j++];
        if( i == size-1 )
          j = m_entries;

        if( i < size-1 && j == m_entries)
          j--;

        result.push_back(j-oldj);
        oldj = j;
      }


      return( result );
    }

    pair< vector<int>,vector< vector<int> > > iterationStatistics::getCounts2(int size) const
    {
      pair< vector<int>,vector< vector<int> > > result;
      int sum = accumulate(m_iters,m_iters+m_entries,0);
      vector< pair<int,int> > sorted;
      for( int i=0;i<m_entries;++i )
        sorted.push_back(make_pair(m_iters[i],i));
      sort(sorted.begin(),sorted.end(),sort_planes);
      double mean = double(sum)/size;
      int iCurr = 0;
      vector<int> dummy;

      int front=0;
      int back=m_entries-1;
      for( int i=0;i<size;++i ) {
        result.second.push_back(dummy);
        iCurr = 0;
        int planes=0;
        while( (front < back && iCurr+sorted[front].first < mean || !planes ) && back-front > size-i) {
          result.second[i].push_back(sorted[front].second);
          iCurr += sorted[front++].first;
          planes++;
        }
        while( (front < back && iCurr < mean || !planes) && back-front > size-i ) {
          if( iCurr+sorted[back].first > mean )
            if( iCurr+sorted[back].first-mean > mean-iCurr )
              break;
          result.second[i].push_back(sorted[back].second);
          iCurr += sorted[back--].first;
          planes++;
        }
        if( i == size-1 ) {
          while( front <= back ) {
            result.second[i].push_back(sorted[back--].second);
            planes++;
          }
        }
        result.first.push_back(planes);
      }
      return( result );
    }

    vector<int> iterationStatistics::getStarts(int size) const
    {
      vector<int> counts = getCounts(size);
      vector<int> result;
      result.push_back(0);
      for( int i=1;i<size;++i )
        result.push_back(result[i-1]+counts[i-1]);
      return( result );
    }

    void iterationStatistics::getDisplacements(vector<int*>& scount, 
        vector<int*>& sdispl,
        vector<int*>& rcount,
        vector<int*>& rdispl,
        int rank, int size, int mult)
    {
      vector<int> c = getCounts(size);
      getDisplacementsInt(scount,sdispl,rcount,rdispl,rank,size,mult,c);
    }

    void iterationStatistics::getDisplacements2(vector<int*>& scount, 
        vector<int*>& sdispl,
        vector<int*>& rcount,
        vector<int*>& rdispl,
        int rank, int size, int mult)
    {
      pair<vector<int>,vector< vector<int> > > c = getCounts2(size);
      getDisplacementsInt(scount,sdispl,rcount,rdispl,rank,size,mult,c.first);
    }

    void iterationStatistics::getDisplacementsInt(vector<int*>& scount, 
        vector<int*>& sdispl,
        vector<int*>& rcount,
        vector<int*>& rdispl,
        int rank, int size,
        int mult, vector<int>& c)
    {
      vector<int> s = getStarts(size);
      for( int i=0;i<size;++i ) {
        m_rcount[i] = c[i];
        m_rdispl[i] = s[i];
      }
      for( int i=0;i<size;++i ) {
        m_scount[i] = m_rcount[rank];
        m_sdispl[i] = m_rdispl[rank];
      }

      /* forward send */
      sdispl[0][0] = rdispl[0][0] = 0;
      scount[0][0] = m_rcount[0]*mult;
      rcount[0][0] = m_rcount[rank]*mult;
      for( int i=1;i<size;++i ) {
        scount[0][i] = m_rcount[i]*mult;
        sdispl[0][i] = sdispl[0][i-1]+scount[0][i-1];
        rcount[0][i] = rcount[0][i-1];
        rdispl[0][i] = rdispl[0][i-1]+rcount[0][i-1];
      }

      /* backward send */
      sdispl[1][0] = rdispl[1][0] = 0;
      scount[1][0] = m_rcount[rank]*mult;
      rcount[1][0] = m_rcount[0]*mult;
      for( int i=1;i<size;++i ) {
        scount[1][i] = m_rcount[rank]*mult;
        sdispl[1][i] = sdispl[1][i-1]+scount[1][i-1];
        rcount[1][i] = m_rcount[i]*mult;
        rdispl[1][i] = rdispl[1][i-1]+rcount[1][i-1];
      }
    }

    void iterationStatistics::exchange()
    {
#ifdef HAS_MPI
      int* iters2 = new int[m_entries];
      memcpy(iters2,m_iters,m_entries*sizeof(int));
      MPI_Allreduce(iters2,m_iters,m_entries,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      delete[] iters2;
#endif
    }

    void iterationStatistics::cleanDisplacements(vector<int*>& scount, 
        vector<int*>& sdispl,
        vector<int*>& rcount,
        vector<int*>& rdispl)
    {
      for( int i=0;i<scount.size();++i ) {
        delete[] scount[i];
        delete[] sdispl[i];
        delete[] rcount[i];
        delete[] rdispl[i];
      }
      scount.clear();
      sdispl.clear();
      rcount.clear();
      rdispl.clear();
    }
  };
};
