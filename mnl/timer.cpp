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

#include "timer.h"

#include <sstream>

#ifdef OPENMP
#include <omp.h>
#endif

#ifdef HAS_MPI
#include <mpi.h>
#endif

namespace mnl {
  namespace utilities {

    Timer::Timer(const std::string& n_name) : m_elapse(0), m_clockstart(0), m_running(false), m_name(n_name)
    {
    }

    Timer::~Timer()
    {
    }

    void Timer::start()
    {
      m_running = true;
#ifdef HAS_MPI
      m_clockstart = MPI_Wtime();
#elif defined(OPENMP)
      m_clockstart = omp_get_wtime();
#else
      m_clockstart = clock();
#endif
    }

    void Timer::pause()
    {
      if( m_running ) {
#ifdef HAS_MPI
        m_elapse += MPI_Wtime()-m_clockstart;
#elif defined(OPENMP)
        m_elapse += omp_get_wtime()-m_clockstart;
#else
        m_elapse += clock()-m_clockstart;
#endif
        m_running = false;
      }			
    }

    std::string Timer::report()
    {
      std::stringstream result;
      result << m_name << " ran for: ";
#if defined(HAS_MPI) || defined(OPENMP)
      result << static_cast<Real>(m_elapse) << "s.";
#else
      result << static_cast<Real>(m_elapse)/CLOCKS_PER_SEC << "s.";
#endif
      return( result.str() );
    }

    Profiler::Profiler() : m_running(false)
    {
      m_timer.push_back(Timer("Global timer"));
    }

    Profiler::~Profiler()
    {
    }

    int Profiler::add(const std::string& n_name)
    {
      m_timer.push_back(Timer(n_name));
      return( m_timer.size()-1 );
    }

    void Profiler::start(int handle)
    {
      //assert( handle < m_timer.size() );

      if( !m_running ) {
        m_running = true;
        m_timer[0].start(); // start global timer
      }
      m_timer[handle].start();
    }

    void Profiler::pause(int handle)
    {
      //assert( handle < m_timer.size() );

      m_timer[handle].pause();
    }

    std::string Profiler::report()
    {
      m_timer[0].pause(); // pause global timer
      m_running = false;
      std::stringstream result;
      Real total=0.f;
      for( int i=1;i<m_timer.size();++i )
        total += m_timer[i].elapsed();

      result << "We used " << (Real)m_timer[0].elapsed() << "s totally." << std::endl;
      result << "Out of this, timers were active " << total << "s." << std::endl;
      for( int i=1;i<m_timer.size();++i )
        result << m_timer[i].report() << " (" << m_timer[i].elapsed()/total*100.f << "%)" << std::endl;

      return( result.str() );
    }

    Profiler g_profiler;
  }
}

