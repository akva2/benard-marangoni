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
#ifndef MNL_TIMER_H_
#define MNL_TIMER_H_

#include "config.h"

#include <vector>
#include <string>
#include <iostream>
#include <ctime>

namespace mnl {
  namespace utilities {
    class Timer {
      public:
        Timer(const std::string& n_name);
        ~Timer();

        void start();
        void pause();

        inline Real elapsed()
        {
#if defined(OPENMP) || defined(HAS_MPI)
          return( static_cast<Real>(m_elapse) );
#else
          return( static_cast<Real>(m_elapse)/CLOCKS_PER_SEC );
#endif
        }

        std::string report();

      protected:
#if defined(OPENMP) || defined(HAS_MPI)
        double m_clockstart;
        double m_elapse;
#else
        clock_t m_clockstart;
        clock_t m_elapse;
#endif
        bool m_running;
        std::string m_name;
    };

    class Profiler {
      public:
        Profiler();
        ~Profiler();

        int add(const std::string& n_name); // returns handle to timer

        void start(int handle);
        void pause(int handle);

        std::string report();

      protected:
        bool m_running;
        std::vector<Timer> m_timer; // global in front!
    };

    extern Profiler g_profiler;
  }
}

#endif

