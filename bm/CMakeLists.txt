PROJECT(bm)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

ENABLE_LANGUAGE(CXX)

SET(SOURCES bigcircle.cpp
            consistent.cpp
            convection.cpp
            divergence.cpp
            gradient.cpp
            laplacian.cpp
            legendrelegendre.cpp
            legendrelegendrew.cpp
            mass.cpp
            quadratic.cpp
            poissonsolver-fe.cpp
            poissonsolver-lg.cpp
            poissonsolver-lgw.cpp
            sem.cpp
            uzawa.cpp)

IF(NOT HDF5_LIBRARIES)
  FIND_PACKAGE(HDF5 REQUIRED)
ENDIF(NOT HDF5_LIBRARIES)

IF(NOT LAPACK_LIBRARIES)
  FIND_PACKAGE(LAPACK REQUIRED)
ENDIF(NOT LAPACK_LIBRARIES)

IF(USE_MPI)
  FIND_PACKAGE(MPI REQUIRED)
ENDIF(USE_MPI)

IF(NOT OPENMP_FOUND)
  FIND_PACKAGE(OpenMP)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DOPENMP")
ENDIF(NOT OPENMP_FOUND)

IF(USE_FFTW)
  FIND_PACKAGE(FFTW REQUIRED)
ENDIF(USE_FFTW)

ADD_LIBRARY(bm ${SOURCES})
