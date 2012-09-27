PROJECT(mnl)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

ENABLE_LANGUAGE(CXX)

SET(SOURCES function.cpp
            geometry.cpp
            gll.cpp
            gordonhall.cpp
            grid.cpp
            hdf5.cpp
            keyboard.cpp
            mapping.cpp
            matrices.cpp
            matrix.cpp
            memtracker.cpp
            range.cpp
            timer.cpp
            util.cpp
            vector.cpp)

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
  SET(SOURCES ${SOURCES} fftw.cpp)
ENDIF(USE_FFTW)

ADD_LIBRARY(mnl ${SOURCES})
