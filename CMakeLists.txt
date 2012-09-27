PROJECT(thesis)

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
  FIND_PACKAGE(HDF5 REQUIRED COMPONENTS CXX)
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

INCLUDE_DIRECTORIES(${HDF5_INCLUDES}
                    ${PROJECT_SOURCE_DIR}
                    ${LAPACK_INCLUDES}
                    ${FFTW_INCLUDES})

ADD_SUBDIRECTORY(mnl)
ADD_SUBDIRECTORY(bm)

ADD_EXECUTABLE(benard-marangoni app/benard-marangoni.cpp)
TARGET_LINK_LIBRARIES(benard-marangoni bm mnl 
                                       ${LAPACK_LIBRARIES}
                                       ${HDF5_LIBRARIES}
                                       ${HDF5_CXX_LIBRARIES})

ADD_EXECUTABLE(benard-marangoni-deformed app/benard-marangoni-deformed.cpp)
TARGET_LINK_LIBRARIES(benard-marangoni-deformed bm mnl 
                                                ${LAPACK_LIBRARIES}
                                                ${HDF5_LIBRARIES}
                                                ${HDF5_CXX_LIBRARIES})

ADD_EXECUTABLE(benard-marangoni-oldschool app/benard-marangoni-oldschool.cpp)
TARGET_LINK_LIBRARIES(benard-marangoni-oldschool bm mnl 
                                                 ${LAPACK_LIBRARIES}
                                                 ${HDF5_LIBRARIES}
                                                 ${HDF5_CXX_LIBRARIES})


ADD_EXECUTABLE(benard-marangoni-periodic app/benard-marangoni-periodic.cpp)
TARGET_LINK_LIBRARIES(benard-marangoni-periodic bm mnl 
                                                ${LAPACK_LIBRARIES}
                                                ${HDF5_LIBRARIES}
                                                ${HDF5_CXX_LIBRARIES})

ADD_EXECUTABLE(benard-marangoni-unsplit app/benard-marangoni-unsplit.cpp)
TARGET_LINK_LIBRARIES(benard-marangoni-unsplit bm mnl 
                                               ${LAPACK_LIBRARIES}
                                               ${HDF5_LIBRARIES}
                                               ${HDF5_CXX_LIBRARIES})

ADD_EXECUTABLE(burgers-2d app/burgers-2d.cpp)
TARGET_LINK_LIBRARIES(burgers-2d bm mnl 
                                 ${LAPACK_LIBRARIES}
                                 ${HDF5_LIBRARIES}
                                 ${HDF5_CXX_LIBRARIES})

ADD_EXECUTABLE(convection-diffusion app/convection-diffusion.cpp)
TARGET_LINK_LIBRARIES(convection-diffusion bm mnl 
                                           ${LAPACK_LIBRARIES}
                                           ${HDF5_LIBRARIES}
                                           ${HDF5_CXX_LIBRARIES})

ADD_EXECUTABLE(convection-diffusion-deformed app/convection-diffusion-deformed.cpp)
TARGET_LINK_LIBRARIES(convection-diffusion-deformed bm mnl 
                                                    ${LAPACK_LIBRARIES}
                                                    ${HDF5_LIBRARIES}
                                                    ${HDF5_CXX_LIBRARIES})

ADD_EXECUTABLE(heatequation-2d app/heatequation-2d.cpp)
TARGET_LINK_LIBRARIES(heatequation-2d bm mnl 
                                 ${LAPACK_LIBRARIES}
                                 ${HDF5_LIBRARIES})

ADD_EXECUTABLE(laplace app/laplace.cpp)
TARGET_LINK_LIBRARIES(laplace bm mnl 
                              ${LAPACK_LIBRARIES}
                              ${HDF5_LIBRARIES}
                              ${HDF5_CXX_LIBRARIES})

ADD_EXECUTABLE(laplace-deformed app/laplace-deformed.cpp)
TARGET_LINK_LIBRARIES(laplace-deformed bm mnl 
                                       ${LAPACK_LIBRARIES}
                                       ${HDF5_LIBRARIES}
                                       ${HDF5_CXX_LIBRARIES})

ADD_EXECUTABLE(laplace-2d app/laplace-2d.cpp)
TARGET_LINK_LIBRARIES(laplace-2d bm mnl 
                                 ${LAPACK_LIBRARIES}
                                 ${HDF5_LIBRARIES}
                                 ${HDF5_CXX_LIBRARIES})

ADD_EXECUTABLE(navierstokes-pc app/navierstokes-pc.cpp)
TARGET_LINK_LIBRARIES(navierstokes-pc bm mnl 
                                      ${LAPACK_LIBRARIES}
                                      ${HDF5_LIBRARIES}
                                      ${HDF5_CXX_LIBRARIES})

ADD_EXECUTABLE(rayleigh-benard app/rayleigh-benard.cpp)
TARGET_LINK_LIBRARIES(rayleigh-benard bm mnl 
                                      ${LAPACK_LIBRARIES}
                                      ${HDF5_LIBRARIES}
                                      ${HDF5_CXX_LIBRARIES})

#ADD_EXECUTABLE(stokes-2d app/stokes-2d.cpp)
#TARGET_LINK_LIBRARIES(stokes-2d bm mnl 
# ${LAPACK_LIBRARIES}
#  ${HDF5_LIBRARIES}
#  ${HDF5_CXX_LIBRARIES})

#ADD_EXECUTABLE(stokes app/stokes.cpp)
#TARGET_LINK_LIBRARIES(stokes bm mnl 
#${LAPACK_LIBRARIES}
#${HDF5_LIBRARIES}
#${HDF5_CXX_LIBRARIES})

ADD_EXECUTABLE(unsteadystokes app/unsteadystokes.cpp)
TARGET_LINK_LIBRARIES(unsteadystokes bm mnl 
                                     ${LAPACK_LIBRARIES}
                                     ${HDF5_LIBRARIES}
                                     ${HDF5_CXX_LIBRARIES})

ADD_EXECUTABLE(unsteadystokes-pc app/unsteadystokes-pc.cpp)
TARGET_LINK_LIBRARIES(unsteadystokes-pc bm mnl 
                                        ${LAPACK_LIBRARIES}
                                        ${HDF5_LIBRARIES}
                                        ${HDF5_CXX_LIBRARIES})

ADD_EXECUTABLE(unsteadystokes-2d app/unsteadystokes-2d.cpp)
TARGET_LINK_LIBRARIES(unsteadystokes-2d bm mnl 
                                        ${LAPACK_LIBRARIES}
                                        ${HDF5_LIBRARIES}
                                        ${HDF5_CXX_LIBRARIES})

ADD_EXECUTABLE(unsteadystokes-pc-2d app/unsteadystokes-pc-2d.cpp)
TARGET_LINK_LIBRARIES(unsteadystokes-pc-2d bm mnl 
                                           ${LAPACK_LIBRARIES}
                                           ${HDF5_LIBRARIES}
                                           ${HDF5_CXX_LIBRARIES})
