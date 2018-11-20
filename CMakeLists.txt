#
# CMakeLists.txt
#

cmake_minimum_required(VERSION 3.12)

project(FEMTOFITTER
        VERSION 0.1.0
        LANGUAGES CXX)

#
# Project Settings
#
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

set(FEMTOFITTER_INCLUDE_DIR "${CMAKE_BINARY_DIR}/include")
set(FEMTOFITTER_HEADER_DEST "${FEMTOFITTER_INCLUDE_DIR}/femtofitter")
EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE BUILD_DATE)

option(TARGET_DOCUMENTATION "Add documentation target - requires Doxygen" OFF)

if (TARGET_DOCUMENTATION)
  find_package(Doxygen REQUIRED)
  set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile)
  add_custom_target(
    doxygen-xml
      COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile_in}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/docs
      COMMENT "Build documentation XML with Doxygen"
      VERBATIM)
endif(TARGET_DOCUMENTATION)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#
# Dependencies
#
find_package(ROOT REQUIRED COMPONENTS PyROOT MathMore Minuit)
include(${ROOT_USE_FILE})
message(STATUS "Found ROOT v${ROOT_VERSION} at ${ROOT_DIR}")
set(ROOTCLING "${ROOT_rootcling_CMD}")

# Build Headers
include_directories(${FEMTOFITTER_INCLUDE_DIR})

# set(INCS
#   ${CMAKE_CURRENT_SOURCE_DIR}/femtofitter/inc/fitter.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/femtofitter/MomentumResolutionCorrector.hpp
# )

set(INCS
  src/FitterGaussOSL.hpp
  src/FitterGaussFull.hpp
  src/FitterLevy.hpp
)
list(TRANSFORM INCS PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -ffast-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -fopenmp")

set(SRCS
  src/Data3D.cpp
  src/Data1D.cpp
  src/Fitter.cpp
  src/CoulombHist.cpp
#   ${CMAKE_CURRENT_SOURCE_DIR}/femtofitter/MomentumResolutionCorrector.cxx
)

ROOT_GENERATE_DICTIONARY(femfit_dict ${INCS} LINKDEF Linkdef.h)
# set_target_properties(femfit_dict PROPERTIES)

set(ROOT_LIBS
  ROOT::Hist
  ROOT::Minuit
)

add_library(FemtoFitter	SHARED femfit_dict ${SRCS})
target_link_libraries(FemtoFitter ${ROOT_LIBS})

add_executable(run-fit-spec	src/bin/run-fit-spec.cc)
target_include_directories(run-fit-spec PRIVATE src)
target_link_libraries(run-fit-spec ROOT::Minuit)

add_executable(quick-fit src/bin/quick-fit.cc)
target_include_directories(quick-fit PRIVATE src)
target_link_libraries(quick-fit FemtoFitter)
