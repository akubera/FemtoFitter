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
message(STATUS "cmake path: ${CMAKE_PREFIX_PATH}")

list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/cmake/modules)
message(STATUS "cmake module path: ${CMAKE_MODULE_PATH}")

#
# Dependencies
#
find_package(ROOT REQUIRED COMPONENTS PyROOT MathMore Minuit vdt)
include(${ROOT_USE_FILE})
message(STATUS "Found ROOT v${ROOT_VERSION} at ${ROOT_DIR}")
set(ROOTCLING "${ROOT_rootcling_CMD}")


# include_directories(${FEMTOFITTER_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR})
include_directories(src)

message(STATUS "VDT: ${VDT_INCLUDE_DIRS}")

if(NOT TARGET Vdt::Vdt)
  add_library(Vdt::Vdt INTERFACE IMPORTED)
  set_property(TARGET Vdt::Vdt PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${VDT_INCLUDE_DIRS}")
  set_property(TARGET Vdt::Vdt PROPERTY INTERFACE_LINK_LIBRARIES "${VDT_LIBRARIES}")
endif(NOT TARGET Vdt::Vdt)

root_generate_dictionary(
  fitinspector_dict
  src/bin/fit-inspector/MainFrame.hpp
  LINKDEF src/bin/fit-inspector/LinkDef.h
  MODULE fitinspector_dict)


include_directories(
  /usr/include/python3.6m
  /usr/include/x86_64-linux-gnu/python3.6m
)


# find_package(Vdt)
# find_package()

add_library(fitinspector STATIC fitinspector_dict src/bin/fit-inspector/MainFrame.cpp)
add_executable(fit-inspector src/bin/fit-inspector/main.cc)
target_link_libraries(fit-inspector fitinspector
    ROOT::Hist
    ROOT::Gui
    ROOT::PyROOT
    # ${ROOT_vdt_LIBRARY}
    # ROOT::vectorDict
    # ROOT::Math
    # Vdt::Vdt
    /usr/lib/x86_64-linux-gnu/libpython3.6m.so.1.0

    # ROOT::VecOps
    # ROOTVecOps
    # ${VDT_LIBRARIES}
    )

    return(0)

# Build Headers
# set(INCS
#   ${CMAKE_CURRENT_SOURCE_DIR}/femtofitter/inc/fitter.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/femtofitter/MomentumResolutionCorrector.hpp
# )

set(INCS
  src/Data1D.hpp
  src/Data3D.hpp
  fitter/FitterGaussOSL.hpp
  # fitter/FitterGaussOSL3.hpp
  # fitter/FitterGaussFull.hpp
  # fitter/FitterLevy3D.hpp
  # fitter/FitterLevyFull.hpp
  fitter/Levy1D.hh
  fsi/FsiKFile.hpp
  fsi/FsiGamov.hpp
  fsi/FsiStatic.hpp
  src/CalculatorFsi.hpp
  src/ParamHints.hpp
)
list(TRANSFORM INCS PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -ffast-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -fopenmp")

set(SRCS
  src/Data3D.cpp
  src/Data1D.cpp
  src/Value.cpp
  fsi/FsiKFile.cpp
  src/CoulombHist.cpp
  src/MinimizerMinuit.cpp
  src/ParamHints.cpp
#   ${CMAKE_CURRENT_SOURCE_DIR}/femtofitter/MomentumResolutionCorrector.cxx
)

ROOT_GENERATE_DICTIONARY(femfit_dict ${INCS}
                         LINKDEF Linkdef.h
                         MODULE FemtoFitter)
# set_target_properties(femfit_dict PROPERTIES)

set(ROOT_LIBS
  ROOT::Hist
  ROOT::Minuit
)

find_library(LIB_FEMTO
             PWGCFfemtoscopy
             HINTS $ENV{ALICE_PHYSICS}/lib)

add_library(FemtoFitter	SHARED femfit_dict ${SRCS})
target_link_libraries(FemtoFitter ${ROOT_LIBS})

add_executable(run-fit-spec	src/bin/run-fit-spec.cc)
target_include_directories(run-fit-spec PRIVATE src)
target_link_libraries(run-fit-spec ROOT::Minuit)

add_executable(quick-fit src/bin/quick-fit.cc)
target_include_directories(quick-fit PRIVATE src)
target_link_libraries(quick-fit FemtoFitter ${LIB_FEMTO})
