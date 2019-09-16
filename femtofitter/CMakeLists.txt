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
find_package(Python3 COMPONENTS Development)

include(${ROOT_USE_FILE})
message(STATUS "Found ROOT v${ROOT_VERSION} at ${ROOT_DIR}")
set(ROOTCLING "${ROOT_rootcling_CMD}")


# include_directories(${FEMTOFITTER_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR})
include_directories(src)

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

include_directories(${Python3_INCLUDE_DIRS})

# Build Headers
set(INCS
  src/Data1D.hpp
  src/Data3D.hpp
  fitter/Fitter1DGauss.hpp
  fitter/Fitter1DGaussLin.hpp
  fitter/Fitter1DGaussPolyBg.hpp
  fitter/Fitter1DLevy.hpp
  fitter/Fitter1DLevyPolyBg.hpp
  fitter/Fitter3DGaussLcms.hpp
  fitter/Fitter3DGaussLcmsOS.hpp
  fitter/Fitter3DGaussLcmsOL.hpp
  fitter/FitterGaussFull.hpp
  fitter/FitterLevy3D.hpp
  fitter/FitterLevyFull.hpp
  fsi/FsiKFile.hpp
  fsi/FsiGamov.hpp
  fsi/FsiStatic.hpp
  mrc/Mrc3DRatio.hpp
  mrc/Mrc3DRatioMixed.hpp
  mrc/Mrc1DRatio.hpp
  mrc/Mrc1DRatioMixed.hpp
  mrc/Mrc1DMatrix.hpp
  mrc/MrcTransform1D.hpp
  mrc/MrcHypercube3D.hpp
  src/CalculatorFsi.hpp
  src/ParamHints.hpp
)

list(TRANSFORM INCS PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/)

set(CMAKE_CXX_FLAGS_DEBUG
    "-O1 -g -fPIC"
    CACHE STRING "Debug compile flags" FORCE)

set(CMAKE_CXX_FLAGS_RELEASE
    "-Ofast -march=native -ffast-math -DNDEBUG"
    CACHE STRING "fast flags" FORCE)

set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -Wshadow -pedantic -fopenmp"
    CACHE STRING "std flags" FORCE)

set(SRCS
  src/Data3D.cpp
  src/Data1D.cpp
  src/Value.cpp
  fsi/FsiKFile.cpp
  fitter/Fitter1DGauss.cpp
  mrc/MrcHypercube3D.cpp
  src/CoulombHist.cpp
  src/MinimizerMinuit.cpp
  src/ParamHints.cpp
)

ROOT_GENERATE_DICTIONARY(femfit_dict ${INCS}
                         LINKDEF Linkdef.h
                         OPTIONS -Wno-inconsistent-missing-override
                         MODULE FemtoFitter)

set(ROOT_LIBS
  ROOT::Hist
  ROOT::Gui
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

find_library(LIB_TENDER
             Tender
             HINTS $ENV{ALICE_PHYSICS}/lib)


#add_library(fitinspector
#            STATIC fitinspector_dict
#            src/bin/fit-inspector/MainFrame.cpp
#            src/bin/fit-inspector/PanelFitResult.cpp
#            src/bin/fit-inspector/PanelCutConfig.cpp
#            src/bin/fit-inspector/Projectionist.cpp)
#target_include_directories(fitinspector PUBLIC $ENV{ALICE_PHYSICS}/include)
#add_executable(fit-inspector src/bin/fit-inspector/main.cc)
#target_link_libraries(fit-inspector fitinspector
#    ROOT::Hist
#    ROOT::Gui
#    ROOT::PyROOT
#    Python3::Python
#    ${LIB_FEMTO}
#${LIB_TENDER}
#    )
