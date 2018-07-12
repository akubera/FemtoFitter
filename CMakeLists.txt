#
# CMakeLists.txt
#

cmake_minimum_required(VERSION 3.4)

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

set(INCS
  ${CMAKE_CURRENT_SOURCE_DIR}/femtofitter/inc/fitter.h
)

ROOT_GENERATE_DICTIONARY(femfit_dict ${INCS})

set(ROOT_LIBS
  ROOT::Hist
  ROOT::Minuit
)

add_library(FemtoFitter	SHARED femfit_dict)
target_link_libraries(FemtoFitter ${ROOT_LIBS})
