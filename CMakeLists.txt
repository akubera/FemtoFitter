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
set(FEMTOFITTER_INCLUDE_DIR "${CMAKE_BINARY_DIR}/include")
set(FEMTOFITTER_HEADER_DEST "${FEMTOFITTER_INCLUDE_DIR}/femtofitter")
set(BUILD_DATE "March 11, 2017")

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
  libCore.so
)

add_library(FemtoFitter	SHARED femfit_dict)
target_link_libraries(FemtoFitter ${ROOT_LIBS})

