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
find_package(ROOT REQUIRED COMPONENTS python mathmore minuit)
message(STATUS "Found ROOT v${ROOT_VERSION} at ${ROOT_DIR}")
set(ROOTCLING "${ROOT_rootcling_CMD}")

# Build Headers
include_directories(${FEMTOFITTER_INCLUDE_DIR})
configure_file(
  "src/femtofitter.hpp.in"
  "${FEMTOFITTER_HEADER_DEST}/femtofitter.hpp")

set(SRCS
  "src/r")

add_library(FemtoFitter	SHARED "src/rect.cpp")

