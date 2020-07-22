cmake_minimum_required(VERSION 3.0)
### Initial definition of cmake variables
set(CMAKE_INSTALL_PREFIX $ENV{HOME} CACHE PATH "install prefix")
set(CMAKE_BUILD_TYPE Debug CACHE STRING "build type")
set(CMAKE_CXX_FLAGS_DEBUG "-g -W -Wall" CACHE STRING "CXX_FLAGS for debug")
set(CMAKE_C_FLAGS_DEBUG "-g -W -Wall" CACHE STRING "C_FLAGS for debug")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -W -Wall" CACHE STRING "CXX_FLAGS for release")
set(CMAKE_C_FLAGS_RELEASE "-O3 -W -Wall" CACHE STRING "C_FLAGS for release")
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

### Definition of project
project(ANKOK)
set(CMAKE_CXX_STANDARD 14)

### External software ###
### ComptonSoft ###
if(NOT DEFINED CS_INSTALL)
  if(DEFINED ENV{CS_INSTALL})
    set(CS_INSTALL $ENV{CS_INSTALL})
  else()
    set(CS_INSTALL $ENV{HOME})
  endif()
endif(NOT DEFINED CS_INSTALL)

set(CMAKE_MODULE_PATH ${CS_INSTALL}/lib/comptonsoft ${CMAKE_MODULE_PATH})
include(UseComptonSoft)

message("-- Install destination prefix: ${CMAKE_INSTALL_PREFIX}")

### Main ###
set(TARGET main)

include_directories(
  /usr/local/include
  $ENV{HOME}/include
  ${Boost_INCLUDE_DIRS}
  ${G4_INC_DIR}
  ${ROOT_INC_DIR}
  ${ANLG4_INC_DIR}
  ${CS_INC_DIR}
  ${ROOT_INC_DIR}/Minuit2
  ${ROOT_USE_FILE}
  /Users/takashima/experiment/ankok/analysis/DualPhaseArTPCAnalyzer/include
  )

link_directories(
  ${ROOT_LIB_DIR}
  ${ANLG4_LIB_DIR}
  )

add_executable(${TARGET}
  ${TARGET}.cpp
  /Users/takashima/experiment/ankok/analysis/DualPhaseArTPCAnalyzer/src/PulseSignalArray.cc
  /Users/takashima/experiment/ankok/analysis/DualPhaseArTPCAnalyzer/src/DPArTPCEventAnalysis.cc
  /Users/takashima/experiment/ankok/analysis/DualPhaseArTPCAnalyzer/src/S2Function.cc
  /Users/takashima/experiment/ankok/analysis/DualPhaseArTPCAnalyzer/src/S2LossFCN.cc
  /Users/takashima/experiment/ankok/analysis/DualPhaseArTPCAnalyzer/src/S2FitParHolder.cc
  /Users/takashima/experiment/ankok/analysis/DualPhaseArTPCAnalyzer/src/Visualizer.cc
  )

target_link_libraries(${TARGET}
  ${Boost_LIBRARIES}
  ${ROOT_LIB}
  ${G4_LIB}
  ${ANLG4_LIB}
#  ${CS_LIB}
  ${ROOT_LIB_DIR}/libMinuit2.so
  )

install(TARGETS ${TARGET}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

### END
