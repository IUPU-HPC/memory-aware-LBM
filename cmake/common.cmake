if(CMAKE_BUILD_TYPE MATCHES Debug)
  message("-- ${CMAKE_CURRENT_SOURCE_DIR} > Debug build.")
elseif(CMAKE_BUILD_TYPE MATCHES Release)
  message("-- ${CMAKE_CURRENT_SOURCE_DIR} > Release build.")

elseif(CMAKE_BUILD_TYPE STREQUAL Stampede)
  message("-- ${CMAKE_CURRENT_SOURCE_DIR} > Stampede uninstructed build")
  set(CMAKE_C_COMPILER  icc)
  set(CMAKE_CXX_COMPILER  icpc)
  #set(CMAKE_C_FLAGS "-O3 -Wall -Wextra -xCORE-AVX2 -axMIC-AVX512" CACHE STRING "cflags")
  #set(CMAKE_C_FLAGS "-O3 -Wall -Wextra -xCORE-AVX2 -axCORE-AVX512,MIC-AVX512" CACHE STRING "cflags")
  set(CMAKE_C_FLAGS "-O3 -Wall -Wextra $ENV{TACC_VEC_FLAGS}" CACHE STRING "cflags")
  #set(CMAKE_C_FLAGS "-O3 -Wall -Wextra -xCORE-AVX2" CACHE STRING "cflags")
  #set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -xCORE-AVX2 -axMIC-AVX512" CACHE STRING "cxxflags")
  set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -xCORE-AVX2 -axCORE-AVX512,MIC-AVX512" CACHE STRING "cxxflags")
  set(PAPI_ROOT $ENV{TACC_PAPI_DIR})


elseif(CMAKE_BUILD_TYPE STREQUAL Stampede_TAU)
  message("-- ${CMAKE_CURRENT_SOURCE_DIR} > Stampede instructed build")

  set(CMAKE_C_COMPILER  tau_cc.sh)
  set(CMAKE_CXX_COMPILER  tau_cxx.sh)
  set(CMAKE_C_FLAGS "-O3 -Wall -Wextra -xCORE-AVX2 -axMIC-AVX512 -fPIC" CACHE STRING "cflags")
  set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -xCORE-AVX2 -axMIC-AVX512 -fPIC" CACHE STRING "cxxflags")

elseif(CMAKE_BUILD_TYPE STREQUAL Bridges)
  set(CMAKE_C_COMPILER  icc)
  set(CMAKE_CXX_COMPILER  icpc)
  message("-- ${CMAKE_CURRENT_SOURCE_DIR} > Bridges uninstructed build")
  set(CMAKE_C_FLAGS "-O3 " CACHE STRING "cflags")
  set(CMAKE_CXX_FLAGS "-O3 " CACHE STRING "cxxflags")
  set(PAPI_ROOT $ENV{PAPI_DIR})

elseif(CMAKE_BUILD_TYPE STREQUAL Bridges_TAU)
  message("-- ${CMAKE_CURRENT_SOURCE_DIR} > Bridges instructed build with tau")
  set(CMAKE_C_COMPILER  tau_cc.sh)
  set(CMAKE_CXX_COMPILER  tau_cxx.sh)
  set(CMAKE_C_FLAGS "-fPIC ${ADD_FLAGS}" CACHE STRING "cflags")
  set(CMAKE_CXX_FLAGS "-fPIC ${ADD_FLAGS}" CACHE STRING "cxxflags")

endif()

message("-- Including transport method in ${TRANSPORT_LIB} is loaded...")
message("-- build type not set")
 

#set(CMAKE_CXX_STANDARD 14)
#set(GCC_COVERAGE_COMPILE_FLAGS "-fPIC -msse3")

#set(CMAKE_CXX_FLAGS_DEBUG "-O2 -g")
#set(CMAKE_C_FLAGS_DEBUG "-O2 -g")
