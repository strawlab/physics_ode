cmake_minimum_required(VERSION 2.4.6)
include(FindPkgConfig)
include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)

# Set the build type.  Options are:
#  Coverage       : w/ debug symbols, w/o optimization, w/ code-coverage
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries
set(ROS_BUILD_TYPE Release)

rosbuild_init()

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

rosbuild_find_ros_package(opende)

##############################################################################

link_directories(${LIBRARY_OUTPUT_PATH})

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( ${PROJECT_SOURCE_DIR}/include/parallel_quickstep )
include_directories( ${opende_PACKAGE_PATH}/opende/include/src )

set(PARALLEL_QUICKSTEP_FLAGS -O3 )#-DTIMING)# -DVERBOSE -DBENCHMARKING -DERROR )
rosbuild_add_boost_directories()
add_definitions(${PARALLEL_QUICKSTEP_FLAGS})

# default to CPU fall back so everyone can compile this package
set(USE_CPU "1")
#set(USE_CUDA "1")
#set(USE_OPENCL "1")
#set(USE_OPENMP "1")

################################################
# Automatically set USE_CUDA to 1 if it is found
################################################
# CMake 2.8.0 or greater required for built-in CUDA module
if( ${CMAKE_MINOR_VERSION} GREATER 6.2 )
  if( ${CMAKE_MINOR_VERSION} LESS 8 )
    set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules/;${CMAKE_MODULE_PATH}")
  endif( ${CMAKE_MINOR_VERSION} LESS 8 )
  find_package(CUDA QUIET)
  if (CUDA_FOUND)
    set(USE_CPU)
    set(USE_CUDA "1")
    MESSAGE(STATUS "CUDA Found, compiling with CUDA support")
  else (CUDA_FOUND)
    MESSAGE(STATUS "CUDA Not Found, using default CPU")
  endif (CUDA_FOUND)
endif( ${CMAKE_MINOR_VERSION} GREATER 6.2 )

if(DEFINED USE_CUDA)

  # CMake 2.8.0 or greater required for built-in CUDA module
  if( ${CMAKE_MINOR_VERSION} LESS 8 )
    set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules/;${CMAKE_MODULE_PATH}")
  endif( ${CMAKE_MINOR_VERSION} LESS 8 )

  find_package(CUDA QUIET REQUIRED)

  include_directories( ${CUDA_INCLUDE_DIRS} )
  add_definitions(-DUSE_CUDA)

  set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)

  #default to oldest stuff
  set(CUDA_TARGET_SM "sm_10")

  if (NOT DEFINED CUDA_TARGET_SM)
    set(CUDA_TARGET_SM "sm_20")
  endif (NOT DEFINED CUDA_TARGET_SM)

  if (CUDA_TARGET_SM STREQUAL "sm_10")
    set(DOUBLE_SUPPORT_ENABLED FALSE)
    set(ATOMIC_SUPPORT_ENABLED FALSE)
    add_definitions(-DCUDA_SM10)
    message(STATUS "using sm_10 without DOUBLE_SUPPORT_ENABLED and ATOMIC_SUPPORT_ENABLED")
  elseif (CUDA_TARGET_SM STREQUAL "sm_11")
    set(DOUBLE_SUPPORT_ENABLED FALSE)
    set(ATOMIC_SUPPORT_ENABLED TRUE)
    add_definitions(-DCUDA_SM11)
    message(STATUS "using sm_11 with ATOMIC_SUPPORT_ENABLED")
  elseif (CUDA_TARGET_SM STREQUAL "sm_12")
    set(DOUBLE_SUPPORT_ENABLED FALSE)
    set(ATOMIC_SUPPORT_ENABLED TRUE)
    message(STATUS "using sm_12 with ATOMIC_SUPPORT_ENABLED")
    add_definitions(-DCUDA_SM12)
  elseif (CUDA_TARGET_SM STREQUAL "sm_13")
    set(DOUBLE_SUPPORT_ENABLED TRUE)
    set(ATOMIC_SUPPORT_ENABLED TRUE)
    message(STATUS "using sm_13 with DOUBLE_SUPPORT_ENABLED and ATOMIC_SUPPORT_ENABLED")
    add_definitions(-DCUDA_SM13)
  elseif (CUDA_TARGET_SM STREQUAL "sm_20")
    set(DOUBLE_SUPPORT_ENABLED TRUE)
    set(ATOMIC_SUPPORT_ENABLED TRUE)
    message(STATUS "using sm_20 with DOUBLE_SUPPORT_ENABLED and ATOMIC_SUPPORT_ENABLED")
    add_definitions(-DCUDA_SM20)
  else()
    message( FATAL_ERROR "Unknown CUDA_TARGET_SM: '${CUDA_TARGET_SM}', must be one of: 'sm_10' 'sm_11' 'sm_12' 'sm_13' 'sm_20'")
  endif()

  message(STATUS "CUDA Target Architecture: ${CUDA_TARGET_SM}")

  if (DOUBLE_SUPPORT_ENABLED)
    add_definitions(-DCUDA_DOUBLESUPPORT)
  else (DOUBLE_SUPPORT_ENABLED)
    add_definitions(-DCUDA_NODOUBLESUPPORT)
  endif (DOUBLE_SUPPORT_ENABLED)

  if (ATOMIC_SUPPORT_ENABLED)
    add_definitions(-DCUDA_ATOMICSUPPORT)
  else (ATOMIC_SUPPORT_ENABLED)
    add_definitions(-DCUDA_NOATOMICSUPPORT)
  endif (ATOMIC_SUPPORT_ENABLED)

  set(CUDA_SOURCE_FILES src/cuda_kernels.cu src/test_lib.cpp )
  set(CUDA_SOLVER_SOURCE_FILES
    src/parallel_stepper.cpp
    src/parallel_solver.cpp
    src/parallel_batch.cpp
    src/parallel_reduce.cpp
    src/parallel_quickstep.cpp
    src/cuda_solver.cpp)

  set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -use_fast_math -O3 "-arch;${CUDA_TARGET_SM}" ${parallel_quickstep_CFLAGS_OTHER})

  cuda_compile(CUDA_GEN_FILES ${CUDA_SOURCE_FILES} SHARED -fPic)
  rosbuild_add_library( parallel_quickstep
    ${CUDA_GEN_FILES}
    ${CUDA_SOURCE_FILES}
    ${CUDA_SOLVER_SOURCE_FILES}
    )

  rosbuild_add_executable( parallel_quickstep_lib_test src/main_for_lib.cpp )
  rosbuild_add_link_flags(parallel_quickstep -fopenmp)
  target_link_libraries(parallel_quickstep ${CUDA_LIBRARIES})
  target_link_libraries(parallel_quickstep_lib_test parallel_quickstep)
  cuda_build_clean_target()

elseif( DEFINED USE_OPENMP )

  add_definitions(-DUSE_OPENMP -fopenmp)

  set(OPENMP_SOLVER_SOURCE_FILES
    src/parallel_stepper.cpp
    src/parallel_solver.cpp
    src/parallel_batch.cpp
    src/parallel_reduce.cpp
    src/parallel_quickstep.cpp
    src/openmp_solver.cpp
    src/openmp_kernels.cpp )

  rosbuild_add_library( parallel_quickstep
    ${OPENMP_SOLVER_SOURCE_FILES}
    )

  rosbuild_add_link_flags(parallel_quickstep)
  rosbuild_add_openmp_flags(parallel_quickstep)

elseif( DEFINED USE_OPENCL )

  set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules/;${CMAKE_MODULE_PATH}")
  find_package( OpenCL REQUIRED )
  add_definitions(-DUSE_OPENCL -fopenmp)
  #include_directories( ${OPENCL_INCLUDE_DIRS} )
  include_directories( ${CUDA_TOOLKIT_ROOT_DIR}/include )

  set(OPENCL_SOLVER_SOURCE_FILES
    src/opencl_utils.cpp
    src/opencl_solver.cpp
    src/opencl_kernels.cpp
    src/parallel_stepper.cpp
    src/parallel_solver.cpp
    src/parallel_batch.cpp
    src/parallel_reduce.cpp
    src/parallel_quickstep.cpp)

  rosbuild_add_library( parallel_quickstep
    ${OPENCL_SOLVER_SOURCE_FILES}
    )

  target_link_libraries(parallel_quickstep ${OPENCL_LIBRARIES})

  rosbuild_add_executable(parallel_quickstep_lib_test src/main_for_lib.cpp src/test_lib.cpp)
  target_link_libraries(parallel_quickstep_lib_test parallel_quickstep)

elseif( DEFINED USE_CPU )

  add_definitions(-DUSE_CPU)

  set(CPU_SOLVER_SOURCE_FILES
    src/parallel_stepper.cpp
    src/parallel_quickstep.cpp)

  rosbuild_add_library( parallel_quickstep
    ${CPU_SOLVER_SOURCE_FILES}
    )

  rosbuild_add_link_flags(parallel_quickstep)
  rosbuild_add_openmp_flags(parallel_quickstep)

endif()


