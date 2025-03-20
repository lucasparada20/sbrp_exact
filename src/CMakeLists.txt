cmake_minimum_required(VERSION 3.2)
project(ExecExact)

# C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Hard code the path to the CPLEX directory
set(CPLEX_DIR "/home/lucasparada20/CPLEX_Studio2211")
message("Using CPLEX_DIR = ${CPLEX_DIR}")

# Include directories for CPLEX
set(CPLEX_INCLUDE_DIRS 
    ${CPLEX_DIR}/cplex/include
    ${CPLEX_DIR}/concert/include
)

# Link directories for CPLEX
set(CPLEX_LIBRARY_DIRS 
    ${CPLEX_DIR}/cplex/lib/x86-64_linux/static_pic
    ${CPLEX_DIR}/concert/lib/x86-64_linux/static_pic
)

# Flags for CPLEX
set(CPLEX_FLAGS ilocplex concert cplex m pthread dl)

message("CPLEX_INCLUDE_DIRS = ${CPLEX_INCLUDE_DIRS}")
message("CPLEX_LIBRARY_DIRS = ${CPLEX_LIBRARY_DIRS}")
message("CPLEX_FLAGS = ${CPLEX_FLAGS}")

# Create a shared interface target for CPLEX 
add_library(CPLEX INTERFACE)
target_include_directories(CPLEX INTERFACE ${CPLEX_INCLUDE_DIRS})
target_link_directories(CPLEX INTERFACE ${CPLEX_LIBRARY_DIRS})
target_link_libraries(CPLEX INTERFACE ${CPLEX_FLAGS})


# Check that all required paths and libraries exist
foreach(dir IN LISTS CPLEX_INCLUDE_DIRS CPLEX_LIBRARY_DIRS)
    if(NOT EXISTS ${dir})
        message(FATAL_ERROR "Required CPLEX path not found: ${dir}")
    endif()
endforeach()

# Source files
# All source files see the cplex headers
set(SOURCES
	main_ex_sbrp.cpp
	Parameters.cpp
	LoadSBRP.cpp
	UbManager.cpp
	CostFunctionSBRP.cpp
	SequentialInsertionSBRP.cpp
	InsRmvMethodSBRP.cpp
	ExactSbrp.cpp
	ExactSbrpMulticut.cpp
	ExactSbrpGraph.cpp
	Network.cpp
	ExactSbrpSep.cpp
	ExactSbrpSepMulticut.cpp
	ExactSbrpCallBacks.cpp
	ExactSbrpCallBackMulticut.cpp
	RouteFeasibility.cpp
	RecourseLowerBound.cpp
	DriverCountBppLb.cpp
	FirstStageDuals.cpp
)

# Target: exec_heur
add_executable(exec_exact ${SOURCES})
target_compile_options(exec_exact PRIVATE) 
target_link_libraries(exec_exact PRIVATE CPLEX)

# Usage message when typing make in the build directory
add_custom_target(usage
    COMMAND ${CMAKE_COMMAND} -E echo "Call: ./exec_exact to see how to use the executable"
)

