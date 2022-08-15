# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# nlohmann/json
set(JSON_Install ON CACHE INTERNAL "")
set(JSON_BuildTests OFF CACHE INTERNAL "")
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/json-3.11.0" "extern/json-3.11.0")
include_directories("${PROJECT_SOURCE_DIR}/extern/json-3.11.0/include")

# argagg
#add_subdirectory("${PROJECT_SOURCE_DIR}/extern/argagg" "extern/argagg")
#include_directories("${PROJECT_SOURCE_DIR}/extern/argagg/include")
