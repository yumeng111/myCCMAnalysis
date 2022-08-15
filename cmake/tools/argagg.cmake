
if(TARGET argagg)
else(TARGET argagg)
    add_subdirectory("${PROJECT_SOURCE_DIR}/extern/argagg" "extern/argagg")
    set(ARGAGG_FOUND ON)
    set(ARGAGG_LIBRARIES "")
    set(ARGAGG_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/extern/argagg/include")
endif(TARGET argagg)
