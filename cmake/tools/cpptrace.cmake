
colormsg(HICYAN "")
colormsg(HICYAN "cpptrace")

find_package(cpptrace QUIET)

tooldef(cpptrace
    include
    cpptrace/cpptrace.hpp
    lib64
    NONE
    cpptrace
)

