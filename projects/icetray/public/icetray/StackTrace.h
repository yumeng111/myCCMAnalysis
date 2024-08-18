#ifndef g4_larsim_StackTrace_H
#define g4_larsim_StackTrace_H

#include <map>
#include <vector>

#include <memory>
#include <string>
#include <cstdlib>
#include <dlfcn.h>
#include <cxxabi.h>
#include <iostream>
#include <typeinfo>
#include <execinfo.h>
#include <cpptrace/cpptrace.hpp>

extern "C" {
    void __cxa_throw(void *ex, void *info, void (*dest)(void *));
}

#endif // g4_larsim_StackTrace_H
