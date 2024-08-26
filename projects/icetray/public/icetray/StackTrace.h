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
#if defined(__clang__)
    void __cxa_throw(void *ex, std::type_info *info, void (_LIBCXXABI_DTOR_FUNC *dest)(void *));
#elif defined(__GNUC__) || defined(__GNUG__)
    void __cxa_throw(void *ex, void           *info, void (                     *dest)(void *));
#elif defined(_MSC_VER)
#endif
}

#endif // g4_larsim_StackTrace_H
