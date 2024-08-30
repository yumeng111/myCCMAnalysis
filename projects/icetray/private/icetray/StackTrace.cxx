#include <icetray/StackTrace.h>

namespace {
    void * last_frames[20];
    size_t last_size;
    std::string exception_name;

    std::string demangle(const char *name) {
        int status;
        std::unique_ptr<char,void(*)(void*)> realname(abi::__cxa_demangle(name, 0, 0, &status), &std::free);
        return status ? "failed" : &*realname;
    }
}

extern "C" {
#if defined(__clang__)
    void __cxa_throw(void *ex, std::type_info *info, void (_LIBCXXABI_DTOR_FUNC *dest)(void *)) {
#elif defined(__GNUC__) || defined(__GNUG__)
    void __cxa_throw(void *ex, void           *info, void (                     *dest)(void *)) {
#elif defined(_MSC_VER)
#endif
        exception_name = demangle(reinterpret_cast<const std::type_info*>(info)->name());
        last_size = backtrace(last_frames, sizeof last_frames/sizeof(void*));
#if defined(__clang__)
        typedef void (*const ThingyConst)(void*,std::type_info*,void(*)(void*)) __attribute__ ((noreturn));
        typedef void (* ThingyNonConst)(void*,std::type_info*,void(*)(void*)) __attribute__ ((noreturn));
        static ThingyConst rethrow = (ThingyNonConst)dlsym(RTLD_NEXT, "__cxa_throw");
#elif defined(__GNUC__) || defined(__GNUG__)
        static void (*const rethrow)(void*,void*,void(*)(void*)) __attribute__ ((noreturn)) = (void (*)(void*,void*,void(*)(void*)))dlsym(RTLD_NEXT, "__cxa_throw");
#elif defined(_MSC_VER)
#endif
        rethrow(ex,info,dest);
    }
}
