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
    void __cxa_throw(void *ex, void *info, void (*dest)(void *)) {
        exception_name = demangle(reinterpret_cast<const std::type_info*>(info)->name());
        last_size = backtrace(last_frames, sizeof last_frames/sizeof(void*));

        cpptrace::generate_trace().print();

        static void (*const rethrow)(void*,void*,void(*)(void*)) __attribute__ ((noreturn)) = (void (*)(void*,void*,void(*)(void*)))dlsym(RTLD_NEXT, "__cxa_throw");
        rethrow(ex,info,dest);
    }
}
