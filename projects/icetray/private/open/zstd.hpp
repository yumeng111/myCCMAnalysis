// (C) Copyright Reimar Döffinger 2018.
// Based on zstd.hpp by:
// (C) Copyright Milan Svoboda 2008.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

// See http://www.boost.org/libs/iostreams for documentation.

#ifndef I3_ZSTD_BOOST_IMPL
#define I3_ZSTD_BOOST_IMPL

#if defined(_MSC_VER)
# pragma once
#endif

#include <cassert>
#include <iosfwd>            // streamsize.
#include <memory>            // allocator, bad_alloc.
#include <new>
#include <boost/config.hpp>  // MSVC, STATIC_CONSTANT, DEDUCED_TYPENAME, DINKUM.
#include <boost/detail/workaround.hpp>
#include <boost/iostreams/constants.hpp>   // buffer size.
#include <boost/iostreams/detail/config/auto_link.hpp>
#include <boost/iostreams/detail/config/dyn_link.hpp>
#include <boost/iostreams/detail/config/wide_streams.hpp>
#include <boost/iostreams/detail/ios.hpp>  // failure, streamsize.
#include <boost/iostreams/filter/symmetric.hpp>
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/iostreams/pipeline.hpp>
#include <boost/type_traits/is_same.hpp>

namespace boost { namespace iostreams {

namespace detail {

class BOOST_IOSTREAMS_DECL zstd_multithread_base {
public:
    typedef char char_type;
protected:
    zstd_multithread_base();
    ~zstd_multithread_base();
    template<typename Alloc>
    void init( const zstd_params& p,
               bool compress,
               unsigned int n_workers,
               zstd_allocator<Alloc>& zalloc )
        {
            bool custom = zstd_allocator<Alloc>::custom;
            do_init( p, compress, n_workers,
                     custom ? zstd_allocator<Alloc>::allocate : 0,
                     custom ? zstd_allocator<Alloc>::deallocate : 0,
                     &zalloc );
        }
    void before( const char*& src_begin, const char* src_end,
                 char*& dest_begin, char* dest_end );
    void after( const char*& src_begin, char*& dest_begin,
                bool compress );
    int deflate(int action);
    int inflate(int action);
    void reset(bool compress, bool realloc);
private:
    void do_init( const zstd_params& p, bool compress,
                  unsigned int n_workers,
                  zstd::alloc_func,
                  zstd::free_func,
                  void* derived );
    void*         cstream_;         // Actual type: ZSTD_CStream *
    void*         dstream_;         // Actual type: ZSTD_DStream *
    void*         in_;              // Actual type: ZSTD_inBuffer *
    void*         out_;             // Actual type: ZSTD_outBuffer *
    int eof_;
    uint32_t level;
    unsigned int n_workers;
};

//
// Template name: zstd_compressor_multithread_impl
// Description: Model of C-Style Filter implementing compression by
//      delegating to the zstd function deflate.
//
template<typename Alloc = std::allocator<char> >
class zstd_compressor_multithread_impl : public zstd_multithread_base, public zstd_allocator<Alloc> {
public:
    zstd_compressor_multithread_impl(const zstd_params& = zstd::default_compression, unsigned int n_workers = 1);
    ~zstd_compressor_multithread_impl();
    bool filter( const char*& src_begin, const char* src_end,
                 char*& dest_begin, char* dest_end, bool flush );
    void close();
};

//
// Template name: zstd_compressor_multithread_impl
// Description: Model of C-Style Filte implementing decompression by
//      delegating to the zstd function inflate.
//
template<typename Alloc = std::allocator<char> >
class zstd_decompressor_multithread_impl : public zstd_multithread_base, public zstd_allocator<Alloc> {
public:
    zstd_decompressor_multithread_impl(const zstd_params&);
    zstd_decompressor_multithread_impl();
    ~zstd_decompressor_multithread_impl();
    bool filter( const char*& begin_in, const char* end_in,
                 char*& begin_out, char* end_out, bool flush );
    void close();
};

} // End namespace detail.

//
// Template name: zstd_compressor
// Description: Model of InputFilter and OutputFilter implementing
//      compression using zstd.
//
template<typename Alloc = std::allocator<char> >
struct multithread_zstd_compressor
    : symmetric_filter<detail::zstd_compressor_multithread_impl<Alloc>, Alloc>
{
private:
    typedef detail::zstd_compressor_multithread_impl<Alloc> impl_type;
    typedef symmetric_filter<impl_type, Alloc>  base_type;
public:
    typedef typename base_type::char_type               char_type;
    typedef typename base_type::category                category;
    multithread_zstd_compressor( const zstd_params& = zstd::default_compression, unsigned int n_workers = 1,
                           std::streamsize buffer_size = default_device_buffer_size );
};
BOOST_IOSTREAMS_PIPABLE(multithread_zstd_compressor, 1)

//
// Template name: zstd_decompressor
// Description: Model of InputFilter and OutputFilter implementing
//      decompression using zstd.
//
template<typename Alloc = std::allocator<char> >
struct multithread_zstd_decompressor
    : symmetric_filter<detail::zstd_decompressor_multithread_impl<Alloc>, Alloc>
{
private:
    typedef detail::zstd_decompressor_multithread_impl<Alloc> impl_type;
    typedef symmetric_filter<impl_type, Alloc>    base_type;
public:
    typedef typename base_type::char_type               char_type;
    typedef typename base_type::category                category;
    multithread_zstd_decompressor( std::streamsize buffer_size = default_device_buffer_size );
    multithread_zstd_decompressor( const zstd_params& p,
                             std::streamsize buffer_size = default_device_buffer_size );
};
BOOST_IOSTREAMS_PIPABLE(multithread_zstd_decompressor, 1)


//----------------------------------------------------------------------------//

//------------------Implementation of zstd_allocator--------------------------//

namespace detail {

//------------------Implementation of zstd_compressor_multithread_impl--------------------//

template<typename Alloc>
zstd_compressor_multithread_impl<Alloc>::zstd_compressor_multithread_impl(const zstd_params& p, unsigned int n_workers)
{ init(p, true, n_workers, static_cast<zstd_allocator<Alloc>&>(*this)); }

template<typename Alloc>
zstd_compressor_multithread_impl<Alloc>::~zstd_compressor_multithread_impl()
{ reset(true, false); }

template<typename Alloc>
bool zstd_compressor_multithread_impl<Alloc>::filter
    ( const char*& src_begin, const char* src_end,
      char*& dest_begin, char* dest_end, bool flush )
{
    before(src_begin, src_end, dest_begin, dest_end);
    int result = deflate(flush ? zstd::flush : zstd::run);
    after(src_begin, dest_begin, true);
    return result != zstd::stream_end;
}

template<typename Alloc>
void zstd_compressor_multithread_impl<Alloc>::close() { reset(true, true); }

//------------------Implementation of zstd_decompressor_multithread_impl------------------//

template<typename Alloc>
zstd_decompressor_multithread_impl<Alloc>::zstd_decompressor_multithread_impl(const zstd_params& p)
{ init(p, false, 1, static_cast<zstd_allocator<Alloc>&>(*this)); }

template<typename Alloc>
zstd_decompressor_multithread_impl<Alloc>::~zstd_decompressor_multithread_impl()
{ reset(false, false); }

template<typename Alloc>
zstd_decompressor_multithread_impl<Alloc>::zstd_decompressor_multithread_impl()
{
    zstd_params p;
    init(p, false, static_cast<zstd_allocator<Alloc>&>(*this));
}

template<typename Alloc>
bool zstd_decompressor_multithread_impl<Alloc>::filter
    ( const char*& src_begin, const char* src_end,
      char*& dest_begin, char* dest_end, bool flush )
{
    before(src_begin, src_end, dest_begin, dest_end);
    int result = inflate(flush ? zstd::finish : zstd::run);
    after(src_begin, dest_begin, false);
    return result != zstd::stream_end;
}

template<typename Alloc>
void zstd_decompressor_multithread_impl<Alloc>::close() { reset(false, true); }

} // End namespace detail.

//------------------Implementation of zstd_compressor-----------------------//

template<typename Alloc>
multithread_zstd_compressor<Alloc>::multithread_zstd_compressor
    (const zstd_params& p, unsigned int n_workers, std::streamsize buffer_size)
    : base_type(buffer_size, p, n_workers) { }

//------------------Implementation of zstd_decompressor-----------------------//

template<typename Alloc>
multithread_zstd_decompressor<Alloc>::multithread_zstd_decompressor
    (std::streamsize buffer_size)
    : base_type(buffer_size) { }

template<typename Alloc>
multithread_zstd_decompressor<Alloc>::multithread_zstd_decompressor
    (const zstd_params& p, std::streamsize buffer_size)
    : base_type(buffer_size, p) { }

//----------------------------------------------------------------------------//

} } // End namespaces iostreams, boost.

#endif // #ifndef I3_ZSTD_BOOST_IMPL
