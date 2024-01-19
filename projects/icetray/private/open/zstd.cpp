// (C) Copyright Reimar Döffinger 2018.
// Based on zstd.cpp by:
// (C) Copyright Milan Svoboda 2008.
// (C) Copyright Jonathan Turkanis 2003.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

// See http://www.boost.org/libs/iostreams for documentation.

// Define BOOST_IOSTREAMS_SOURCE so that <boost/iostreams/detail/config.hpp>
// knows that we are building the library (possibly exporting code), rather
// than using it (possibly importing code).
#define BOOST_IOSTREAMS_SOURCE

#include <iostream>
#include <zstd.h>

#include <boost/throw_exception.hpp>
#include <boost/iostreams/detail/config/dyn_link.hpp>
#include <boost/iostreams/filter/zstd.hpp>
#include "zstd.hpp"

namespace boost { namespace iostreams {

//------------------Implementation of zstd_multithread_base-------------------------------//

namespace detail {

zstd_multithread_base::zstd_multithread_base()
    : cstream_(ZSTD_createCStream()), dstream_(ZSTD_createDStream()), in_(new ZSTD_inBuffer), out_(new ZSTD_outBuffer), eof_(0)
    { }

zstd_multithread_base::~zstd_multithread_base()
{
    ZSTD_freeCStream(static_cast<ZSTD_CStream *>(cstream_));
    ZSTD_freeDStream(static_cast<ZSTD_DStream *>(dstream_));
    delete static_cast<ZSTD_inBuffer*>(in_);
    delete static_cast<ZSTD_outBuffer*>(out_);
}

void zstd_multithread_base::before( const char*& src_begin, const char* src_end,
                        char*& dest_begin, char* dest_end )
{
    ZSTD_inBuffer *in = static_cast<ZSTD_inBuffer *>(in_);
    ZSTD_outBuffer *out = static_cast<ZSTD_outBuffer *>(out_);
    in->src = src_begin;
    in->size = static_cast<size_t>(src_end - src_begin);
    in->pos = 0;
    out->dst = dest_begin;
    out->size = static_cast<size_t>(dest_end - dest_begin);
    out->pos = 0;
}

void zstd_multithread_base::after(const char*& src_begin, char*& dest_begin, bool)
{
    ZSTD_inBuffer *in = static_cast<ZSTD_inBuffer *>(in_);
    ZSTD_outBuffer *out = static_cast<ZSTD_outBuffer *>(out_);
    src_begin = reinterpret_cast<const char*>(in->src) + in->pos;
    dest_begin = reinterpret_cast<char*>(out->dst) + out->pos;
}

int zstd_multithread_base::deflate(int action)
{
    ZSTD_CStream *s = static_cast<ZSTD_CStream *>(cstream_);
    ZSTD_inBuffer *in = static_cast<ZSTD_inBuffer *>(in_);
    ZSTD_outBuffer *out = static_cast<ZSTD_outBuffer *>(out_);
    // Ignore spurious extra calls.
    // Note size > 0 will trigger an error in this case.
    if (eof_ && in->size == 0) return zstd::stream_end;
    ZSTD_EndDirective directive = action == zstd::finish ? ZSTD_e_end : (action == zstd::flush ? ZSTD_e_flush : ZSTD_e_continue);
    size_t result = ZSTD_compressStream2(s, out, in, directive);
    zstd_error::check BOOST_PREVENT_MACRO_SUBSTITUTION(result);
    if (action != zstd::run)
    {
        while((action == zstd::finish && result != 0) or (action == zstd::flush && in->pos < in->size)) {
            result = ZSTD_compressStream2(s, out, in, directive);
            zstd_error::check BOOST_PREVENT_MACRO_SUBSTITUTION(result);
        }
        eof_ = action == zstd::finish && result == 0;
        return result == 0 ? zstd::stream_end : zstd::okay;
    }
    return zstd::okay;
}

int zstd_multithread_base::inflate(int action)
{
    ZSTD_DStream *s = static_cast<ZSTD_DStream *>(dstream_);
    ZSTD_inBuffer *in = static_cast<ZSTD_inBuffer *>(in_);
    ZSTD_outBuffer *out = static_cast<ZSTD_outBuffer *>(out_);
    // need loop since iostream code cannot handle short reads
    do {
        size_t result = ZSTD_decompressStream(s, out, in);
        zstd_error::check BOOST_PREVENT_MACRO_SUBSTITUTION(result);
    } while (in->pos < in->size && out->pos < out->size);
    return action == zstd::finish && in->size == 0 && out->pos == 0 ? zstd::stream_end : zstd::okay;
}

void zstd_multithread_base::reset(bool compress, bool realloc)
{
    ZSTD_inBuffer *in = static_cast<ZSTD_inBuffer *>(in_);
    ZSTD_outBuffer *out = static_cast<ZSTD_outBuffer *>(out_);
    if (realloc)
    {
        memset(in, 0, sizeof(*in));
        memset(out, 0, sizeof(*out));
        eof_ = 0;

        zstd_error::check BOOST_PREVENT_MACRO_SUBSTITUTION(
            compress ?
                ZSTD_initCStream(static_cast<ZSTD_CStream *>(cstream_), level) :
                ZSTD_initDStream(static_cast<ZSTD_DStream *>(dstream_))
        );
        if(compress)
            ZSTD_CCtx_setParameter(static_cast<ZSTD_CStream *>(cstream_), ZSTD_c_nbWorkers, 4);
    }
}

void zstd_multithread_base::do_init
    ( const zstd_params& p, bool compress, unsigned int n_workers,
      zstd::alloc_func, zstd::free_func,
      void* )
{
    ZSTD_inBuffer *in = static_cast<ZSTD_inBuffer *>(in_);
    ZSTD_outBuffer *out = static_cast<ZSTD_outBuffer *>(out_);

    memset(in, 0, sizeof(*in));
    memset(out, 0, sizeof(*out));
    eof_ = 0;

    level = p.level;
    zstd_error::check BOOST_PREVENT_MACRO_SUBSTITUTION(
        compress ?
            ZSTD_initCStream(static_cast<ZSTD_CStream *>(cstream_), level) :
            ZSTD_initDStream(static_cast<ZSTD_DStream *>(dstream_))
    );
    if(compress)
        ZSTD_CCtx_setParameter(static_cast<ZSTD_CStream *>(cstream_), ZSTD_c_nbWorkers, n_workers);
}

} // End namespace detail.

//----------------------------------------------------------------------------//

} } // End namespaces iostreams, boost.
