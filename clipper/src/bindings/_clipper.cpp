#include "../molc.h"

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-cif.h>
#include <clipper/clipper-cns.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-phs.h>

using namespace clipper;

#define DEFAULT_NEW(FNAME, CNAME) \
extern "C" EXPORT void* \
FNAME##_new() \
{ \
    try { \
        auto ptr = new CNAME(); \
        return ptr; \
    } catch (...) { \
        molc_error(); \
        return 0; \
    } \
}

#define DEFAULT_DEL(FNAME, CNAME) \
extern "C" EXPORT void \
FNAME##_delete(void* ptr) \
{ \
    try { \
        auto p = static_cast<CNAME *>(ptr); \
        delete p; \
    } catch (...) { \
        molc_error(); \
    } \
}

DEFAULT_NEW(mtzfile, CCP4MTZfile)
DEFAULT_DEL(mtzfile, CCP4MTZfile)

extern "C" EXPORT void
mtzfile_open_read(void *handler, const char *fname)
{
    auto h = static_cast<CCP4MTZfile *>(handler);
    error_wrap(h, &CCP4MTZfile::open_read, String(fname));
}

extern "C" EXPORT void
mtzfile_close_read(void *handler)
{
    auto h = static_cast<CCP4MTZfile *>(handler);
    error_wrap(h, &CCP4MTZfile::close_read);
}

extern "C" EXPORT void
mtzfile_open_write(void *handler, const char *fname)
{
    auto h = static_cast<CCP4MTZfile *>(handler);
    error_wrap(h, &CCP4MTZfile::open_write, String(fname));
}

extern "C" EXPORT void
mtzfile_close_write(void *handler)
{
    auto h = static_cast<CCP4MTZfile *>(handler);
    error_wrap(h, &CCP4MTZfile::close_write);
}
