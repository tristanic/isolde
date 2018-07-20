#pragma once
#ifdef _WIN32
# ifdef CLIPPER_EXPORTS
#  define CLIPPER_IMEX __declspec(dllexport)
# else
#  define CLIPPER_IMEX __declspec(dllimport)
# endif
#else
# if (__GNUC__ > 4) || (__GNUC__ == 4 && (defined(__APPLE__) || __GNUC_MINOR__ >= 3))
#  define CLIPPER_IMEX __attribute__((visibility("default")))
# else
#  define CLIPPER_IMEX
# endif
#endif
