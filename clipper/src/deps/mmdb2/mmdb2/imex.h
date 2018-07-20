#ifdef _WIN32
# ifdef MMDB_EXPORT
#  define MMDB_IMEX __declspec(dllexport)
# else
#  define MMDB_IMEX __declspec(dllimport)
# endif
#else
# if (__GNUC__ > 4) || (__GNUC__ == 4 && (defined(__APPLE__) || __GNUC_MINOR__ >= 3))
#  define MMDB_IMEX __attribute__((visibility("default")))
# else
#  define MMDB_IMEX
# endif
#endif
