#ifndef MOTIF_OMP_H
#define MOTIF_OMP_H
/* Allow Apple Clang to build a serial CPU reference alongside the Metal path.
 * Builds with -fopenmp continue to use the real OpenMP runtime. */
#ifdef _OPENMP
#include <omp.h>
#else
#include <time.h>
typedef int omp_lock_t;
static inline int omp_get_max_threads(void) { return 1; }
static inline double omp_get_wtime(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}
static inline void omp_init_lock(omp_lock_t *lock) { *lock = 0; }
static inline void omp_destroy_lock(omp_lock_t *lock) { (void)lock; }
static inline void omp_set_lock(omp_lock_t *lock) { (void)lock; }
static inline void omp_unset_lock(omp_lock_t *lock) { (void)lock; }
#endif
#endif
