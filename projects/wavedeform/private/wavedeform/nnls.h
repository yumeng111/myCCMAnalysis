
#include <cholmod.h>

#include "timer.h"

/*
 * Lawson-Hanson NNLS
 *
 * Algorithm NNLS from "Solving Least Squares Problems", Charles Lawson and
 *  Richard Hanson. Prentice-Hall, 1974.
 */
cholmod_dense *
nnls_lawson_hanson(cholmod_sparse *A, cholmod_dense *y, double tolerance,
    unsigned min_iterations, unsigned max_iterations, unsigned npos,
    int normaleq, int solve_with_normaleq, int verbose, cholmod_common *c, DurationTimer & timer, size_t & iterations);


/* Solve Ap(x) = y by transforming to normal equations */
cholmod_dense *
solve_by_normaleq(cholmod_sparse *Ap, cholmod_dense *y,
    int nP, cholmod_common *c);
