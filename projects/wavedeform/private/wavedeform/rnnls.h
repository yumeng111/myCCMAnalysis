#include <cholmod.h>
#include "timer.h"

cholmod_dense *
rnnls(cholmod_sparse *A, cholmod_dense *y, double tolerance,
      unsigned max_iterations, int verbose, cholmod_common *c, DurationTimer & timer);
