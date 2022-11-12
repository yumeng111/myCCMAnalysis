
#include <cholmod.h>

cholmod_dense*
lspcg(cholmod_sparse *A, cholmod_dense *y, double tolerance,
      unsigned max_iterations, cholmod_common *c);
