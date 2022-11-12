
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>

#include <cholmod.h>

/*
 *  Include <complex> here, before <SuiteSparseQR_C.h>, because older versions
 *  of SuiteSparse include <complex> within an extern "C" block
 */

#include <complex>
#include <vector>
#include <set>
#include <algorithm>
#include <SuiteSparseQR_C.h>

#include "rnnls.h"

double SD_ONE[2] = {1., 0.};
double SD_ZERO[2] = {0., 0.};
double SD_MINUS_ONE[2] = {-1., 0.};


typedef struct {

  /* Allocate these objects */
  cholmod_dense *x, *Aty, *rowwork, *colwork;
  cholmod_factor *L;
  cholmod_sparse *update;
  long *P, *Z;

  /* Simple data */
  long nP, nZ;
  long last_free;

  /* Hold on to useful external pointers */
  cholmod_sparse *A;
  cholmod_dense *y;
  cholmod_common *c;

} rnnls_context;


rnnls_context*
rnnls_allocate_context(cholmod_sparse *A,
                       cholmod_dense *y,
                       cholmod_common *c) {

  rnnls_context *cxt = (rnnls_context*)malloc(sizeof(rnnls_context));
  cxt->x = cholmod_l_zeros(A->ncol, 1, CHOLMOD_REAL, c);

  /* Precompute Aty: This is the RHS if the system we're solving */
  cxt->Aty = cholmod_l_zeros(A->ncol, 1, CHOLMOD_REAL, c);
  cholmod_l_sdmult(A, 1, SD_ONE, SD_ZERO, y, cxt->Aty, c);

  cxt->rowwork = cholmod_l_zeros(A->ncol, 1, CHOLMOD_REAL, c);
  cxt->colwork = cholmod_l_zeros(A->nrow, 1, CHOLMOD_REAL, c);

  cxt->update = cholmod_l_allocate_sparse(A->ncol, 1, A->ncol,
        1, 1, 0, CHOLMOD_REAL, c);

  cxt->L = cholmod_l_allocate_factor(A->ncol, c);

  cxt->P = static_cast<long*>(malloc(A->ncol * sizeof(long)));
  cxt->Z = static_cast<long*>(malloc(A->ncol * sizeof(long)));
  cxt->A = A;
  cxt->y = y;
  cxt->c = c;

  /* Start with all elements actively constrained */
  cxt->nZ = A->ncol;
  cxt->nP = 0;
  for (long i = 0; i < cxt->nZ; i++) {
    cxt->Z[i] = i;
  }

  cxt->last_free = -1;
  return cxt;
}


void
rnnls_free_context(rnnls_context *cxt) {

  free(cxt->Z);
  free(cxt->P);
  cholmod_l_free_factor(&(cxt->L), cxt->c);
  cholmod_l_free_sparse(&(cxt->update), cxt->c);
  cholmod_l_free_dense(&(cxt->colwork), cxt->c);
  cholmod_l_free_dense(&(cxt->rowwork), cxt->c);
  cholmod_l_free_dense(&(cxt->Aty), cxt->c);
  cholmod_l_free_dense(&(cxt->x), cxt->c);
  free(cxt);
}


long
rnnls_next_free(double tolerance,
                rnnls_context *cxt) {

  /* Calculate negative gradient */
  cholmod_l_copy_dense2(cxt->y, cxt->colwork, cxt->c);
  // y - Ax -> cxt->colwork
  cholmod_l_sdmult(cxt->A, 0, SD_MINUS_ONE,
                   SD_ONE, cxt->x, cxt->colwork, cxt->c);
  // At(y - Ax) -> cxt->rowwork
  cholmod_l_sdmult(cxt->A, 1, SD_ONE, SD_ZERO,
                   cxt->colwork, cxt->rowwork, cxt->c);

  /* Find the index of the member in passive set with largest
   * negative gradient
   */

  // If we have no passive members, we can't free any member
  if (cxt->nZ == 0) {
    return -1;
  }

  double *wx = ((double *)((cxt->rowwork)->x));
  long *Z = cxt->Z;
  long *P = cxt->P;

  // Find the passive member with the largest gradient
  double wmax = wx[Z[0]];
  long zidx = 0;
  for (long i = 1; i < cxt->nZ; i++) {
    if (wx[Z[i]] > wmax && cxt->last_free != Z[i]) {
      zidx = i;
      wmax = wx[Z[zidx]];
    }
  }

  // If no members have positive gradient, we're done
  if (wmax <= 0) {
    return -1;
  }

  if (wmax < tolerance) {
    /* Check if any passive coefficients need to be reduced
     * due to clipped ringing */
    if (cxt->nP == 0) {
      return -1;
    }

    double wpmin = wx[P[0]];
    for (long i = 1; i < cxt->nP; i++) {
      if (wx[P[i]] < wpmin) {
        wpmin = wx[P[i]];
      }
    }

    if (-wpmin < tolerance) {
      // All passive members are within tolerance
      return -1;
    }
  }

  // Member Z[zidx] should be freed next
  return zidx;
}

void
rnnls_set_passive(long zidx,
                  rnnls_context *cxt) {

  cholmod_sparse *A = cxt->A;
  cholmod_sparse *update = cxt->update;
  long *P = cxt->P;
  long *Z = cxt->Z;

  /* Move coefficient Z[zidx] into the passive set P */
  cxt->last_free = Z[zidx];
  P[(cxt->nP)++] = Z[zidx];
  (cxt->nZ)--;
  for (long i = zidx; i < cxt->nZ; ++i) {
    Z[i] = Z[i+1];
  }

  /* Calculate the factorization update */
  /* First, extract column k from A.  Avoid a costly submatrix call */
  // A_k --> cxt->colwork
  memset((cxt->colwork)->x, '\0', sizeof(double) * (cxt->colwork)->nrow);
  long colstart = ((long *)(A->p))[cxt->last_free];
  long colend = ((long *)(A->p))[cxt->last_free + 1];
  if (!(A->packed)) {
    colend = colstart + ((long *)(A->nz))[cxt->last_free];
  }
  for (long idx = colstart; idx < colend; ++idx) {
    ((double*)((cxt->colwork)->x))[((long *)(A->i))[idx]] =
                                          ((double *)(A->x))[idx];
  }

  /* Now calculate AtA for column k */
  // AtA_k --> cxt->rowwork
  cholmod_l_sdmult(A, 1, SD_ONE, SD_ZERO, cxt->colwork, cxt->rowwork, cxt->c);

  /* Insert the entries of AtA for the passive set into the update matrix
   * NB: Reuse update matrix since allocation of cholmod_sparse is slow
   */
  ((long*)(update->p))[0] = 0;
  ((long*)(update->p))[1] = cxt->nP;
  for (long i = 0; i < cxt->nP; ++i) {
    ((long*)(update->i))[i] = P[i];
    ((double*)(update->x))[i] = ((double*)((cxt->rowwork)->x))[P[i]];
  }

  /* Update the factorization */
  cholmod_l_rowadd(cxt->last_free, cxt->update, cxt->L, cxt->c);
}

void
rnnls_set_active(long pidx,
                 rnnls_context *cxt) {

  long *P = cxt->P;
  long *Z = cxt->Z;

  /* Move coefficient P[pidx] into the active set Z */
  long kcol = P[pidx];
  Z[(cxt->nZ)++] = P[pidx];
  (cxt->nP)--;
  for (long i = pidx; i < cxt->nP; ++i) {
    P[i] = P[i+1];
  }
  cholmod_l_rowdel(kcol, NULL, cxt->L, cxt->c);
}


int
rnnls_solve(rnnls_context *cxt) {

  cholmod_dense *x = cxt->x;
  double *xx = (double *)(x->x);
  long *P = cxt->P;

  for (;;) {

    /* Step 6b: Update solution */
    cholmod_dense *p = cholmod_l_solve(CHOLMOD_A, cxt->L, cxt->Aty, cxt->c);
    double *px = (double *)(p->x);

    /* Step 6c: Sanity check solution */
    for (long i = 0; i < cxt->nP; ++i) {
      if (std::isnan(px[P[i]]) || std::isinf(px[P[i]])) {

        // This should not happen with reasonable optimal_tolerance.
        // Solve (expensively) by QR factorization of the subproblem.
        // Note that c->status does not alert us of potential problems.
        cholmod_sparse *Ap = cholmod_l_submatrix(cxt->A, NULL, -1, P,
                                                 cxt->nP, 1, 1, cxt->c);
        cholmod_dense *pp =
            SuiteSparseQR_C_backslash_default(Ap, cxt->y, cxt->c);
        for (int i = 0; i < cxt->nP; ++i) {
          px[P[i]] = ((double *)(pp->x))[i];
        }
        cholmod_l_free_sparse(&Ap, cxt->c);
        cholmod_l_free_dense(&pp, cxt->c);
      }
    }

    /*
     * Step 7: Check if any coefficients need be constrained
     */
    long ridx;
    for (ridx = 0; ridx < cxt->nP; ++ridx) {
      if (px[P[ridx]] <= 0) {
        break;
      }
    }

    if (ridx == cxt->nP) {
      /*
       * All were positive. Cycle back for the next.
       */
      memset(xx, '\0', sizeof(double) * x->nrow);
      for (long i = 0; i < cxt->nP; ++i) {
        xx[P[i]] = px[P[i]];
      }
      cholmod_l_free_dense(&p, cxt->c);
      return 0;
    }

    /*
     * Step 8-9: Compute q, alpha
     */
    double alpha = 2; /* All computed values <= 1 */
    long qmax = -1;
    for (long i = 0; i < cxt->nP; i++) {
      if (px[P[i]] > 0) {
        continue;
      }

      double qtemp = xx[P[i]] / (xx[P[i]] - px[P[i]]);

      if (qtemp < alpha && qtemp != 0) {
        qmax = P[i];
        alpha = qtemp;
      } else if (cxt->last_free == P[i]) {
        /* Anti-cycling advice from LH */
        alpha = 0;
        qmax = P[i];
        break;
      }
    }

    if (qmax < 0) {
      fprintf(stderr, "%s line %d: Math has failed\n", __FILE__, __LINE__);
      exit(1);
    }

    /*
     * Step 10: Recompute x
     */
    for (long i = 0; i < cxt->nP; i++) {
      xx[P[i]] += alpha * (px[P[i]] - xx[P[i]]);
    }

    /* Avoid rounding errors above */
    xx[qmax] = 0;
    cholmod_l_free_dense(&p, cxt->c);

    /*
     * Step 11: Move coefficients equal to zero to the
     * active set.
     */
    for (long i = 0; i < cxt->nP; ++i) {
      if (xx[P[i]] > 0) {
        continue;
      }

      rnnls_set_active(i, cxt);
      i--;
    }

    /* If alpha = 0, we've reached equilibrium */
    if (alpha == 0) {
      break;
    }
  }

  /* Exit to the caller in equilibrium */
  return -1;
}


/*
 *  Sort members of P by amplitude, ascending.
 */
void
rnnls_sort_P_by_amplitude(rnnls_context *cxt) {
  std::vector<std::pair<double, long> > pNew;
  double *xx = (double *)((cxt->x)->x);
  long *P = cxt->P;
  for (long i = 0; i < cxt->nP; ++i) {
    pNew.push_back(std::pair<double, long>(xx[P[i]], P[i]));
  }
  std::sort(pNew.begin(), pNew.end());
  for (long i = 0; i < cxt->nP; ++i) {
    P[i] = (pNew[i]).second;
  }
}

/*
 *  Find the index in Z of member idx
 */
int
rnnls_find_in_Z(long idx,
                rnnls_context *cxt) {

  for (long i = cxt->nZ - 1; i >= 0; --i) {
    if ((cxt->Z)[i] == idx) {
      return i;
    }
  }

  return -1;
}


/*
 *  Perform one forward iteration of NNLS
 */
int
rnnls_iterate(double tolerance,
              rnnls_context *cxt) {

  // Find the member of the active set to free next
  long next_free = rnnls_next_free(tolerance, cxt);
  if (next_free == -1) {
    // We've converged.
    return -1;
  }

  // Free member at index Z[next_free] and re-solve
  rnnls_set_passive(next_free, cxt);
  if (rnnls_solve(cxt) < 0) {
    // Solution is cycling or cannot be improved.
    return -1;
  }

  return 0;
}


/*
 *  Note: Here we avoid explicitly calculating AtA because of CPU cost.
 *  Calculate only the entries of AtA we need.  Since we solve in the
 *  normal equations, the problem must be sufficiently well-conditioned.
 */
cholmod_dense*
rnnls(cholmod_sparse *A, cholmod_dense *y, double tolerance,
      unsigned max_iterations, int verbose, cholmod_common *c) {

  rnnls_context *cxt = rnnls_allocate_context(A, y, c);

  // Find the optimal NNLS solution
  for (unsigned n = 0; n < max_iterations || max_iterations == 0; ++n) {

    // Machine precision limits our ability find the optimal solution.
    // Use a stopping tolerance suggested by Adlers
    double tol = 1.e2 * A->ncol * DBL_EPSILON;
    if (n == 0) {
      // If we're below tolerance at the first iteration, we can just quit
      tol = tolerance;
    }

    // Do an iteration of NNLS
    if (rnnls_iterate(tol, cxt) == -1) {
      break;
    }
  }

  // We've found the optimal solution.  Now remove passive members in order
  // of increasing amplitude.  Put them back if we cannot stay under tolerance
  // NB: we need to reset last_free since there is now no danger of cycling
  //     and any passive member can reasonably be removed.
  cxt->last_free = -1;
  cholmod_dense *ret = cholmod_l_copy_dense(cxt->x, c);

  // Track members we haven't checked.  Positions of members in P may change.
  std::set<long> Pset(cxt->P, cxt->P + cxt->nP);
  rnnls_sort_P_by_amplitude(cxt);
  for (;;) {

    // Find next member to test.
    long next_active = -1;
    for (long i = 0; i < cxt->nP; ++i) {
      if (Pset.count((cxt->P)[i]) > 0) {
        next_active = i;
        Pset.erase((cxt->P)[i]);
        break;
      }
    }

    // We're done when there is no member left to test.  Also, we already
    // know that if we have one passive member and we remove it, we'll be
    // above tolerance.
    if (next_active == -1 || cxt->nP == 1) {
      break;
    }

    // Remove the member at P[next_active] and solve the system again
    long member_index = (cxt->P)[next_active];
    rnnls_set_active(next_active, cxt);
    if (rnnls_solve(cxt) < 0) {
      // This should never happen as we're only removing passive members
      fprintf(stderr, "%s line %d: Solve failed\n", __FILE__, __LINE__);
      exit(1);
    }

    // Check tolerance
    if (rnnls_next_free(tolerance, cxt) == -1) {

      // We're below tolerance.  Accept the current solution and re-sort P
      cholmod_l_copy_dense2(cxt->x, ret, c);
      rnnls_sort_P_by_amplitude(cxt);

    } else {

      // We're above tolerance, so put this member back. No need to re-sort P
      rnnls_set_passive(rnnls_find_in_Z(member_index, cxt), cxt);
      cxt->last_free = -1;
      cholmod_l_copy_dense2(ret, cxt->x, c);
    }
  }

  rnnls_free_context(cxt);
  return ret;
}
