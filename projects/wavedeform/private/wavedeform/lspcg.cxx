/*
 * Copyright (c) 2018- IceCube Collaboration. All rights reserved.
 *
 * Authors:
 *    Jim Braun	 <jbraun@icecube.wisc.edu>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <cholmod.h>

/*
 *  Include <complex> here, before <SuiteSparseQR_C.h>, because older versions
 *  of SuiteSparse include <complex> within an extern "C" block
 */

#include <complex>
#include <SuiteSparseQR_C.h>

#include "lspcg.h"


#define EPSILON 0.00001

double
truncated_ratio(double s1, double s2) {
  if (s2 == 0.) {
    return 0.;
  }
  double r = s1 / s2;
  if (r < 0.) {
    r = 0.;
  }
  return r;
}

void
pcg_gradient(cholmod_sparse *A, cholmod_dense *y,
             cholmod_dense *s, cholmod_dense *g,
             cholmod_common *c) {

  /*
   *  Calculate the gradient vector g_i = sum_j[(y_j - s_j)A_ij]
   *  A: Basis matrix
   *  y: Data vector (len A->nrow)
   *  s: Projected solution vector (len A->nrow)
   *  g: Gradient vector to fill (len A->ncol)
   */
  cholmod_dense *diff = cholmod_l_zeros(A->nrow, 1, CHOLMOD_REAL, c);
  for (int i = 0; i < y->nrow; ++i) {
    ((double *)(diff->x))[i] = ((double *)(y->x))[i] - ((double *)(s->x))[i];
  }
  double alpha[2] = {1., 1.};
  double beta[2] = {0., 0.};
  cholmod_l_sdmult(A, 1 /* transpose */, alpha, beta, diff, g, c);
  cholmod_l_free_dense(&diff, c);
}


void
pcg_precondition(cholmod_sparse *A, cholmod_dense *x,
                 cholmod_dense *s, cholmod_dense *g,
                 cholmod_dense *v, cholmod_common *c) {
  /*
   *  Calculate the preconditioned gradient vector g'_i = C_ii * g_i
   *  C_ii = (x_i + EPSILON) / (sum_j[s_jA_ij] + EPSILON)
   *  A: Basis matrix
   *  x: Result vector (len A->ncol)
   *  s: Projected solution vector (len A->nrow)
   *  g: Gradient vector (len A->ncol)
   *  v: Preconditioned gradient vector to fill (len A->ncol)
   */
  // Step 1: Calculate the C_ii denominator:
  cholmod_dense *d = cholmod_l_zeros(A->ncol, 1, CHOLMOD_REAL, c);
  double alpha[2] = {1., 1.};
  double beta[2] = {0., 0.};
  cholmod_l_sdmult(A, 1 /* transpose */, alpha, beta, s, d, c);
  // Apply the preconditioning
  for (int i = 0; i < g->nrow; ++i) {
    double cii = (((double *)(x->x))[i] + EPSILON) /
                          (((double *)(d->x))[i] + EPSILON);
    ((double *)(v->x))[i] = ((double *)(g->x))[i] * cii;
  }
  cholmod_l_free_dense(&d, c);
}


double
pcg_stepsize(cholmod_dense *y, cholmod_dense *s, cholmod_dense *mu) {
  /*
   *  Calculate the pcg step size
   *  alpha = (sum_j[(y_j - s_j)mu_j] / sum_j[mu_j*mu_j]
   *  y: Data vector (len A->nrow)
   *  s: Projected solution vector (len A->nrow)
   *  mu: Projected search direction vector (len A->nrow)
   */
  double s1, s2;
  s1 = s2 = 0.;
  for (int i = 0; i < y->nrow; ++i) {
    s1 += (((double *)(y->x))[i] - ((double *)(s->x))[i]) *
                                              ((double *)(mu->x))[i];
    s2 += ((double *)(mu->x))[i] * ((double *)(mu->x))[i];
  }
  return truncated_ratio(s1, s2);
}


cholmod_dense*
lspcg(cholmod_sparse *A, cholmod_dense *y, double tolerance,
      unsigned max_iterations, cholmod_common *c) {

  /*
   *  Solve the non-negative least-squares problem y = Ax using PCG
   */

  // Step 1: Initialize solution vector
  cholmod_dense *x = cholmod_l_ones(A->ncol, 1, CHOLMOD_REAL, c);
  double alpha[2] = {1., 1.};
  double beta[2] = {0., 0.};
  cholmod_dense *temp = cholmod_l_zeros(A->nrow, 1, CHOLMOD_REAL, c);
  cholmod_l_sdmult(A, 0, alpha, beta, x, temp, c);
  double s1, s2;
  s1 = s2 = 0.;
  for (int i = 0; i < A->nrow; ++i) {
    s1 += ((double *)(y->x))[i];
    s2 += ((double *)(temp->x))[i];
  }
  if (s2 > 0) {
    s1 /= s2;
    for (int i = 0; i < A->ncol; ++i) {
      ((double *)(x->x))[i] *= s1;
    }
  }
  cholmod_l_free_dense(&temp, c);

  // Step 2: Initialize projected data vector
  cholmod_dense *s = cholmod_l_zeros(A->nrow, 1, CHOLMOD_REAL, c);
  cholmod_l_sdmult(A, 0, alpha, beta, x, s, c);

  // Steps 3 & 4: Initialize gradient and calculate
  // preconditioned gradient direction
  cholmod_dense *g = cholmod_l_zeros(A->ncol, 1, CHOLMOD_REAL, c);
  pcg_gradient(A, y, s, g, c);
  cholmod_dense *v = cholmod_l_zeros(A->ncol, 1, CHOLMOD_REAL, c);
  pcg_precondition(A, x, s, g, v, c);

  // Step 5: Initialize the search direction
  cholmod_dense *d = cholmod_l_copy_dense(v, c);

  // Steps 6-19: Iterate
  cholmod_dense *mu = cholmod_l_zeros(A->nrow, 1, CHOLMOD_REAL, c);
  cholmod_dense *muhat = cholmod_l_zeros(A->nrow, 1, CHOLMOD_REAL, c);
  cholmod_dense *dhat = cholmod_l_zeros(A->ncol, 1, CHOLMOD_REAL, c);
  cholmod_dense *xhat = cholmod_l_zeros(A->ncol, 1, CHOLMOD_REAL, c);
  cholmod_dense *vnew = cholmod_l_zeros(A->ncol, 1, CHOLMOD_REAL, c);
  cholmod_dense *gnew = cholmod_l_zeros(A->ncol, 1, CHOLMOD_REAL, c);
  for (unsigned k = 0; k < max_iterations; ++k) {

    // Step 6: Forward project the search direction
    cholmod_l_sdmult(A, 0, alpha, beta, d, mu, c);

    // Step 7: Calculate the step size:
    // In the least-squares case we avoid the line search, so this is fast.
    double a = pcg_stepsize(y, s, mu);

    // Step 8: Create the truncated intermediate solution vector
    for (int i = 0; i < A->ncol; ++i) {
      ((double *)(xhat->x))[i] =
                ((double *)(x->x))[i] + a*((double *)(d->x))[i];
      if (((double *)(xhat->x))[i] < 0.) {
        ((double *)(xhat->x))[i] = 0.;
      }
    }

    // Step 9: Create the intermediate search direction
    for (int i = 0; i < A->ncol; ++i) {
      ((double *)(dhat->x))[i] =
                ((double *)(xhat->x))[i] - ((double *)(x->x))[i];
    }

    // Step 10: Forward project the intermediate search direction
    cholmod_l_sdmult(A, 0, alpha, beta, dhat, muhat, c);

    // Step 11: Calculate the new step size
    a = pcg_stepsize(y, s, muhat);

    // Step 11a: Limit the step size to [0, 1]
    if (a > 1.) {
      a = 1.;
    }

    // Step 12: Update the solution vector
    for (int i = 0; i < A->ncol; ++i) {
      ((double *)(x->x))[i] += a * ((double *)(dhat->x))[i];
    }

    // Step 13: Update the projected data vector
    for (int i = 0; i < A->nrow; ++i) {
      ((double *)(s->x))[i] += a * ((double *)(muhat->x))[i];
    }

    // Steps 14 & 15: Update the gradient and preconditioned gradient vectors
    pcg_gradient(A, y, s, gnew, c);
    pcg_precondition(A, x, s, gnew, vnew, c);

    // Step 16: Calculate the Polak-Ribiere value
    s1 = s2 = 0.;
    for (int i = 0; i < A->ncol; ++i) {
      s1 += (((double *)(vnew->x))[i] - ((double *)(v->x))[i]) *
                                              ((double *)(gnew->x))[i];
      s2 += ((double *)(g->x))[i] * ((double *)(v->x))[i];
    }
    double b = truncated_ratio(s1, s2);

    // Step 16a: Copy the new gradient and search direction
    for (int i = 0; i < A->ncol; ++i) {
      ((double *)(v->x))[i] = ((double *)(vnew->x))[i];
      ((double *)(g->x))[i] = ((double *)(gnew->x))[i];
    }

    // Step 17: Update the search direction
    for (int i = 0; i < A->ncol; ++i) {
      ((double *)(d->x))[i] *= b;
      ((double *)(d->x))[i] += ((double *)(v->x))[i];
    }

    // Step 18: Ensure new direction is ascending
    s1 = 0.;
    for (int i = 0; i < A->ncol; ++i) {
      s1 += ((double *)(d->x))[i] * ((double *)(g->x))[i];
    }
    if (s1 <= 0.) {
      // Bad search direction.  Set search direction to gradient
      for (int i = 0; i < A->ncol; ++i) {
        ((double *)(d->x))[i] = ((double *)(g->x))[i];
      }
    }

    // Step 19: Check stopping criteria:
    double maxgrad = 0.;
    for (int i = 0; i < A->ncol; ++i) {
      if (fabs(((double *)(g->x))[i]) > maxgrad) {
        maxgrad = fabs(((double *)(g->x))[i]);
      }
    }
    if (maxgrad < tolerance) {
      break;
    }
  }

  // Done.  Clean up.
  cholmod_l_free_dense(&s, c);
  cholmod_l_free_dense(&g, c);
  cholmod_l_free_dense(&v, c);
  cholmod_l_free_dense(&d, c);
  cholmod_l_free_dense(&mu, c);
  cholmod_l_free_dense(&muhat, c);
  cholmod_l_free_dense(&dhat, c);
  cholmod_l_free_dense(&xhat, c);
  cholmod_l_free_dense(&vnew, c);
  cholmod_l_free_dense(&gnew, c);
  return x;
}
