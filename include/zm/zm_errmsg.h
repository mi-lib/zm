/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_errmsg.h
 * \brief error messages.
 * \author Zhidao
 */

#ifndef __ZM_ERRMSG_H__
#define __ZM_ERRMSG_H__

/* NOTE: never include this header file in user programs. */

/* warning message */

#define ZM_WARN_BESSEL_UNDEF "undefined for negative integer"
#define ZM_WARN_SIZINV_VEC   "invalid vector size, ignored"
#define ZM_WARN_SIZINV_MAT   "invalid matrix size, ignored"

#define ZM_WARN_VEC_SIZMIS   "vector has only %d components, while specified size is %d"
#define ZM_WARN_MAT_SIZMIS   "matrix has only %d components, while specified size is %d"

#define ZM_WARN_VECLIST_EMPTY "empty vector list assigned"

#define ZM_WARN_LE_ZEROPIVOT "cannot sweep out by zero pivot"

#define ZM_WARN_NLE_NOROOT   "probably, equation has no root"

#define ZM_WARN_SEQ_EMPTY    "empty vector sequence assigned"
#define ZM_WARN_SEQ_STEP     "too large step index"

#define ZM_WARN_PEX_SIZMIS   "polynomial expression has only %d components, while specified size is %d"

#define ZM_WARN_INSUFFICIENT_SAMPLES "maybe insufficient number of samples %d for %d parameters"

#define ZM_WARN_OPT_BADINI   "probably inadequate initial guess"

#define ZM_WARN_ODE_GEAR1    "invalid step number %d, modified to 1"
#define ZM_WARN_ODE_GEAR2    "step number %d over 6 may cause instability, modified"

#define ZM_WARN_GRAPH_DUPSPC "duplicate specification of a node"
#define ZM_WARN_GRAPH_DUPCON "duplicate node connection, overwritten"

#define ZM_MCA_NOSILHOUETTE  "silhouette score not computed, need to call zMClusterSilhouetteScore()"

/* error messages */

#define ZM_ERR_ZERODIV       "cannot devide by tiny value"
#define ZM_ERR_OUTOFRANGE    "specified index out of range"

#define ZM_ERR_NULLVEC       "null vector assigned"
#define ZM_ERR_SIZMIS_VEC    "size mismatch of vectors"
#define ZM_ERR_SIZMIS_MAT    "size mismatch of matrices"
#define ZM_ERR_SIZMIS_MATVEC "size mismatch of matrix and vector"
#define ZM_ERR_NONSQR_MAT    "not a square matrix"

#define ZM_ERR_SIZUNFOUND_VEC "vector size not specified"
#define ZM_ERR_SIZUNFOUND_MAT "matrix size not specified"

#define ZM_ERR_INV_ROW       "invalid row specified"
#define ZM_ERR_INV_COL       "invalid column specified"
#define ZM_ERR_INV_INDEX     "invalid index specified"

#define ZM_ERR_ZERONORM      "cannot normalize zero vector"

#define ZM_ERR_CVEC_CONJPAIR_UNABLE "unable to pair co-conjugate complex numbers."

#define ZM_ERR_LE_SINGULAR   "matrix is singular"
#define ZM_ERR_LE_CHOLESKY   "matrix is not positive-semidefinite"

#define ZM_ERR_LOG_INVALID   "invalid base of logarithm"

#define ZM_ERR_PE_DEFL       "equation deflated"

#define ZM_ERR_SEQ_DT_UNFOUND "delta time not specified"

#define ZM_ERR_PEX_DIMUNFOUND "dimension not specified"
#define ZM_ERR_PEX_INVDIM    "invalid dimension of a polynomial"
#define ZM_ERR_PEX_DIMMIS    "dimension mismatch"
#define ZM_ERR_PEX_EQ_SIZMIS "shortage of vector size %d for answers of a %s-order polynomial equation."

#define ZM_ERR_IP_INVTERM    "invalid term for polynomial curve"
#define ZM_ERR_IP_SIZMIS     "size mismatch of sample point vector"
#define ZM_ERR_IP_INVTYPE    "unknown edge type"

#define ZM_ERR_NURBS_INVORDER     "invalid order for NURBS"
#define ZM_ERR_NURBS_INVDIFFORDER "invalid order of differentiation"

#define ZM_ERR_INVALID_NUMSAMP  "too many samples %d required out of %d"

#define ZM_ERR_OPT_NOEVAL       "no evaluator assigned"
#define ZM_ERR_OPT_INI          "unable to set the initial point"
#define ZM_ERR_OPT_STEP         "unable to compute step length"
#define ZM_ERR_OPT_UNSOLVE      "optimal solution doesn't exist"
#define ZM_ERR_OPT_INF          "infinite solution exists"
#define ZM_ERR_OPT_NONCONVEX    "cannot solve non-convex QP programs"

#define ZM_ERR_MCA_INVSIZ       "invalid size of vector for clustering: %d"
#define ZM_ERR_MCA_EMPTY        "empty set assigned"
#define ZM_ERR_MCA_NOCOREFUNC   "method not assigned to find core of a cluster"
#define ZM_ERR_MCA_NOERRORFUNC  "method not assigned to compute error between samples"
#define ZM_ERR_MCA_NODISTFUNC   "method not assigned to compute distance between samples"

#define ZM_ERR_GNG_NONEMPTY_UNIT "not an empty GNG unit"

#define ZM_ERR_OSCIL_INVTERM "invalid natural term of the oscillator"

#define ZM_ERR_FATAL         "fatal error! - please report to the author"

#endif /* __ZM_ERRMSG_H__ */
