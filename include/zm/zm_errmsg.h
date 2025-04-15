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

#define ZM_WARN_ITERATION                  "iteration unconverged within %d steps"

#define ZM_WARN_BESSEL_UNDEF               "undefined for negative integer"
#define ZM_WARN_VEC_INVALIDSIZE            "invalid vector size, ignored"
#define ZM_WARN_MAT_INVALIDSIZE            "invalid matrix size, ignored"

#define ZM_WARN_VEC_SIZEMISMATCH           "vector has only %d components, while specified size is %d"
#define ZM_WARN_MAT_SIZEMISMATCH           "matrix has only %d components, while specified size is %d"

#define ZM_WARN_VECLIST_EMPTY              "empty vector list assigned"

#define ZM_WARN_LE_ZEROPIVOT               "cannot sweep out by zero pivot"

#define ZM_WARN_NLE_NOROOT                 "probably, equation has no root"

#define ZM_WARN_SEQ_EMPTY                  "empty vector sequence assigned"
#define ZM_WARN_SEQ_STEP                   "too large step index"

#define ZM_WARN_PEX_SIZEMISMATCH           "polynomial expression has only %d components, while specified size is %d"

#define ZM_WARN_INSUFFICIENT_SAMPLES       "maybe insufficient number of samples %d for %d parameters"

#define ZM_WARN_OPT_BADINI                 "probably inadequate initial guess"

#define ZM_WARN_ODE_GEAR1                  "invalid step number %d, modified to 1"
#define ZM_WARN_ODE_GEAR2                  "step number %d over 6 may cause instability, modified"

#define ZM_WARN_GRAPH_DUPNODE              "duplicate specification of a node"
#define ZM_WARN_GRAPH_DUPCONNECTION        "duplicate node connection, overwritten"

#define ZM_WARN_MVA_NOSILHOUETTE           "silhouette score not computed, need to call zMClusterSilhouetteScore()"

/* error messages */

#define ZM_ERR_ZERODIV                     "cannot devide by tiny value"
#define ZM_ERR_OUTOFRANGE                  "specified index out of range"

#define ZM_ERR_VEC_NULL                    "null vector assigned"
#define ZM_ERR_VEC_SIZEMISMATCH            "size mismatch of vectors"
#define ZM_ERR_MAT_SIZEMISMATCH            "size mismatch of matrices"
#define ZM_ERR_MAT_SIZEMISMATCH_VEC        "size mismatch of matrix and vector"
#define ZM_ERR_MAT_SIZEMISMATCH_CVEC       "size mismatch of matrix and complex vector"
#define ZM_ERR_MAT_NOTSQR                  "not a square matrix"
#define ZM_ERR_MAT_NOTSYMMETRIC            "not a symmetric matrix"

#define ZM_ERR_VEC_SIZENOTFOUND            "vector size not specified"
#define ZM_ERR_MAT_SIZENOTFOUND            "matrix size not specified"

#define ZM_ERR_INVALID_ROW                 "invalid row specified"
#define ZM_ERR_INVALID_COL                 "invalid column specified"
#define ZM_ERR_INVALID_INDEX               "invalid index specified"

#define ZM_ERR_VEC_ZERONORM                "cannot normalize zero vector"

#define ZM_ERR_CVEC_CONJPAIR_UNABLE        "unable to pair co-conjugate complex numbers."

#define ZM_ERR_MAT_SINGULAR                "matrix is singular"
#define ZM_ERR_MAT_NOTPOSITIVESEMIDEFINITE "matrix is not positive-semidefinite"

#define ZM_ERR_INVALID_LOGBASE             "invalid base of logarithm"

#define ZM_ERR_FFT_SIZEMISMATCH_VEC        "size mismatch of data vectors %d and %d"
#define ZM_ERR_FFT_SIZEMISMATCH_MAT        "size mismatch of data matrices %dx%d and %dx%d"

#define ZM_ERR_PE_DEFLATED                 "equation deflated"

#define ZM_ERR_SEQ_DT_NOTFOUND             "delta time not specified"

#define ZM_ERR_PEX_DIMNOTFOUND             "dimension not specified"
#define ZM_ERR_PEX_INVALID_DIM             "invalid dimension of a polynomial"
#define ZM_ERR_PEX_DIMMISMATCH             "dimension mismatch"
#define ZM_ERR_PEX_EQ_SIZEMISMATCH         "shortage of vector size %d for answers of a %s-order polynomial equation."

#define ZM_ERR_IP_INVALID_TERM             "invalid term for polynomial curve"
#define ZM_ERR_IP_SIZEMISMATCH             "size mismatch of sample point vector"
#define ZM_ERR_IP_INVALID_EDGETYPE         "unknown edge type"

#define ZM_ERR_IP_TRVELPROF_NEGATIVEDISTANCE   "cannot travel negative distance %g"
#define ZM_ERR_IP_TRVELPROF_NONPOSITIVEMAXVEL  "not positive maximum velocity %g"
#define ZM_ERR_IP_TRVELPROF_INVALID_INITIALVEL "initial velocity %g larger than the maximum velocity %g"
#define ZM_ERR_IP_TRVELPROF_INVALID_TERMVEL    "terminal velocity %g larger than the maximum velocity %g"
#define ZM_ERR_IP_TRVELPROF_NONPOSITIVEACC     "not positive acceleration %g"

#define ZM_ERR_NURBS_INVALID_ORDER         "invalid order for NURBS %d was specified, has to be less than %d"
#define ZM_ERR_NURBS_INVALID_DIFFNUM       "invalid number of differentiation %d was specified, has to be within [0-%d]"
#define ZM_ERR_NURBS_ORDERMISMATCH         "size mismatch of B-spline orders (%d - %d)"
#define ZM_ERR_NURBS_KNOTNUMMISMATCH       "mismatch of number of knots of B-spline parameters (%d - %d)"
#define ZM_ERR_NURBS_KNOTALREADY           "knot already allocated"
#define ZM_ERR_NURBS_CPALREADY             "control point already allocated"
#define ZM_ERR_NURBS_INVALID_CP            "invalid index of control point specified"
#define ZM_ERR_NURBS_INVALID_KNOTSIZE      "the size of knot is %d, must be larger than %d"

#define ZM_ERR_INVALID_NUMSAMP             "too many samples %d required out of %d"

#define ZM_ERR_OPT_NOEVALUATOR             "no evaluator assigned"
#define ZM_ERR_OPT_INI                     "unable to set the initial point"
#define ZM_ERR_OPT_STEP                    "unable to compute step length"
#define ZM_ERR_OPT_UNSOLVABLE              "optimal solution doesn't exist"
#define ZM_ERR_OPT_INF                     "infinite solution exists"
#define ZM_ERR_OPT_NONCONVEX               "cannot solve non-convex QP programs"

#define ZM_ERR_MVA_INVALID_SIZE            "invalid size of vector for clustering: %d"
#define ZM_ERR_MVA_EMPTYSET                "empty set assigned"
#define ZM_ERR_MVA_NOCOREFUNC              "method not assigned to find core of a cluster"
#define ZM_ERR_MVA_NOERRORFUNC             "method not assigned to compute error between samples"
#define ZM_ERR_MVA_NODISTFUNC              "method not assigned to compute distance between samples"

#define ZM_ERR_GNG_NONEMPTY_UNIT           "not an empty GNG unit"

#define ZM_ERR_OSCIL_INVALID_TERM          "invalid natural term of the oscillator"

#define ZM_ERR_FATAL                       "fatal error! - please report to the author"

#endif /* __ZM_ERRMSG_H__ */
