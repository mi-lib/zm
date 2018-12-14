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

#define ZM_WARN_LE_ZEROPIVOT "cannot sweep out by zero pivot"

#define ZM_WARN_NLE_NOROOT   "probably, equation has no root"

#define ZM_WARN_SEQ_EMPTY    "sequence list is empty"
#define ZM_WARN_SEQ_STEP     "too large step index"

#define ZM_WARN_OPT_BADINI   "probably inadequate initial guess"

#define ZM_WARN_ODE_GEAR1    "invalid step number %d, modified to 1"
#define ZM_WARN_ODE_GEAR2    "step number %d over 6 may cause instability, modified"

#define ZM_WARN_GRAPH_DUPSPC "duplicate specification of a node"
#define ZM_WARN_GRAPH_DUPCON "duplicate node connection, overwritten"

/* error messages */

#define ZM_ERR_ZERODIV       "cannot devide by tiny value"

#define ZM_ERR_STAT_ILLS     "illegal series assigned"

#define ZM_ERR_NULLVEC       "null vector assigned"
#define ZM_ERR_SIZMIS_VEC    "size mismatch of vector"
#define ZM_ERR_SIZMIS_MAT    "size mismatch of matrix"
#define ZM_ERR_SIZMIS_MATVEC "size mismatch of matrix and vector"
#define ZM_ERR_NONSQR_MAT    "not a square matrix"

#define ZM_ERR_INV_ROW       "invalid row specified"
#define ZM_ERR_INV_COL       "invalid column specified"
#define ZM_ERR_INV_INDEX     "invalid index specified"

#define ZM_ERR_ZERONORM      "cannot normalize zero vector"

#define ZM_ERR_LE_SINGULAR   "matrix is singular"
#define ZM_ERR_LE_CHOLESKY   "matrix is not positive-semidefinite"

#define ZM_ERR_LOG_INVALID   "invalid base of logarithm"

#define ZM_ERR_PE_DEFL       "equation deflated"

#define ZM_ERR_PEX_DIMMIS    "dimension mismatch"

#define ZM_ERR_IP_INVTERM    "invalid term for polynomial curve"
#define ZM_ERR_IP_SIZMIS     "size mismatch of sample point vector"
#define ZM_ERR_IP_INVTYPE    "unknown edge type"

#define ZM_ERR_NURBS_INVDIM  "invalid curve dimension"
#define ZM_ERR_NURBS_INVODR  "invalid order of differentiation"

#define ZM_ERR_OPT_NOEVAL    "no evaluator assigned"
#define ZM_ERR_OPT_INI       "unable to set the initial point"
#define ZM_ERR_OPT_STEP      "unable to compute step length"
#define ZM_ERR_OPT_UNSOLVE   "optimal solution doesn't exist"
#define ZM_ERR_OPT_INF       "infinite solution exists"

#define ZM_ERR_MCA_NOMEAN    "mean computation method not assigned"

#define ZM_ERR_OSCIL_INVTERM "invalid natural term of oscillator"

#define ZM_ERR_FATAL         "fatal error! - please report to the author"

#endif /* __ZM_ERRMSG_H__ */
