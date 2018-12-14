/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_nle_ss - nonlinear equation: successive substitution method.
 */

#ifndef __ZM_NLE_SS_H__
#define __ZM_NLE_SS_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief successive substitution method.
 *
 * zSSSolve() solves an equation \a x=\a f(\a x) by
 * successive substitution method. The iteration will
 * converge when the maximum eigenvalue of the Jacobian
 * matrix of \a f is less than 1.
 * \a iter is the maximum iteration number. If zero is
 * given for \a iter, Z_MAX_ITER_NUM is applied instead.
 * The function forms as \a f(\a x, \a y, \a util),
 * where \a util is a programmer's utility, and the
 * function value is stored into \a y.
 * \return
 * a pointer \a x is returned.
 */
__EXPORT zVec zSSSolve(zVec (* f)(zVec,zVec,void*) ,zVec x, void *util, int iter);

__END_DECLS

#endif /* __ZM_NLE_SS_H__ */
