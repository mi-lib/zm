/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_nm - optimization tools: Nelder-Mead's downhill simplex method.
 */

#ifndef __ZM_OPT_NM_H__
#define __ZM_OPT_NM_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief solve an optimization problem by Nelder-Mead method.
 *
 * zOptSolveNM() solves an optimization problem by Nelder-Mead method, which
 * is also known as the downhill simplex/polytope method proposed by Nelder
 * and Mead (1965).
 */
__ZM_EXPORT int zOptSolveNM(double (* f)(zVec,void*), void *util, zVec min, zVec max, int iter, double tol, zVec ans, double *eval);

__END_DECLS

#endif /* __ZM_OPT_NM_H__ */
