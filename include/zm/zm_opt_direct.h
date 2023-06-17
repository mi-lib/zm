/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_direct - optimization tools: DIRECT method.
 */

#ifndef __ZM_OPT_DIRECT_H__
#define __ZM_OPT_DIRECT_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* DIviding RECTangle method by Jones, Perttunen and Stuckman (1993).
 * reference:
 * D. R. Jones, C. D. Perttunen and B. E. Stuckman, Lipschizian Optimization Without
 * the Lipschitz Constant, Journal of Optimization Theory and Application, Vol. 79,
 * No. 1, pp. 157--181, 1993.
 */
__ZM_EXPORT int zOptSolveDIRECT(double (* f)(zVec,void*), void *util, zVec min, zVec max, int iter, double tol, zVec ans, double *eval);

#define ZOPT_DIRECT_MAX_ITER_NUM 1000

__END_DECLS

#endif /* __ZM_OPT_DIRECT_H__ */
