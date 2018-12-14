/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_nle_se - nonlinear equation: single nonlinear equation solver.
 */

#ifndef __ZM_NLE_SE_H__
#define __ZM_NLE_SE_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* Bisec: bisection method */
__EXPORT double zNLE_Bisec(double (* f)(double,void*), double x1, double x2, void *priv, double tol, int iter);
/* Secant: Secant method */
__EXPORT double zNLE_Secant(double (* f)(double,void*), double x1, double x2, void *priv, double tol, int iter);
/* RF: Regula-Falsi method */
__EXPORT double zNLE_RF(double (* f)(double,void*), double x1, double x2, void *priv, double tol, int iter);
/* VDB: Van Wijngaarden-Dekker-Brent method */
__EXPORT double zNLE_VDB(double (* f)(double,void*), double x1, double x2, void *priv, double tol, int iter);

__END_DECLS

#endif /* __ZM_NLE_SE_H__ */
