/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_pex_eq - polynomial expression class: polynomial equation solver.
 */

#ifndef __ZM_PEX_EQ_H__
#define __ZM_PEX_EQ_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief numerical solution of polynomial equation.
 *
 * zPexBHDST(), zPexBH() and zPexDKA() solve a polynomial equation
 * \a a=0 based on Bairstow-Hitchcock's method and Durand-Kerner-Aberth's
 * method, respectively.
 *
 * All answers are put into an array pointed by \a ans as complex numbers.
 * \a ans must have enough size.
 * \a tol is a tolerance to terminate iteration.
 * \a iter is the maximum number of iteration. When it runs over \a iter
 * times, it is aborted. If zero is given for \a iter, Z_MAX_ITER_NUM
 * defined in zm_misc.h is used instead.
 * zPexBHDST() is a destructive version of zPexBD(), which destroys \a a
 * in the computation process.
 * \note
 * Since Bairstow-Hitchcock's method is an iteration, zPexBH() may not
 * converge due to bad initialization.
 * zPexDKA() is advantageous to it from this viewpoint.
 * \return
 * zPexBHDST() and zPexBH() return a pointer \a ans if succeeding.
 * Otherwise, the null pointer is returned.
 */
#define ZM_PEX_EQ_TOL ( 1.0e-7 )
__ZM_EXPORT zCVec zPexBHDST(zPex a, zCVec ans, double tol, int iter);
__ZM_EXPORT zCVec zPexBH(zPex a, zCVec ans, double tol, int iter);
__ZM_EXPORT zCVec zPexDKA(zPex a, zCVec ans, double tol, int iter);

__END_DECLS

#endif /* __ZM_PEX_EQ_H__ */
