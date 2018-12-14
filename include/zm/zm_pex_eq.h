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
 * zPexBHDST(), zPexBH() and zPexDKA() solve a
 * polynomial equation 'a'=0. 'zPexBHDST()' and 'zPexBH()'
 * are based on Bairstow-Hitchcock's method, while
 * 'zPexDKA()' is based on Durand-Kerner-Aberth's method.
 * #
 * All answers are put into the array pointed by 'ans'
 * ('ans' should have enough size) as complex numbers.
 * 'tol' is a tolerance to finish iteration.
 * 'iter' is the maximum number of iteration. When it
 * does iteration over 'iter' times, it is aborted.
 * If zero is given for 'iter', Z_MAX_ITER_NUM defined
 * in 'zm_misc.h' is chosen instead.
 * 'zPexBHDST()' destroys 'a' in the course of computation.
 * [KNOWN PROBLEMS]
 * Since Bairstow-Hitchcock's method is an iteration,
 * 'zPexBH()' may not converge due to bad initialization.
 * 'zPexDKA()' is slightly advantageous from this viewpoint.
 * \return
 * 'zPexBHDST()' and 'zPexBH()' return a pointer 'ans'
 * if succeeding. Or, the null pointer is returned when
 * failing to allocate working memory.
 */
#define ZM_PEX_EQ_TOL ( 1.0e-7 )
__EXPORT zComplex *zPexBHDST(zPex a, zComplex *ans, double tol, int iter);
__EXPORT zComplex *zPexBH(zPex a, zComplex *ans, double tol, int iter);
__EXPORT zComplex *zPexDKA(zPex a, zComplex *ans, double tol, int iter);

__END_DECLS

#endif /* __ZM_PEX_EQ_H__ */
