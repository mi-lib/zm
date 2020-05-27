/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_complex_pe - complex number class: polynomial equation solver.
 */

#ifndef __ZM_COMPLEX_PE_H__
#define __ZM_COMPLEX_PE_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief solve quadratic and cubic equation.
 *
 * zQESolve() finds two answers including duplicate answer of
 * a quadratic equation \a a x^2 + \a b x + \a c = 0.
 * \a a must be nonzero.
 * The answers are put into an array pointed by \a ans.
 *
 * zCESolve() finds three answers including duplicate answer of
 * a cubic equation \a a x^3 + \a b x^2 + \a c x + \a d = 0
 * by Cardano's formula (originally derived by Fontana=Tartaglia).
 * \a a must be nonzero.
 * The answers are put into an array pointed by \a ans.
 * \return
 * zQESolve() and zCESolve() return \a ans if succeeding.
 * If \a a is zero, the null pointer is returned.
 * \notes
 * The size of the array pointed by \a ans must be larger than
 * 2 for zQESolve() and 3 for zCESolve(), respectively.
 * If not, anything might happen.
 */
__EXPORT zComplex *zQESolve(double a, double b, double c, zComplex ans[]);
__EXPORT zComplex *zCESolve(double a, double b, double c, double d, zComplex ans[]);

__END_DECLS

#endif /* __ZM_COMPLEX_PE_H__ */
