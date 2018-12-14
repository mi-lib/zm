/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_complex_pe - complex number class: polynomial equation solver.
 */

#ifndef __ZM_COMPLEX_PE_H__
#define __ZM_COMPLEX_PE_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* METHOD:
 * zQESolve, zCardano - solve quadratic and cubic equation.
 * [SYNOPSIS]
 * zComplex *zQESolve(double a, double b, double c, zComplex ans[]);
 * zComplex *zCardano(double a, double b, double c, double d, zComplex ans[]);
 * [DESCRIPTION]
 * 'zQESolve()' solves a quadratic equation
 *  'a' x^2 + 'b' x + 'c' = 0.
 * 'a' must be nonzero.
 * 2 answers of the quadratic equation (including duplicate
 * answer) are put into the array 'ans'.
 * #
 * 'zCardano()' solves a cubic equation
 *  'a' x^3 + 'b' x^2 + 'c' x + 'd' = 0
 * by Cardano s formula, which is originally derived by
 * Fontana=Tartaglia.
 * 'a' must be nonzero.
 * 3 answers of the equation (including duplicate answer) are
 * put into the array 'ans'.
 * [RETURN VALUE]
 * 'zQESolve()' and 'zCardano()' return the address 'ans',
 * when 'a' is nonzero, or the null pointer, otherwise.
 * [NOTES]
 * At least, the size of the array beginning with 'ans' must be
 * more than 2 and 3 for 'zQESolve()' and 'zCardano()',
 * respectively.
 * If size is not enough, anything might happen.
 */
__EXPORT zComplex *zQESolve(double a, double b, double c, zComplex ans[]);
__EXPORT zComplex *zCardano(double a, double b, double c, double d, zComplex ans[]);

__END_DECLS

#endif /* __ZM_COMPLEX_PE_H__ */
