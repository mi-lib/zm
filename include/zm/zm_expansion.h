/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_expansion.h
 * \brief arithmetics of expansions.
 * \author Zhidao
 */

#ifndef __ZM_EXPANSION_H__
#define __ZM_EXPANSION_H__

#include <zm/zm_misc.h>

__BEGIN_DECLS

/*! \brief arithmetics of expansions.
 *
 * An expansion is an array of floating-point numbers called components, which represents their sum
 * exactly. The components are nonoverlapping, i.e., the highest nonzero bit of a component is lower
 * than the lowest nonzero bit of the next one, and are sorted in increasing order of magnitude.
 * A zero expansion is represented by a single zero component. Expansions produced by the following
 * functions have no zero components except this case. The sign of such an expansion is given by
 * zSgnExpansion().
 * They enable to compute the exact signs of polynomials of floating-point numbers, e.g., determinants.
 * See J. R. Shewchuk, Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric
 * Predicates, Discrete & Computational Geometry, 18(3):305-363, 1997.
 *
 * zGrowExpansion() adds a floating-point number \a b to an expansion \a e with \a n components, and
 * puts the result into \a h. \a h has to have at least \a n+1 elements. \a h can be the same with
 * \a e.
 *
 * zScaleExpansion() multiplies an expansion \a e with \a n components by a floating-point number
 * \a b, and puts the result into \a h. \a h has to have at least 2\a n elements, and has to be
 * different from \a e.
 *
 * zAddExpansionTripleProd() adds the product of three floating-point numbers \a a, \a b and \a c to an
 * expansion \a e with \a n components. \a e has to have at least \a n+4 elements. Note that the
 * product can be subtracted by giving -\a a instead of \a a, since the negation is exact.
 *
 * \a n can be zero for an empty expansion in all functions.
 * The same assumptions with zAdd() and zMul() are made. In particular, the results are exact unless
 * overflow or underflow happens.
 * \return
 * zGrowExpansion() and zScaleExpansion() return the number of components of \a h.
 * zAddExpansionTripleProd() returns the number of components of the updated \a e.
 */
__ZM_EXPORT int zGrowExpansion(int n, const double e[], double b, double h[]);
__ZM_EXPORT int zScaleExpansion(int n, const double e[], double b, double h[]);
__ZM_EXPORT int zAddExpansionTripleProd(int n, double e[], double a, double b, double c);

/*! \brief sign of an expansion.
 *
 * zSgnExpansion() returns the sign of an expansion \a e with \a n components, namely, 1 if it is
 * positive, -1 if negative, and 0 if zero or empty (\a n = 0). \a e must not have zero components
 * except a single zero component, which is satisfied by expansions produced by the above functions.
 */
#define zSgnExpansion(n,e) ( (n) > 0 ? zSgn( (e)[(n)-1] ) : 0 )

__END_DECLS

#endif /* __ZM_EXPANSION_H__ */
