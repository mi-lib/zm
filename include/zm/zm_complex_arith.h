/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_complex_arith - complex number class: arithmetics.
 */

#ifndef __ZM_COMPLEX_ARITH_H__
#define __ZM_COMPLEX_ARITH_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief squared absolute value of a complex number
 * \retval the squared absolute value of a complex number \a c.
 */
__EXPORT double zComplexSqrAbs(zComplex *c);

/*! \brief absolute value of a complex number
 * \retval the absolute value of a complex number \a c.
 */
#define zComplexAbs(c) sqrt(zComplexSqrAbs(c))

/*! \brief argument angle of a complex number
 * \retval the argument angle of a complex number \a c,
 * which is in radians between -PI and PI.
 */
__EXPORT double zComplexArg(zComplex *c);

/*! \brief the complex conjugate of a complex number.
 * zComplexConj() computes the complex conjugate of a complex
 * number \a c. The result is stored where is pointed by \a cc.
 * \retval a pointer \a cc.
 */
__EXPORT zComplex *zComplexConj(zComplex *c, zComplex *cc);

/* 'zComplexAdd()' adds the 2 given complex numbers,
 * 'c1' and 'c2'. The result is put into 'c'.
 * #
 * 'zComplexSub()' subtracts the given complex number 'c2' from the
 * other 'c1'. The result is put into 'c'.
 * #
 * 'zComplexRev()' reverses the given complex number 'c'.
 * The result is put into 'c'.
 * #
 * 'zComplexMul()' multiply the given complex number 'c'
 * by a real number 'k'. The result is put into 'ec'.
 * #
 * 'zComplexDiv()' divide the given complex number 'c' by a real
 * number 'k'. The result is put into 'rc'.
 * #
 * 'zComplexCMul()' multiplies the 2 given complex number 'c1' and
 * 'c2'. The result is put into 'c'.
 * #
 * 'zComplexCDiv()' divides the given complex number 'c1' by the
 * other 'c2'. The result is put into 'c'.
 * #
 * 'zComplexPow()' calculates the complex number 'c' raised to the
 * power of 'z', which is a real number. The result is put into 'pc'.
 * #
 * 'zComplexCPow()' calculates the complex number 'c' raised to the
 * power of 'z', which is a complex number. The result is put into
 * 'pc'.
 * #
 * 'zComplexLog()' calculates the base 'base' - a real number -
 * logarithm of the given complex number 'c'. The result is put into
 * 'lc'.
 * #
 * 'zComplexCLog()' calculates the base 'base' - a complex number -
 * logarithm of the given complex number 'c'. The result is put into
 * 'lc'.
 * [RETURN VALUE]
 * Each of all these functions returns a pointer to the result.
 * #
 * 'zComplexDiv()' returns the null pointer if 'k' is 0.
 * [NOTES]
 * A number raised to the power of a certain complex number, or a
 * logarithm of a certain complex number  could not be determined
 * uniquely indeed, because one complex has non-unique argument
 * values. 'zComplexPow()', 'zComplexCPow()', 'zComplexLog()' and
 * 'zComplexCLog()' calculate only a solution with the argument
 * value of 'c' between -PI and PI.
 */
__EXPORT zComplex *zComplexAdd(zComplex *c1, zComplex *c2, zComplex *c);
__EXPORT zComplex *zComplexSub(zComplex *c1, zComplex *c2, zComplex *c);
__EXPORT zComplex *zComplexRev(zComplex *c, zComplex *rc);
__EXPORT zComplex *zComplexMul(zComplex *c, double k, zComplex *ec);
__EXPORT zComplex *zComplexDiv(zComplex *c, double k, zComplex *rc);
__EXPORT zComplex *zComplexCMul(zComplex *c1, zComplex *c2, zComplex *c);
__EXPORT zComplex *zComplexCDiv(zComplex *c1, zComplex *c2, zComplex *c);
__EXPORT zComplex *zComplexPow(zComplex *c, double z, zComplex *pc);
__EXPORT zComplex *zComplexPowRef(zComplex *c, double z, zComplex *ref, zComplex *pc);
__EXPORT zComplex *zComplexCPow(zComplex *c, zComplex *z, zComplex *pc);
__EXPORT zComplex *zComplexLog(zComplex *c, double base, zComplex *lc);
__EXPORT zComplex *zComplexCLog(zComplex *c, zComplex *base, zComplex *lc);

/*! \brief normalize a complex number.
 *
 * zComplexNormalize() normalizes a complex number \a c so as
 * to magnify its absolute value to be 1.
 * The result is stored where is pointed by \a nc.
 * \retval a pointer \a nc
 * \retval the null pointer if the absolute value of \a c is 0
 */
__EXPORT zComplex *zComplexNormalize(zComplex *c, zComplex *nc);

__END_DECLS

#endif /* __ZM_COMPLEX_ARITH_H__ */
