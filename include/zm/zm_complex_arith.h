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

/*! \brief arithmatics of complex numbers.
 *
 * zComplexAdd() adds two complex numbers \a c1 and \a c2.
 * The result is put into \a c.
 *
 * zComplexSub() subtracts a complex number \a c2 from another
 * \a c1. The result is put into \a c.
 *
 * zComplexRev() reverses a complex number \a c.
 * The result is put into \a rc.
 *
 * zComplexMul() multiplies a complex number \a c by a real
 * number \a k. The result is put into \a ec.
 *
 * zComplexDiv() divides a complex number \a c by a real number
 * \a k. The result is put into \a rc.
 *
 * zComplexCMul() multiplies a complex number \a c1 by another
 * \a c2. The result is put into \a c.
 *
 * zComplexCMulConj() multiplies a complex number \a c1 by the
 * conjugate of another \a c2. The result is put into \a c.
 *
 * zComplexCDiv() divides a complex number \a c1 by \another
 * \a c2. The result is put into \a c.
 *
 * zComplexPow() calculates a complex number \a c raised to the
 * power of a real number \a z. The result is put into \a pc.
 *
 * zComplexCPow() calculates a complex number \a c raised to the
 * power of another complex number \a z. The result is put into
 * \a pc.
 *
 * zComplexLog() calculates the base \a base (a real number)
 * logarithm of a complex number \a c. The result is put into
 * \a lc.
 *
 * zComplexCLog() calculates the base \a base (a complex number)
 * logarithm of another complex number \a c. The result is put
 * into \a lc.
 * \return
 * These functions returns a pointer to the result.
 *
 * zComplexDiv() returns the null pointer if \a k is zero.
 * \notes
 * A number raised to the power of a complex number and a logarithm
 * of a complex number are indeed not unique because one complex
 * has non-unique argument angles.
 * zComplexPow(), zComplexCPow(), zComplexLog() and zComplexCLog()
 * calculate only one solution with the argument angle of \a c
 * between -pi and pi.
 */
__EXPORT zComplex *zComplexAdd(zComplex *c1, zComplex *c2, zComplex *c);
__EXPORT zComplex *zComplexSub(zComplex *c1, zComplex *c2, zComplex *c);
__EXPORT zComplex *zComplexRev(zComplex *c, zComplex *rc);
__EXPORT zComplex *zComplexMul(zComplex *c, double k, zComplex *ec);
__EXPORT zComplex *zComplexDiv(zComplex *c, double k, zComplex *rc);
__EXPORT zComplex *zComplexCMul(zComplex *c1, zComplex *c2, zComplex *c);
__EXPORT zComplex *zComplexCMulConj(zComplex *c1, zComplex *c2, zComplex *c);
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
