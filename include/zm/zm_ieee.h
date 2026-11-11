/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_ieee.h
 * \brief IEEE conforming irregular values.
 * \author Zhidao
 */

#ifndef __ZM_IEEE_H__
#define __ZM_IEEE_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \cond __ieee_fp_t only required for my HUGE_VAL and NAN */
typedef union{ ubyte c[8]; double d; } __ieee_fp_t;
/*! \endcond */

/*! \brief the \a i-th byte of an IEEE 754 double-precision floating-point number \a x.
 *
 * The bytes are indexed in the little endian order regardless of the byte order of the machine,
 * namely, the 7th byte includes the sign bit and the upper seven bits of the exponent part, the 6th
 * byte includes the lower four bits of the exponent part and the upper four bits of the mantissa part,
 * and the 0th to 5th bytes include the rest of the mantissa part.
 *
 * Since it refers to the address of \a x, \a x has to be an lvalue of double type, i.e., a variable
 * (not an expression nor a constant) that is not declared with the register storage class. The same
 * applies to all the following macros, some of which refer to the address of \a x more than once,
 * so that \a x should not include side effects such as a[i++]. Those macros do not depend on the
 * compiler options for floating-point arithmetics, since they only consist of integer operations.
 */
#if __BYTE_ORDER == __BIG_ENDIAN
#define _zIEEEDoubleByte(x,i) ( ((const ubyte *)&(x))[7-(i)] )
#else
#define _zIEEEDoubleByte(x,i) ( ((const ubyte *)&(x))[i] )
#endif

/*! \brief check if the exponent bits of an IEEE 754 double-precision floating-point number \a x are
 * all one, namely, \a x is either infinity or not-a-number. */
#define _zIEEEExpIsMax(x) ( ( _zIEEEDoubleByte(x,7) & 0x7f ) == 0x7f && ( _zIEEEDoubleByte(x,6) & 0xf0 ) == 0xf0 )

/*! \brief check if the mantissa bits of an IEEE 754 double-precision floating-point number \a x are
 * all zero. */
#define _zIEEEMantIsZero(x) ( ( ( _zIEEEDoubleByte(x,6) & 0x0f ) | _zIEEEDoubleByte(x,5) | _zIEEEDoubleByte(x,4) | \
  _zIEEEDoubleByte(x,3) | _zIEEEDoubleByte(x,2) | _zIEEEDoubleByte(x,1) | _zIEEEDoubleByte(x,0) ) == 0 )

/*! \brief check if the sign bit of an IEEE 754 double-precision floating-point number \a x is one. */
#define _zIEEESignIsNeg(x) ( ( _zIEEEDoubleByte(x,7) & 0x80 ) != 0 )

/*! \brief HUGE_VAL conforming to IEEE 754 floating point. */
__ZM_EXPORT const __ieee_fp_t zm_huge_val;
#ifndef HUGE_VAL
#define HUGE_VAL ( zm_huge_val.d )
#endif /* HUGE_VAL */

/*! \brief check if the value \a x is infinity.
 *
 * _zIsInf() is a macro version of zIsInf(), which inspects the bytes of \a x directly without a
 * function call. \a x has to be an lvalue of double type without side effects (see _zIEEEDoubleByte()).
 * \return
 * zIsInf() and _zIsInf() return 1 if \a x is the positive infinity, -1 if \a x is the negative
 * infinity, and 0 otherwise.
 */
#define _zIsInf(x) ( _zIEEEExpIsMax(x) && _zIEEEMantIsZero(x) ? ( _zIEEESignIsNeg(x) ? -1 : 1 ) : 0 )
__ZM_EXPORT int zIsInf(double x);

/*! \brief Not-a-number conforming to IEEE 754 floating point. */
__ZM_EXPORT const __ieee_fp_t zm_nan_val;
#ifndef NAN
#define NAN ( zm_nan_val.d )
#endif /* NAN */

/*! \brief check if the value \a x is not-a-number.
 *
 * _zIsNan() is a macro version of zIsNan(), which inspects the bytes of \a x directly without a
 * function call. \a x has to be an lvalue of double type without side effects (see _zIEEEDoubleByte()).
 * \return
 * zIsNan() and _zIsNan() return 1 if \a x is not-a-number, and 0 otherwise.
 */
#define _zIsNan(x) ( _zIEEEExpIsMax(x) && !_zIEEEMantIsZero(x) ? 1 : 0 )
__ZM_EXPORT int zIsNan(double x);

/*! \brief check if the value is a finite number.
 *
 * zIsFinite() checks if a value \a x is a finite number, namely, neither infinity nor not-a-number.
 *
 * _zIsFinite() is a macro version of zIsFinite(), which inspects the bytes of \a x directly without
 * a function call. \a x has to be an lvalue of double type (see _zIEEEDoubleByte()).
 * \return
 * zIsFinite() and _zIsFinite() return 1 if \a x is a finite number, and 0 otherwise.
 */
#define _zIsFinite(x) ( _zIEEEExpIsMax(x) ? 0 : 1 )
__ZM_EXPORT int zIsFinite(double x);

__END_DECLS

#endif /* __ZM_IEEE_H__ */
