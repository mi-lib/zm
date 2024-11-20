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

/*! \brief HUGE_VAL conforming to IEEE 754 floating point. */
__ZM_EXPORT const __ieee_fp_t zm_huge_val;
#ifndef HUGE_VAL
#define HUGE_VAL ( zm_huge_val.d )
#endif /* HUGE_VAL */

/*! \brief check if the value \a x is infinity. */
__ZM_EXPORT int zIsInf(double x);

/*! \brief Not-a-number conforming to IEEE 754 floating point. */
__ZM_EXPORT const __ieee_fp_t zm_nan_val;
#ifndef NAN
#define NAN ( zm_nan_val.d )
#endif /* NAN */

/*! \brief check if the value \a x is not-a-number. */
__ZM_EXPORT int zIsNan(double x);

/*! \brief check if the value is a finite number. */
__ZM_EXPORT int zIsFinite(double x);

__END_DECLS

#endif /* __ZM_IEEE_H__ */
