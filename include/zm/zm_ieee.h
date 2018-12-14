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
typedef union{ unsigned char c[8]; double d; } __ieee_fp_t;
/*! \endcond */

/*! \brief HUGE_VAL, if not conforming to C99. */
#ifndef HUGE_VAL
/*! \cond __huge_val only required for my HUGE_VAL */
extern __ieee_fp_t __huge_val;
/*! \endcond */
#define __ZM_NEED_HUGE_VAL
#define HUGE_VAL ( __huge_val.d )
#endif /* HUGE_VAL */

/*! \brief check if the value \a x is HUGE_VAL. */
__EXPORT int zIsInf(double x);

/* define NAN, if not conforming to C99. */
#ifndef NAN
/*! \cond __nan_val only required for my NAN */
extern __ieee_fp_t __nan_val;
/*! \endcond */
#define __ZM_NEED_NAN
#define NAN ( __nan_val.d )
#endif /* NAN */

/*! \brief check if the value \a x is NAN. */
__EXPORT int zIsNan(double x);

__END_DECLS

#endif /* __ZM_IEEE_H__ */
