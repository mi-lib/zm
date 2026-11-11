/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ieee - IEEE conforming irregular values.
 */

#include <zm/zm_misc.h>

/* HUGE_VAL conforming to IEEE 754 floating point. */
const __ieee_fp_t zm_huge_val = {
# if __BYTE_ORDER == __BIG_ENDIAN
  { 0x7f, 0xf0, 0, 0, 0, 0, 0, 0 }
# endif
# if __BYTE_ORDER == __LITTLE_ENDIAN
  { 0, 0, 0, 0, 0, 0, 0xf0, 0x7f }
# endif
};

/* check if the value is infinity. */
int zIsInf(double x)
{
  return _zIsInf( x );
}

/* Not-a-number conforming to IEEE 754 floating point. */
const __ieee_fp_t zm_nan_val = {
# if __BYTE_ORDER == __BIG_ENDIAN
  { 0x7f, 0xf8, 0, 0, 0, 0, 0, 0 }
# endif
# if __BYTE_ORDER == __LITTLE_ENDIAN
  { 0, 0, 0, 0, 0, 0, 0xf8, 0x7f }
# endif
};

/* check if the value is not-a-number. */
int zIsNan(double x)
{
  return _zIsNan( x );
}

/* check if the value is a finite number. */
int zIsFinite(double x)
{
  return _zIsFinite( x );
}
