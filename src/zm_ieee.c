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
  __ieee_fp_t val;

  val.d = x;
  if( val.c[6] != 0xf0 ||
      ( val.c[5] | val.c[4] | val.c[3] | val.c[2] | val.c[1] | val.c[0] ) != 0 )
    return 0;
  return val.c[7] == 0x7f ? 1 : ( val.c[7] == 0xff ? -1 : 0 );
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
  __ieee_fp_t val;

  val.d = x;
  return ( ( val.c[7] & 0x7f ) == 0x7f && ( val.c[6] & 0xf0 ) == 0xf0 &&
           ( ( val.c[6] & 0xf ) | val.c[5] | val.c[4] | val.c[3] |
             val.c[2] | val.c[1] | val.c[0] ) ) ? 1 : 0;
}

/* check if the value is a finite number. */
int zIsFinite(double x)
{
  return zIsInf( x ) != 1 && zIsInf( x ) != -1 && !zIsNan( x ) ? 1 : 0;
}
