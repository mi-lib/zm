/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ieee - IEEE conforming irregular values.
 */

#include <zm/zm_misc.h>

/* define HUGE_VAL, if not. conforming to C99. */
#ifdef __ZM_NEED_HUGE_VAL
# include <endian.h>
__ieee_fp_t __huge_val = {
# if __BYTE_ORDER == __BIG_ENDIAN
  { 0x7f, 0xf0, 0, 0, 0, 0, 0, 0 }
# endif
# if __BYTE_ORDER == __LITTLE_ENDIAN
  { 0, 0, 0, 0, 0, 0, 0xf0, 0x7f }
# endif
};
#endif /* __ZM_NEED_HUGE_VAL */

/* zIsInf
 * - check if the value is HUGE_VAL.
 */
int zIsInf(double x)
{
  __ieee_fp_t val;

  val.d = x;
  return ( ( val.c[7] & 0x7f ) == 0x7f && val.c[6] == 0xf0 &&
           ( val.c[5] | val.c[4] | val.c[3] |
             val.c[2] | val.c[1] | val.c[0] ) == 0 ) ? 1 : 0;
}

/* define NAN, if not. conforming to C99. */
#ifdef __ZM_NEED_NAN
# include <endian.h>
__ieee_fp_t __nan_val = {
# if __BYTE_ORDER == __BIG_ENDIAN
  { 0x7f, 0xf8, 0, 0, 0, 0, 0, 0 }
# endif
# if __BYTE_ORDER == __LITTLE_ENDIAN
  { 0, 0, 0, 0, 0, 0, 0xf8, 0x7f }
# endif
};
#endif /* __ZM_NEED_NAN */

/* zIsNan
 * - check if the value is NAN.
 */
int zIsNan(double x)
{
  __ieee_fp_t val;

  val.d = x;
  return ( ( val.c[7] & 0x7f ) == 0x7f && ( val.c[6] & 0xf0 ) == 0xf0 &&
           ( ( val.c[6] & 0xf ) | val.c[5] | val.c[4] | val.c[3] |
             val.c[2] | val.c[1] | val.c[0] ) ) ? 1 : 0;
}
