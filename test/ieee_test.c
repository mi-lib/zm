#include <zm/zm.h>

int main(void)
{
#ifndef __WINDOWS__
  double val;

  val = 1.0 / 0.0;
  zAssert( zIsInf, zIsInf( HUGE_VAL ) == 1 &&
                   zIsInf(-HUGE_VAL ) == -1 &&
                   zIsInf( log( 0.0 ) ) == -1 &&
                   zIsInf( val ) == 1 &&
                   zIsInf( NAN ) == 0 &&
                   zIsInf( 0 ) == 0 );
  zAssert( zIsNan, zIsNan( NAN ) &&
                   zIsNan( sqrt( -1 ) ) &&
                   zIsNan( log( -1 ) ) &&
                   !zIsNan( HUGE_VAL ) &&
                   !zIsNan( 0 ) );
  zAssert( zIsFinite, zIsFinite( 1.0 ) &&
                      zIsFinite(-1.0 ) &&
                      !zIsFinite( HUGE_VAL ) &&
                      !zIsFinite(-HUGE_VAL ) &&
                      !zIsFinite( sqrt( -1 ) ) );
#endif /* __WINDOWS__ */
  return EXIT_SUCCESS;
}
