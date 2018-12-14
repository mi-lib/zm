#include <zm/zm.h>

int main(void)
{
  double val;

  val = 1.0 / 0.0;
  zAssert( zIsInf, zIsInf( HUGE_VAL ) &&
                   zIsInf(-HUGE_VAL ) &&
                   zIsInf( log( 0.0 ) ) &&
                   zIsInf( val ) &&
                   !zIsInf( NAN ) &&
                   !zIsInf( 0 ) );
  zAssert( zIsNan, zIsNan( NAN ) &&
                   zIsNan( sqrt( -1 ) ) &&
                   zIsNan( log( -1 ) ) &&
                   !zIsNan( HUGE_VAL ) &&
                   !zIsNan( 0 ) );
  return EXIT_SUCCESS;
}
