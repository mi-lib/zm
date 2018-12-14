#include <zm/zm_vec.h>

int main(void)
{
  zVec v;
  double val;
  int im;

  zVecWrite( ( v = zVecCreateList( 10,
    1.0, 3.0, 5.0, -3.0, -4.0, 2.0, 1.0, -2.0, 4.0, 0.0 ) ) );
  val = zVecMax( v, &im );
  printf( "max = %f (%d)\n", val, im );
  val = zVecMin( v, &im );
  printf( "min = %f (%d)\n", val, im );
  val = zVecAbsMax( v, &im );
  printf( "|max| = %f (%d)\n", val, im );
  val = zVecAbsMin( v, &im );
  printf( "|min| = %f (%d)\n", val, im );
  return 0;
}
