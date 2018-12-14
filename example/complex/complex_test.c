#include <zm/zm_complex.h>

int main(void)
{
  zComplex c1, c2, c3, c;

  zComplexCreate( &c1, 1, 0 );
  printf( "C1 = " ); zComplexWrite( &c1 ); zEndl();
  printf( " abs = %f\n", zComplexAbs(&c1) );
  printf( " arg = %f\n", zRad2Deg(zComplexArg(&c1)) );

  zComplexPolar( &c2, 1, zDeg2Rad(45) );
  printf( "C2 = " ); zComplexWrite( &c2 ); zEndl();
  printf( " abs = %f\n", zComplexAbs(&c2) );
  printf( " arg = %f\n", zRad2Deg(zComplexArg(&c2)) );

  zComplexConj( &c2, &c3 );
  printf( "C2* = " ); zComplexWrite( &c3 ); zEndl();
  zComplexCMul( &c2, &c3, &c );
  printf( "C2.C3 = " ); zComplexWrite( &c ); zEndl();
  zComplexCDiv( &c2, &c3, &c );
  printf( "C2/C3 = " ); zComplexWrite( &c ); zEndl();

  zComplexCreate( &c1, 0.5, 0.5*sqrt(3) );
  printf( "C1 = " ); zComplexWrite( &c1 ); zEndl();

  zComplexPow( &c1, 6, &c3 );
  printf( "C1^6 = " ); zComplexWrite( &c3 ); zEndl();

  zComplexCreate( &c, 1, 1 );
  printf( "C = " ); zComplexWrite( &c ); zEndl();
  zComplexCPow( &c1, &c, &c2 );
  printf( "C1^C = " ); zComplexWrite( &c2 ); zEndl();

  zComplexCLog( &c2, &c1, &c3 );
  printf( "log_C1 C2 = " ); zComplexWrite( &c3 ); zEndl();

  zComplexNormalize( &c, &c );
  printf( "C/|C| = " ); zComplexWrite( &c ); zEndl();

  return 0;
}

/*
 * valid output:
 *
 * 1
 *  abs = 1
 *  arg = 0
 * 0.7070 + 0.7070 i
 *  abs = 1
 *  arg = 45
 * 0.7070 - 0.7070 i
 * 1
 * 0 + i
 * 0.5 + 0.8660 i
 * 1
 * 0.1757 + 0.3039 i
 * 1 + i
 * 0.7070 + 0.7070 i
 */
