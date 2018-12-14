#define DEBUG
#include <zm/zm_opt.h>

#define TEST 1

int main(void)
{
  zMat m;
  zVec q, z;

#if TEST == 1
  double marray[] = {
    4,-1, 1, 2, 1,
   -1, 4, 2, 1, 1,
   -1,-2, 0, 0, 0,
   -2,-1, 0, 0, 0,
   -1,-1, 0, 0, 0 };
  double qarray[] = { -1, -2, 1, 1, 10 };
  double ans[] = { 1.0/4.0, 3.0/8.0, 3.0/8.0, 0, 0 };
  int d = 5;

#elif TEST == 2
  double marray[] = {
    1,-1, 1,-1, 2,
   -1, 2, 1, 2, 1,
   -1,-1, 0, 0, 0,
    1,-2, 0, 0, 0,
   -2,-1, 0, 0, 0 };
  double qarray[] = { -2, -6, 2, 2, 3 };
  double ans[] = { 2.0/3.0, 4.0/3.0, 28.0/9, 4.0/9, 0 };
  int d = 5;

#elif TEST == 3
  double marray[] = {
    -2.0, 0.0, 2.5, 5.0, 3.0,
     0.0,-2.0, 5.0, 6.0, 2.0,
    -2.5,-5.0, 0.0, 0.0, 0.0,
    -5.0,-6.0, 0.0, 0.0, 0.0,
    -3.0,-2.0, 0.0, 0.0, 0.0 };
  double qarray[] = { -1.0, -1.0, 350, 450, 240 };
  double ans[] = { 80.0, 0.0, 0.0, 0.0, 161.0/3 };
  int d = 5;

#elif TEST == 4
  /* in this case, answer is [ 3  2 ] */
  double marray[] = {
    1, 0.5, 1, 1,
    0.5, 1, 2, 1,
    -1, -2, 0, 0,
    -1, -1, 0, 0 };
  double qarray[] = { -10, -11, 7, 5 };
  double ans[] = { 3.0, 2.0, 3.0/2, 9.0/2 };
  int d = 4;

#elif TEST == 5
  double marray[] = {
    1,-1,
    1, 0 };
  double qarray[] = { -2, 20 };
  double ans[] = { 2.0, 0.0 };
  int d = 2;

#elif TEST == 6
  double marray[] = {
    1,-1,-1,-1,
   -1, 1,-1,-1,
    1, 1, 2, 0,
    1, 1, 0, 2 };
  double qarray[] = { 3, 5,-9,-5 };
  double ans[] = { 2, 1, 3, 1 };
  int d = 4;

#elif TEST == 7
  /* unable-to-solve case */
  double marray[] = {
    -1, 0,-3,
     1,-2,-5,
    -2,-1,-2 };
  double qarray[] = { -3,-2,-1 };
  double ans[] = { 0, 0, 0 };
  int d = 3;

#elif TEST == 8
  double marray[] = {
    2, 0,
    0, 2 };
  double qarray[] = { 1, -1 };
  double ans[] = { 0, 0.5 };
  int d = 2;

#else
  double marray[] = {
    1, 0,-1,
    0, 1,-1,
    1, 1, 0 };
  double qarray[] = { 0, 0,-1 };
  double ans[] = { 0.5, 0.5, 0.5 };
  int d = 3;

#endif

  z = zVecAlloc( d );
  m = zMatCloneArray( marray, d, d );
  q = zVecCloneArray( qarray, d );
  zMatWrite( m );
  zVecWrite( q );
  printf( "true ans =\n%d ( ", sizeof(ans)/sizeof(double) );
  zRawVecWrite( ans, sizeof(ans)/sizeof(double) );
  printf( "\n*** result ***\n" );
  zLCPSolveIP( m, q, NULL, z );
  zVecWrite( z );
  zMatFree( m );
  zVecFree( q );
  return 0;
}
