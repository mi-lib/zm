#define DEBUG

#include <zm/zm_opt.h>

#define TEST 8

int main(void)
{
  zMat q, a;
  zVec c, b, x;
  double cost;

#if TEST == 1
  double qarray[] = { 4, -1, -1, 4 };
  double aarray[] = { -1, -2, -2, -1, -1, -1 };
  double carray[] = { -1, -2 };
  double barray[] = { -1, -1, -10 };
  double ans[] = { 1.0/4, 3.0/8 };
  int d = 2, dc = 3;

#elif TEST == 2
  /* in this case, answer is [ 2/3  4/3 ] */
  double qarray[] = { 1, -1, -1, 2 };
  double aarray[] = { -1, -1, 1, -2, -2, -1 };
  double carray[] = { -2, -6 };
  double barray[] = { -2, -2, -3 };
  double ans[] = { 2.0/3, 4.0/3 };
  int d = 2, dc = 3;

#elif TEST == 3
  double qarray[] = {
   -2, 0,
    0,-2 };
  double carray[] = {
   -1,
   -1 };
  double aarray[] = {
   -2.5,-5.0,
   -5.0,-6.0,
   -3.0,-2.0, };
  double barray[] = {
   -350,
   -450,
   -240 };
  double ans[] = { 80, 0 };
  int d = 2, dc = 3;

#elif TEST == 4
  double qarray[] = { 1, 0.5, 0.5, 1 };
  double aarray[] = { -1, -2, -1, -1 };
  double carray[] = { -10, -11 };
  double barray[] = { -7, -5 };
  double ans[] = { 3, 2 };
  int d = 2, dc = 2;

#elif TEST == 5
  double qarray[] = {
    1 };
  double aarray[] = {
    1 };
  double carray[] = { -2 };
  double barray[] = { -20 };
  double ans[] = { 2 };
  int d = 1, dc = 1;

#elif TEST == 6
  double qarray[] = {
    1,-1,
   -1, 1 };
  double aarray[] = {
    1,-1 };
  double carray[] = { 2, -2 };
  double barray[] = { -20 };
  double ans[] = { 0, 2 };
  int d = 2, dc = 1;

#elif TEST == 7
  double qarray[] = {
    6,-2, 1,
   -2, 6, 5,
    1, 5, 8 };
  double aarray[] = {
    4,-2, 1,
    1, 0, 5,
    2, 4,-2 };
  double carray[] = { 0, -2, 4 };
  double barray[] = { 15, 17, 8 };
  double ans[] = { 3.8, 1.42, 2.64 };
  int d = 3, dc = 3;

#elif TEST == 8
  double qarray[] = { 4, 1, 1, 2 };
  double carray[] = { -3, -4 };
  double aarray[] = { 1, 2, -1, -2 };
  double barray[] = { 1, -1 };
  double ans[] = { 2.0/7.0, 5.0/14.0 };
  int d = 2, dc = 2;

#else
  /* in this case, answer is [1/2 1/2] */
  double qarray[] = {
    1, 0,
    0, 1 };
  double carray[] = {
    0, 0 };
  double aarray[] = {
    1, 1 };
  double barray[] = {
    1 };
  double ans[] = { 0.5, 0.5 };
  int d = 2, dc = 1;

#endif

  x = zVecAlloc( d );
  q = zMatCloneArray( qarray, d, d );
  c = zVecCloneArray( carray, d );
  a = zMatCloneArray( aarray, dc, d );
  b = zVecCloneArray( barray, dc );
  zMatWrite( q );
  zVecWrite( c );
  zMatWrite( a );
  zVecWrite( b );
  printf( "true ans =\n%d ( ", sizeof(ans)/sizeof(double) );
  zRawVecWrite( ans, sizeof(ans)/sizeof(double) );
  printf( "\n*** result ***\n" );
  zQPSolveLemke( q, c, a, b, x, &cost );
  printf( "\nresult = %f\n", cost );
  zVecWrite( x );
  return 0;
}
