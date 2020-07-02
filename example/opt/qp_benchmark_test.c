#include <zm/zm_opt.h>

void alloc_sample(int d, double qarray[], double carray[], int dc, double aarray[], double barray[], double ansarray[], zMat *q, zVec *c, zMat *a, zVec *b, zVec *ans, zVec *x)
{
  *q = zMatCloneArray( qarray, d, d );
  *c = zVecCloneArray( carray, d );
  *a = zMatCloneArray( aarray, dc, d );
  *b = zVecCloneArray( barray, dc );
  *ans = zVecCloneArray( ansarray, d );
  *x = zVecAlloc( d );
}

void gen_sample1(zMat *q, zVec *c, zMat *a, zVec *b, zVec *ans, zVec *x)
{
  double qarray[] = { 4, -1, -1, 4 };
  double aarray[] = { -1, -2, -2, -1, -1, -1 };
  double carray[] = { -1, -2 };
  double barray[] = { -1, -1, -10 };
  double ansarray[] = { 1.0/4, 3.0/8 };
  alloc_sample( 2, qarray, carray, 3, aarray, barray, ansarray, q, c, a, b, ans, x );
}

void gen_sample2(zMat *q, zVec *c, zMat *a, zVec *b, zVec *ans, zVec *x)
{
  /* in this case, answer is [ 2/3  4/3 ] */
  double qarray[] = { 1, -1, -1, 2 };
  double aarray[] = { -1, -1, 1, -2, -2, -1 };
  double carray[] = { -2, -6 };
  double barray[] = { -2, -2, -3 };
  double ansarray[] = { 2.0/3, 4.0/3 };
  alloc_sample( 2, qarray, carray, 3, aarray, barray, ansarray, q, c, a, b, ans, x );
}

void gen_sample3(zMat *q, zVec *c, zMat *a, zVec *b, zVec *ans, zVec *x)
{
  double qarray[] = { -2, 0, 0,-2 };
  double carray[] = { -1,-1 };
  double aarray[] = {
   -2.5,-5.0,
   -5.0,-6.0,
   -3.0,-2.0, };
  double barray[] = { -350, -450, -240 };
  double ansarray[] = { 80, 0 };
  alloc_sample( 2, qarray, carray, 3, aarray, barray, ansarray, q, c, a, b, ans, x );
}

void gen_sample4(zMat *q, zVec *c, zMat *a, zVec *b, zVec *ans, zVec *x)
{
  double qarray[] = { 1, 0.5, 0.5, 1 };
  double aarray[] = { -1, -2, -1, -1 };
  double carray[] = { -10, -11 };
  double barray[] = { -7, -5 };
  double ansarray[] = { 3, 2 };
  alloc_sample( 2, qarray, carray, 2, aarray, barray, ansarray, q, c, a, b, ans, x );
}

void gen_sample5(zMat *q, zVec *c, zMat *a, zVec *b, zVec *ans, zVec *x)
{
  double qarray[] = { 1 };
  double aarray[] = { 1 };
  double carray[] = { -2 };
  double barray[] = { -20 };
  double ansarray[] = { 2 };
  alloc_sample( 1, qarray, carray, 1, aarray, barray, ansarray, q, c, a, b, ans, x );
}

void gen_sample6(zMat *q, zVec *c, zMat *a, zVec *b, zVec *ans, zVec *x)
{
  double qarray[] = { 1,-1,-1, 1 };
  double aarray[] = { 1,-1 };
  double carray[] = { 2, -2 };
  double barray[] = { -20 };
  double ansarray[] = { 0, 2 };
  alloc_sample( 2, qarray, carray, 1, aarray, barray, ansarray, q, c, a, b, ans, x );
}

void gen_sample7(zMat *q, zVec *c, zMat *a, zVec *b, zVec *ans, zVec *x)
{
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
  double ansarray[] = { 3.8, 1.42, 2.64 };
  alloc_sample( 3, qarray, carray, 3, aarray, barray, ansarray, q, c, a, b, ans, x );
}

void gen_sample8(zMat *q, zVec *c, zMat *a, zVec *b, zVec *ans, zVec *x)
{
  double qarray[] = { 4, 1, 1, 2 };
  double carray[] = { -3, -4 };
  double aarray[] = { 1, 2, -1, -2 };
  double barray[] = { 1, -1 };
  double ansarray[] = { 2.0/7.0, 5.0/14.0 };
  alloc_sample( 2, qarray, carray, 2, aarray, barray, ansarray, q, c, a, b, ans, x );
}

void gen_sample9(zMat *q, zVec *c, zMat *a, zVec *b, zVec *ans, zVec *x)
{
  /* in this case, answer is [1/2 1/2] */
  double qarray[] = { 1, 0, 0, 1 };
  double carray[] = { 0, 0 };
  double aarray[] = { 1, 1 };
  double barray[] = { 1 };
  double ansarray[] = { 0.5, 0.5 };
  alloc_sample( 2, qarray, carray, 1, aarray, barray, ansarray, q, c, a, b, ans, x );
}


void (* gen_sample[])(zMat*,zVec*,zMat*,zVec*,zVec*,zVec*) = {
  gen_sample1,
  gen_sample2,
  gen_sample3,
  gen_sample4,
  gen_sample5,
  gen_sample6,
  gen_sample7,
  gen_sample8,
  gen_sample9,
  NULL,
};

#define TOL (1.0e-7)

int main(void)
{
  zMat q, a;
  zVec c, b, x, ans;
  double cost, err;
  int i;

  for( i=0; gen_sample[i]; i++ ){
    gen_sample[i]( &q, &c, &a, &b, &ans, &x );
    printf( "[sample #%d] (cond.=%d, var.=%d)\t", i+1, zMatRowSize(a), zMatRowSize(q) );
    printf( "Lemke ... " );
    if( !zQPSolveLemke( q, c, a, b, x, &cost ) ){
      printf( "unsolvable\t" );
    } else{
      if( zIsTol( ( err = cost - zQuadraticValue( q, c, x ) ), TOL ) )
        printf( "OK\t" );
      else
        printf( "failed (%.10g)\t", err );
    }
    printf( "IP ... " );
    if( !zQPSolveIP( q, c, a, b, x, &cost ) ){
      printf( "unsolvable\n" );
    } else{
      if( zIsTol( ( err = cost - zQuadraticValue( q, c, x ) ), TOL ) )
        printf( "OK\n" );
      else
        printf( "failed (%.10g)\n", err );
    }
    zMatFreeAO( 2, q, a );
    zVecFreeAO( 4, c, b, ans, x );
  }
  return 0;
}
