#include <zm/zm_opt.h>

#define TOL (1.0e-7)

void alloc_sample(int d, double marray[], double qarray[], double ansarray[], zMat *m, zVec *q, zVec *ans, zVec *z, zVec *err)
{
  *m = zMatCloneArray( marray, d, d );
  *q = zVecCloneArray( qarray, d );
  *ans = zVecCloneArray( ansarray, d );
  *z = zVecAlloc( d );
  *err = zVecAlloc( d );
}

void gen_sample1(zMat *m, zVec *q, zVec *ans, zVec *z, zVec *err)
{
  double marray[] = {
    4,-1, 1, 2, 1,
   -1, 4, 2, 1, 1,
   -1,-2, 0, 0, 0,
   -2,-1, 0, 0, 0,
   -1,-1, 0, 0, 0 };
  double qarray[] = { -1, -2, 1, 1, 10 };
  double ansarray[] = { 1.0/4.0, 3.0/8.0, 3.0/8.0, 0, 0 };
  alloc_sample( 5, marray, qarray, ansarray, m, q, ans, z, err );
}

void gen_sample2(zMat *m, zVec *q, zVec *ans, zVec *z, zVec *err)
{
  double marray[] = {
    1,-1, 1,-1, 2,
   -1, 2, 1, 2, 1,
   -1,-1, 0, 0, 0,
    1,-2, 0, 0, 0,
   -2,-1, 0, 0, 0 };
  double qarray[] = { -2, -6, 2, 2, 3 };
  double ansarray[] = { 2.0/3.0, 4.0/3.0, 28.0/9, 4.0/9, 0 };
  alloc_sample( 5, marray, qarray, ansarray, m, q, ans, z, err );
}

void gen_sample3(zMat *m, zVec *q, zVec *ans, zVec *z, zVec *err)
{
  double marray[] = {
    -2.0, 0.0, 2.5, 5.0, 3.0,
     0.0,-2.0, 5.0, 6.0, 2.0,
    -2.5,-5.0, 0.0, 0.0, 0.0,
    -5.0,-6.0, 0.0, 0.0, 0.0,
    -3.0,-2.0, 0.0, 0.0, 0.0 };
  double qarray[] = { -1.0, -1.0, 350, 450, 240 };
  double ansarray[] = { 80.0, 0.0, 0.0, 0.0, 161.0/3 };
  alloc_sample( 5, marray, qarray, ansarray, m, q, ans, z, err );
}

void gen_sample4(zMat *m, zVec *q, zVec *ans, zVec *z, zVec *err)
{
  /* in this case, answer is [ 3  2 ] */
  double marray[] = {
    1, 0.5, 1, 1,
    0.5, 1, 2, 1,
    -1, -2, 0, 0,
    -1, -1, 0, 0 };
  double qarray[] = { -10, -11, 7, 5 };
  double ansarray[] = { 3.0, 2.0, 3.0/2, 9.0/2 };
  alloc_sample( 4, marray, qarray, ansarray, m, q, ans, z, err );
}

void gen_sample5(zMat *m, zVec *q, zVec *ans, zVec *z, zVec *err)
{
  double marray[] = {
    1,-1,
    1, 0 };
  double qarray[] = { -2, 20 };
  double ansarray[] = { 2.0, 0.0 };
  alloc_sample( 2, marray, qarray, ansarray, m, q, ans, z, err );
}

void gen_sample6(zMat *m, zVec *q, zVec *ans, zVec *z, zVec *err)
{
  double marray[] = {
    1,-1,-1,-1,
   -1, 1,-1,-1,
    1, 1, 2, 0,
    1, 1, 0, 2 };
  double qarray[] = { 3, 5,-9,-5 };
  double ansarray[] = { 2, 1, 3, 1 };
  alloc_sample( 4, marray, qarray, ansarray, m, q, ans, z, err );
}

void gen_sample7(zMat *m, zVec *q, zVec *ans, zVec *z, zVec *err)
{
  /* unsolvable case */
  double marray[] = {
    -1, 0,-3,
     1,-2,-5,
    -2,-1,-2 };
  double qarray[] = { -3,-2,-1 };
  double ansarray[] = { 0, 0, 0 };
  alloc_sample( 3, marray, qarray, ansarray, m, q, ans, z, err );
}

void gen_sample8(zMat *m, zVec *q, zVec *ans, zVec *z, zVec *err)
{
  double marray[] = {
    2, 0,
    0, 2 };
  double qarray[] = { 1, -1 };
  double ansarray[] = { 0, 0.5 };
  alloc_sample( 2, marray, qarray, ansarray, m, q, ans, z, err );
}

void gen_sample9(zMat *m, zVec *q, zVec *ans, zVec *z, zVec *err)
{
  double marray[] = {
    1, 0,-1,
    0, 1,-1,
    1, 1, 0 };
  double qarray[] = { 0, 0,-1 };
  double ansarray[] = { 0.5, 0.5, 0.5 };
  alloc_sample( 3, marray, qarray, ansarray, m, q, ans, z, err );
}


void (* gen_sample[])(zMat*,zVec*,zVec*,zVec*,zVec*) = {
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

int main(void)
{
  zMat m;
  zVec q, z, ans, err;
  int i;

  for( i=0; gen_sample[i]; i++ ){
    gen_sample[i]( &m, &q, &ans, &z, &err );
    printf( "[sample #%d] (var.=%d)\t", i+1, zVecSize(z) );
    printf( "Lemke ... " );
    if( !zLCPSolveLemke( m, q, NULL, z ) ){
      printf( "unsolvable\t" );
    } else{
      zVecSub( ans, z, err );
      if( zVecIsTol( err, TOL ) )
        printf( "OK\t" );
      else
        printf( "failed (%.10g)\t", zVecNorm(err) );
    }
    printf( "IP ... " );
    if( !zLCPSolveIP( m, q, NULL, z ) ){
      printf( "unsolvable\n" );
    } else{
      zVecSub( ans, z, err );
      if( zVecIsTol( err, TOL ) )
        printf( "OK\n" );
      else
        printf( "failed (%.10g)\n", zVecNorm(err) );
    }
    zMatFree( m );
    zVecFreeAO( 4, q, ans, z, err );
  }
  return 0;
}
