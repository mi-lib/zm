#include <zm/zm_opt.h>

#define TOL (1.0e-7)

void alloc_sample(int row, int col, double aa[], double ba[], double ca[], double ansa[], zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  *a = zMatCloneArray( aa, row, col );
  *b = zVecCloneArray( ba, row );
  *c = zVecCloneArray( ca, col );
  *ans = zVecCloneArray( ansa, col );
  *x = zVecAlloc( col );
}

void gen_sample1(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
    3.0, 1.0, 1.0, 0.0, 0.0,
    2.5, 2.0, 0.0, 1.0, 0.0,
    1.0, 2.0, 0.0, 0.0, 1.0 };
  double b_arr[] = {
    9.0, 12.5, 8.0 };
  double c_arr[] = {
   -3.0, -2.0, 0.0, 0.0, 0.0 };
  double ans_arr[] = {
    2.0, 3.0, 0.0, 0.0, 0.0 };
  alloc_sample( 3, 5, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

void gen_sample2(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
    1.0, 2.0, 1.0, 0.0, 0.0,
    3.0, 4.0, 0.0, 1.0, 0.0,
    3.0, 1.0, 0.0, 0.0, 1.0,
  };
  double b_arr[] = {
    8.0, 18.0, 15.0,
  };
  double c_arr[] = {
    -2.0, -3.0, 0.0, 0.0, 0.0,
  };
  double ans_arr[] = { 2.0, 3.0 };
  alloc_sample( 3, 5, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

void gen_sample3(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
   -1.0, 2.0,-2.0, 1.0, 0.0,
    2.0, 3.0, 1.0, 0.0, 1.0,
  };
  double b_arr[] = {
   -8.0, 5.0,
  };
  double c_arr[] = {
    5.0, 3.0,-1.0, 0.0, 0.0,
  };
  double ans_arr[] = { 0.0, 0.0, 5.0 };
  alloc_sample( 2, 5, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

void gen_sample4(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
    1.0,-1.0, 1.0, 1.0, 0.0, 0.0,
    3.0, 2.0, 4.0, 0.0, 1.0, 0.0,
    3.0, 2.0, 0.0, 0.0, 0.0, 1.0,
  };
  double b_arr[] = {
    20.0, 42.0, 30.0,
  };
  double c_arr[] = {
    -5.0,-4.0,-6.0, 0.0, 0.0, 0.0,
  };
  double ans_arr[] = { 0.0, 15.0, 3.0 };
  alloc_sample( 3, 6, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

void gen_sample5(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
    1.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 2.0, 0.0,-7.0, 0.0, 1.0, 0.0,
    0.0,-1.0, 1.0,-2.0, 0.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
  };
  double b_arr[] = {
    740.0, 0.0, -0.5, 9.0,
  };
  double c_arr[] = {
   -1.0,-1.0,-3.0, 0.5, 0.0, 0.0, 0.0,
  };
  double ans_arr[] = { 0.0, 3.325, 4.725, 0.95 };
  alloc_sample( 4, 7, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

void gen_sample6(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
    6.0,-1.0, 1.0, 0.0,
   -3.0, 4.0, 0.0, 1.0,
  };
  double b_arr[] = {
    2.0, 8.0,
  };
  double c_arr[] = {
   -2.0, 4.0, 0.0, 0.0,
  };
  double ans_arr[] = { 1.0/3.0, 0.0, 0.0, 9.0 };
  alloc_sample( 2, 4, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

void gen_sample7(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
    1.0, 4.0, 3.0, 1.0, 0.0,
   -1.0, 2.0,-3.0, 0.0, 1.0,
  };
  double b_arr[] = {
    12.0, 4.0,
  };
  double c_arr[] = {
    1.0, -3.0, 1.0, 0.0, 0.0,
  };
  double ans_arr[] = { 0.0, 8.0/3.0, 4.0/9.0 };
  alloc_sample( 2, 5, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

void gen_sample8(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
   -2.5,-3.0,-5.0, 1.0, 0.0, 0.0,
   -2.5,-2.0,-3.0, 0.0, 1.0, 0.0,
   -3.0,-1.0,-2.0, 0.0, 0.0, 1.0,
  };
  double b_arr[] = {
   -200.0,-160.0,-120.0,
  };
  double c_arr[] = {
    9.0, 5.0, 8.0, 0.0, 0.0, 0.0,
  };
  double ans_arr[] = { 10.0, 0.0, 45.0 };
  alloc_sample( 3, 6, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

void gen_sample9(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
    2.0, 1.0, 4.0, 5.0,
    1.0, 2.0, 4.0, 3.0,
  };
  double b_arr[] = {
    34.0, 22.0,
  };
  double c_arr[] = {
    1.0, 1.0, 3.0, 3.0,
  };
  double ans_arr[] = { 46.0/3.0, 10.0/3.0, 0.0, 0.0 };
  alloc_sample( 2, 4, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

void gen_sample10(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
    6.0, 4.0, 1.0, 0.0, 0.0,
    2.0, 3.0, 0.0, 1.0, 0.0,
    2.0, 1.0, 0.0, 0.0, 1.0,
  };
  double b_arr[] = {
    38.0, 21.0, 12.0,
  };
  double c_arr[] = {
   -4.0,-3.0, 0.0, 0.0, 0.0,
  };
  double ans_arr[] = { 3.0, 5.0 };
  alloc_sample( 3, 5, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

void gen_sample11(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
    0.5,-5.5,-2.5, 9.0, 1.0, 0.0, 0.0,
    0.5,-1.5,-0.5, 1.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
  };
  double b_arr[] = {
    0.0, 0.0, 1.0,
  };
  double c_arr[] = {
   -10.0, 57.0, 9.0, 24.0, 0.0, 0.0, 0.0,
  };
  double ans_arr[] = { 1.0, 0.0, 1.0, 0.0, 2.0, 0.0, 0.0 };
  alloc_sample( 3, 7, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

void gen_sample12(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
    5.0, 0.0, 6.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 2.0, 8.0, 0.0, 1.0, 0.0, 0.0,
    7.0, 0.0,15.0, 0.0, 0.0, 1.0, 0.0,
    3.0,11.0, 0.0, 0.0, 0.0, 0.0, 1.0,
  };
  double b_arr[] = {
    8.0, 5.0, 10.0, 7.0,
  };
  double c_arr[] = {
    7.0, 12.0, 3.0, 0.0, 0.0, 0.0, 0.0,
  };
  double ans_arr[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  alloc_sample( 3, 7, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

void gen_sample13(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
    5.0, 1.0, 2.0, 1.0, 0.0, 0.0,
    2.0, 2.0, 6.0, 0.0, 1.0, 0.0,
    2.0, 6.0, 4.0, 0.0, 0.0, 1.0,
  };
  double b_arr[] = {
    20.0, 30.0, 40.0,
  };
  double c_arr[] = {
   -7.0,-5.0,-4.0, 0.0, 0.0, 0.0,
  };
  double ans_arr[] = { 20.0/7.0, 40.0/7.0, 0.0, 0.0, 0.0, 0.0 };
  alloc_sample( 3, 6, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

void gen_sample14(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
    1.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 2.0, 0.0,-7.0, 0.0, 1.0, 0.0,
    0.0,-1.0, 1.0,-2.0, 0.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
  };
  double b_arr[] = {
    740.0, 0.0, -0.5, 9.0,
  };
  double c_arr[] = {
   -1.0, -1.0, -3.0, 0.5, 0.0, 0.0, 0.0,
  };
  double ans_arr[] = { 0.0, 3.325, 4.725, 0.95, 0.0, 0.0, 0.0 };
  alloc_sample( 4, 7, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

void gen_sample15(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
          1.0,       1.0,       1.0,        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    10.0/24.0,  5.0/12.0, 15.0/14.0, 120.0/18.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    10.0/24.0, 25.0/12.0, 35.0/14.0,  20.0/18.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
         -1.0,       0.0,       0.0,        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
          0.0,      -1.0,       0.0,        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
          0.0,       0.0,      -1.0,        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
          0.0,       0.0,       0.0,       -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
  };
  double b_arr[] = {
    10.0, 24.0, 12.0,-1.0,-1.0,-1.0,-1.0
  };
  double c_arr[] = {
    25.1, 25.4, 20.4, 17.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  };
  double ans_arr[] = { 4.0621052631, 1.0, 1.971052632, 2.966842105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  alloc_sample( 7, 10, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

#if 0
void gen_sample16(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
    -1.9,  -0.4,  -0.6,  -0.7,  -0.6,  -2.9,  -0.6,  1.0, 0.0, 0.0,
    -0.08, -0.02, -0.03, -0.05, -0.04, -0.04, -0.05, 0.0, 1.0, 0.0,
   -16.0, -11.0,  -5.0, -15.0,  -4.0,  -1.0,  -5.0,  0.0, 0.0, 1.0,
  };
  double b_arr[] = {
    -10.0, -0.5, -50.0,
  };
  double c_arr[] = {
    136.5, 157.5, 273.0, 136.5, 157.5, 50.0, 147.0, 0.0, 0.0, 0.0,
  };
  double ans_arr[] = { 2.679, 0.0, 0.0, 0.0, 0.0, 7.143, 0.0, 0.0, 0.0, 0.0 };
  alloc_sample( 3, 10, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

void gen_sample17(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
    4.5, 2.0, 1.5, 1.2, 8.0,-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 3.0, 0.5, 1.0, 2.0, 0.0,-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.5, 0.0, 4.0, 0.0, 0.0,-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 4.0, 8.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
  };
  double b_arr[] = {
    0.0, 0.0, 0.0, 1500.0, 150.0, 100.0, 200.0,
  };
  double c_arr[] = {
   -45.0, -60.0, -15.0, -20.0, -150.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  };
  double ans_arr[] = {
    18.8, 20.0, 0.0, 21.2, 0.0, 150.0, 100.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  };
  alloc_sample( 7, 15, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}
#endif

void gen_sample18(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.4, 0.0, 1.0, 0.4, 0.0, 0.0, 6.4,
    0.9, 0.0, 0.0,-0.35,1.0, 0.0, 6.4,
    0.4, 0.0, 0.0,-0.1, 0.0, 1.0, 2.4,
  };
  double b_arr[] = {
    1.0, 6.4, 6.4, 2.4,
  };
  double c_arr[] = {
   -0.2, 0.0, 0.0, -11.0/80.0, 0.0, 0.0, -6.0/5.0,
  };
  double ans_arr[] = {
    8.0, 1.0, 0.0, 8.0, 2.0, 0.0, 0.0,
  };
  alloc_sample( 4, 7, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

void gen_sample19(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
    1.0, 2.0, 1.0, 0.0, 0.0,
    3.0, 2.0, 0.0, 1.0, 0.0,
    2.0, 1.0, 0.0, 0.0, 1.0,
  };
  double b_arr[] = {
    2.0, 3.0, 2.0,
  };
  double c_arr[] = {
   -1.0, -1.0, 0.0, 0.0, 0.0,
  };
  double ans_arr[] = {
    0.5, 0.75, 0.25, 0.25, 0.0,
  };
  alloc_sample( 3, 5, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
}

void gen_sample20(zMat *a, zVec *b, zVec *c, zVec *ans, zVec *x)
{
  double a_arr[] = {
   -3.0,-5.0,-1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
   -6.0,-2.0,-4.0, 1.0, 0.0, 1.0, 0.0, 0.0,
   -1.0,-4.0,-3.0, 1.0, 0.0, 0.0, 1.0, 0.0,
   -4.0,-2.0,-5.0, 1.0, 0.0, 0.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  };
  double b_arr[] = {
    0.0, 0.0, 0.0, 0.0, 1.0,
  };
  double c_arr[] = {
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
  };
  double ans_arr[] = {
    1.0/8.0, 1.0/2.0, 3.0/8.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  };
  alloc_sample( 5, 8, a_arr, b_arr, c_arr, ans_arr, a, b, c, ans, x );
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
  gen_sample10,
  gen_sample11,
  gen_sample12,
  gen_sample13,
  gen_sample14,
  gen_sample15,
  gen_sample18,
  gen_sample19,
  gen_sample20,
  NULL,
};

int main(int argc, char *argv[])
{
  zMat a;
  zVec b, c, ans, x;
  double cost, err;
  int i;

  for( i=0; gen_sample[i]; i++ ){
    gen_sample[i]( &a, &b, &c, &ans, &x );
    printf( "[sample #%d] (cond.=%d, var.=%d)\t", i+1, zMatRowSize(a), zMatColSize(a) );
    printf( "Simplex ... " );
    if( !zLPSolveSimplex( a, b, c, x, &cost ) ){
      printf( "unsolvable\t" );
    } else{
      if( zIsTol( ( err = cost - zVecInnerProd( ans, c ) ), TOL ) )
        printf( "OK\t" );
      else
        printf( "failed (%.10g)\t", err );
    }
    printf( "PDIP_PC ... " );
    if( !zLPSolvePDIP_PC( a, b, c, x, &cost ) ){
      printf( "unsolvable\n" );
    } else{
      if( zIsTol( ( err = cost - zVecInnerProd( ans, c ) ), TOL ) )
        printf( "OK\n" );
      else
        printf( "failed (%.10g)\n", err );
    }
    zMatFree( a );
    zVecFreeAO( 4, b, c, ans, x );
  }
  return 0;
}
