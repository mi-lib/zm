#include <zm/zm_misc.h>

#define N1 100
#define N2  20

int main(int argc, char *argv[])
{
  int n1, n2;

  n1 = argc > 1 ? atof(argv[1]) : N1;
  n2 = argc > 2 ? atof(argv[2]) : N2;
  printf( "%d_P_%d = %f\n", n1, n2, zPermutation( n1, n2 ) );
  printf( "%d! = %f\n", n1, zFactorial( n1 ) );
  printf( "%d_C_%d = %f\n", n1, n2, zCombination( n1, n2 ) );
  return 0;
}
