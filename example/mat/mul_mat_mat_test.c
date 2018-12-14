#include <sys/time.h>
#include <zm/zm_mat.h>

#define N1 500
#define N2 600
#define N3 400

int main(void)
{
  zMat m1, m2, m;
  register int i, j;
  struct timeval tv1, tv2;

  zRandInit();
  m1 = zMatAlloc( N1, N2 );
  m2 = zMatAlloc( N2, N3 );
  m  = zMatAlloc( N1, N3 );
  for( i=0; i<N1; i++ )
    for( j=0; j<N2; j++ )
      zMatSetElem( m1, i, j, zRandF(-100,100) );
  for( i=0; i<N2; i++ )
    for( j=0; j<N3; j++ )
      zMatSetElem( m2, i, j, zRandF(-100,100) );

  gettimeofday( &tv1, NULL );
  zMulMatMat( m1, m2, m );
  gettimeofday( &tv2, NULL ); \
  printf( "%ld\n", tv2.tv_sec*1000000+tv2.tv_usec - tv1.tv_sec*1000000-tv1.tv_usec );

  zMatFree( m1 );
  zMatFree( m2 );
  zMatFree( m );
  return 0;
}
