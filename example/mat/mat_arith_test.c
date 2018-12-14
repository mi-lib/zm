#include <zm/zm_mat.h>

#define FILENAME "matrix.dat"

int main(void)
{
  zMat m1, m2, m;
  FILE *fp;

  if( !( fp = fopen( FILENAME, "r" ) ) ){
    ZOPENERROR( FILENAME );
    exit( 1 );
  }
  m1 = zMatFRead( fp );
  m2 = zMatFRead( fp );
  printf( "m1: " ); zMatWrite( m1 );
  printf( "m2: " ); zMatWrite( m2 );
  fclose( fp );

  m = zMatAlloc( zMatRowSize(m1), zMatColSize(m2) );
  printf( "m1 + m2: " ); zMatWrite( zMatAdd( m1, m2, m ) );
  printf( "m1 - m2: " ); zMatWrite( zMatSub( m1, m2, m ) );
  printf( "-m1    : " ); zMatWrite( zMatRev( m1, m ) );
  printf( "m1 * 2 : " ); zMatWrite( zMatMul( m1, 2, m ) );
  printf( "m1 / 2 : " ); zMatWrite( zMatDiv( m1, 2, m ) );

  zMatFree( m1 );
  zMatFree( m2 );
  zMatFree( m );
  return 0;
}
