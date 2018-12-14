#include <zm/zm_mat.h>

#define ROW 10
#define COL 10

zMat mat_create_rand(int row, int col)
{
  zMat m;

  m = zMatAlloc( ROW, COL );
  zMatRand( m, -10, 10 );
  return m;
}

zMat mat_create_ut_rand(int row, int col)
{
  register int i, j;
  zMat m;

  m = zMatAlloc( ROW, COL );
  for( i=0; i<row; i++ )
    for( j=i; j<col; j++ )
      zMatSetElem( m, i, j, zRandF(-10,10) );
  return m;
}

zMat mat_create_tri_rand(int row, int col)
{
  register int i, j;
  zMat m;

  m = zMatAlloc( ROW, COL );
  for( i=0; i<row; i++ )
    for( j=zMax(i-1,0); j<i+2 && j<col; j++ )
      zMatSetElem( m, i, j, zRandF(-10,10) );
  return m;
}

zMat mat_create_sym_rand(int row, int col)
{
  register int i, j;
  zMat m;

  m = zMatAlloc( ROW, COL );
  for( i=0; i<row; i++ )
    for( j=i; j<col; j++ )
      zMatSetElem( m, j, i, zMatSetElem( m, i, j, zRandF(-10,10) ) );
  return m;
}

int main(void)
{
  zMat m;

  zRandInit();

  printf( "ordinary matrix\n" );
  m = mat_create_rand( ROW, COL );
  zMatImg( m ); zEndl();
  zMatFree( m );

  printf( "upper triangle matrix\n" );
  m = mat_create_ut_rand( ROW, COL );
  zMatImg( m ); zEndl();
  zMatFree( m );

  printf( "tridiagonal matrix\n" );
  m = mat_create_tri_rand( ROW, COL );
  zMatImg( m ); zEndl();
  zMatFree( m );

  printf( "symmetric matrix\n" );
  m = mat_create_sym_rand( ROW, COL );
  zMatImg( m ); zEndl();
  zMatFree( m );
  return 0;
}
