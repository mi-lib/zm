#include <zm/zm_mat.h>

int main(void)
{
  zMat m;
  double val[] = {
    1, 2, 3, 4,
    5, 6, 7, 8,
    9,10,11,12,
   13,14,15,16,
  };

  m = zMatCloneArray( val, 4, 4 );

  printf( "m: " ); zMatWrite( m );
  printf( "swap row( 0, 2 ): " ); zMatWrite( zMatSwapRow( m, 0, 2 ) );
  printf( "swap row( 1, 0 ): " ); zMatWrite( zMatSwapRow( m, 1, 0 ) );
  printf( "swap row( 3, 1 ): " ); zMatWrite( zMatSwapRow( m, 3, 1 ) );
  printf( "swap row( 4, 0 ): " ); zMatWrite( zMatSwapRow( m, 4, 0 ) );
  printf( "swap col( 0, 2 ): " ); zMatWrite( zMatSwapCol( m, 0, 2 ) );
  printf( "swap col( 3, 1 ): " ); zMatWrite( zMatSwapCol( m, 3, 1 ) );
  printf( "swap col( 0, 1 ): " ); zMatWrite( zMatSwapCol( m, 0, 1 ) );
  printf( "swap col( 4, 0 ): " ); zMatWrite( zMatSwapCol( m, 4, 0 ) );

  zMatFree( m );
  return 0;
}
