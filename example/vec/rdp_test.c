#include <zm/zm_vec.h>

#define N 1000

int main(void)
{
  zVecList list;
  zVec v;
  double angle;
  int i;
  FILE *fp;

  zListInit( &list );
  v = zVecAlloc( 2 );
  for( i=0; i<N; i++ ){
    angle = zPIx2 * i / N;
    zVecSetElem( v, 0, cos(angle) );
    zVecSetElem( v, 1, sin(angle) );
    zVecListInsertHead( &list, v );
  }
  fp = fopen( "0", "w" );
  zVecListFPrint( fp, &list );
  fclose( fp );

  zVecListRDP( &list, NULL, NULL, 0.01 );
  fp = fopen( "1", "w" );
  zVecListFPrint( fp, &list );
  fclose( fp );

  zVecListRDP( &list, NULL, NULL, 0.02 );
  fp = fopen( "2", "w" );
  zVecListFPrint( fp, &list );
  fclose( fp );

  zVecListRDP( &list, NULL, NULL, 0.04 );
  fp = fopen( "3", "w" );
  zVecListFPrint( fp, &list );
  fclose( fp );

  zVecListRDP( &list, NULL, NULL, 0.08 );
  fp = fopen( "4", "w" );
  zVecListFPrint( fp, &list );
  fclose( fp );

  zVecListDestroy( &list );
  return 0;
}
