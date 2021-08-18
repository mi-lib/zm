#include <zm/zm.h>

int zm_ave_usage(char *argv)
{
  eprintf( "Usage: %s <filename>\n", argv );
  return 0;
}

bool zm_ave_sample(FILE *fp, zVecList *list)
{
  char buf[BUFSIZ];
  double val;
  int i, nd = 0;
  zVec v;

  zListInit( list );
  if( !fgets( buf, BUFSIZ, fp ) ){
    ZRUNERROR( "empty data file" );
    return false;
  }
  while( *zSDouble( buf, &val ) ) nd++;
  rewind( fp );
  while( fgets( buf, BUFSIZ, fp ) ){
    if( !( v = zVecAlloc( nd ) ) ) break;
    for( i=0; i<nd; i++ ){
      zSDouble( buf, &val );
      zVecSetElemNC( v, i, val );
    }
    zVecAddrListInsertHead( list, v );
  }
  return true;
}

bool zm_ave_mean_cov(zVecList *list)
{
  zVec max, min, mean;
  zMat cov;
  bool ret = true;
  int i;

  max = zVecAlloc( zVecSize(zListHead(list)->data) );
  min = zVecAlloc( zVecSize(zListHead(list)->data) );
  mean = zVecAlloc( zVecSize(zListHead(list)->data) );
  cov = zMatAllocSqr( zVecSize(zListHead(list)->data) );
  if( !max || !min || !mean || !cov ){
    ret = false;
    goto TERMINATE;
  }
  zVecListMinMax( list, min, max );
  zVecListMeanCov( list, mean, cov );

  eprintf( "min :" ); zVecDataPrint( min );
  eprintf( "max :" ); zVecDataPrint( max );
  eprintf( "mean:" ); zVecDataPrint( mean );
  eprintf( "var :" );
  for( i=0; i<zMatRowSizeNC(cov); i++ )
    printf( " %.10g", sqrt( zMatElemNC(cov,i,i) ) );
  printf( "\n" );

 TERMINATE:
  zVecFreeAO( 3, max, min, mean );
  zMatFree( cov );
  return ret;
}

int main(int argc, char *argv[])
{
  FILE *fp;
  zVecList list;

  if( argc <= 1 ) return zm_ave_usage( argv[0] );
  if( !( fp = fopen( argv[1], "r" ) ) ){
    ZOPENERROR( argv[1] );
    return 1;
  }
  if( zm_ave_sample( fp, &list ) ){
    eprintf( "%d data x %d sample\n", zVecSize(zListHead(&list)->data), zListSize(&list) );
    zm_ave_mean_cov( &list );
    zVecListDestroy( &list );
  }
  fclose( fp );
  return 0;
}
