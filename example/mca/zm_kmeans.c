#include <zm/zm_mca.h>

#define ZM_WARN_VECLIST_EMPTY "empty vector list assigned"

zVec zVecDataFScan(FILE *fp)
{
  char buf[BUFSIZ], valstr[BUFSIZ];
  int size, dim, i;
  char *sp;
  zVec v;

  if( !fgets( buf, BUFSIZ, fp ) ) return NULL;
  for( dim=0, size=strlen(buf), sp=buf; ; dim++ ){
    if( !*( sp = zSTokenSkim( sp, valstr, size ) ) ) break;
    size -= strlen( valstr );
  }
  if( dim == 0 ) return NULL;
  if( !( v = zVecAlloc( dim ) ) ) return NULL;
  for( i=0; i<dim; i++ )
    zSDouble( buf, &zVecElemNC(v,i) );
  return v;
}

zVecList *zVecListFScan(FILE *fp, zVecList *list)
{
  zVec v;

  zListInit( list );
  while( ( v = zVecDataFScan( fp ) ) )
    zVecListInsertHead( list, v );
  if( zListIsEmpty( list ) )
    ZRUNWARN( ZM_WARN_VECLIST_EMPTY );
  return list;
}




bool kmeans_try(zMCluster *mc, zVecList *sample, int k, int testnum)
{
  zMCluster *test_mc;
  double *score;
  int i, selected;
  bool ret = true;

  if( !( test_mc = zAlloc( zMCluster, testnum ) ) ){
    ZALLOCERROR();
    ret = false;
    goto TERMINATE;
  }
  if( !( score = zAlloc( double, testnum ) ) ){
    ZALLOCERROR();
    ret = falase;
    goto TERMINATE;
  }
  for( i=0; i<testnum; i++ ){
    zMClusterInit( &test_mc[i], zVecSizeNC(zListTail(&sample)->data), zVecSizeNC(zListTail(&sample)->data) );
    zMClusterKMeans( &test_mc[i], sample, k );
    score[i] = zMClusterSilhouetteScore( &test_mc[i] );
  }
  zDataMax( score, testnum, &selected );



 TERMINATE:
  for( i=0; i<testnum; i++ )
    zMClusterDestroy( &test_mc[i] );
  return ret;
}

int main(int argc, char *argv[])
{
  zVecList sample;
  zMCluster mc;
  FILE *fp;

  int nc;
  double score;


  if( argc > 1 ){
    if( !( fp = fopen( argv[1], "r" ) ) ){
      ZOPENERROR( argv[1] );
      return EXIT_FAILURE;
    }
  } else
    fp = stdin;
  zVecListFScan( fp, &sample );

  zRandInit();





  zVecListDestroy( &sample );
  if( fp != stdin ) fclose( fp );
  return EXIT_SUCCESS;
}
