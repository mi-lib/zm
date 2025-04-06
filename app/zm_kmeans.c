#include <zm/zm.h>

zVecList sample;
zVecMCluster mc;

enum{
  ZM_KMEANS_HELP = 0,
  ZM_KMEANS_VERBOSE,
  ZM_KMEANS_KMIN,
  ZM_KMEANS_KMAX,
  ZM_KMEANS_TRY_N,
  ZM_KMEANS_CLUSTER_FILENAME,
  ZM_KMEANS_SILHOUETTE_FILENAME,
};
zOption option[] = {
  { "h",    "help",    NULL, "display this messages", NULL, false },
  { "v", "verbose",    NULL, "display this messages", NULL, false },
  { "kmin", "kmin",    "<integer value>", "minimum number of clusters to be tested", (char *)"2", false },
  { "kmax", "kmax",    "<integer value>", "maximum number of clusters to be tested", (char *)"10", false },
  { "n",     "try",    "<integer value>", "number of trials of clustering for each k", (char *)"3", false },
  { "c", "cluster",    "<file name>", "common name for files of clusters", (char *)"", false },
  { "s", "silhouette", "<file name>", "common name for silhouettes", (char *)"s", false },
  { NULL, NULL, NULL, NULL, NULL, false },
};

int zm_kmeans_usage(void)
{
  ZECHO( "usage: zm_kmeans <option> [filename]" );
  ZECHO( "<option>" );
  zOptionHelp( option );
  exit( 0 );
}

bool zm_kmeans_init(int argc, char *argv[])
{
  zStrAddrList arglist;
  FILE *fp;

  if( !zOptionRead( option, argv, &arglist ) ) return false;
  if( option[ZM_KMEANS_HELP].flag )
    zm_kmeans_usage();
  if( zListIsEmpty( &arglist ) ){
    zVecListScan( &sample );
  } else{
    if( !( fp = fopen( zListTail(&arglist)->data, "r" ) ) ){
      ZOPENERROR( zListTail(&arglist)->data );
      return false;
    }
    zVecListFScan( fp, &sample );
    fclose( fp );
  }
  zStrAddrListDestroy( &arglist );
  return true;
}

double kmeans_try_k(zVecMCluster *mc, zVecList *sample, int k, int testnum)
{
  zVecMCluster *test_mc;
  double *score = NULL;
  int i, selected;
  double score_min = -HUGE_VAL;

  if( !( test_mc = zAlloc( zVecMCluster, testnum ) ) ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  if( !( score = zAlloc( double, testnum ) ) ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  for( i=0; i<testnum; i++ ){
    zVecMClusterInit( &test_mc[i], zVecSizeNC(zListTail(sample)->data) );
    zVecMClusterKMeans( &test_mc[i], sample, k );
    score[i] = zVecMClusterMeanSilhouette( &test_mc[i] );
  }
  score_min = zDataMax( score, testnum, &selected );
  zVecMClusterMove( &test_mc[selected], mc );

 TERMINATE:
  for( i=0; i<testnum; i++ )
    zVecMClusterDestroy( &test_mc[i] );
  free( test_mc );
  free( score );
  return score_min;
}

bool kmeans_check_silhouette(zVecMCluster *mc, double score)
{
  zVecClusterListCell *cp;

  zListForEach( zVecMClusterClusterList(mc), cp )
    if( zVecClusterMaxSilhouette(&cp->data) < score ) return false;
  return true;
}

int kmeans_try(zVecMCluster *mc, zVecList *sample, int kmin, int kmax, int testnum)
{
  zVecMCluster *test_mc;
  double score;
  double evenness, evenness_min = HUGE_VAL;
  int n, i, ik = -1;

  n = kmax - kmin + 1;
  if( !( test_mc = zAlloc( zVecMCluster, n ) ) ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  for( i=0; i<n; i++ ){
    if( option[ZM_KMEANS_VERBOSE].flag ) printf( "trying k=%d ... ", kmin + i );
    score = kmeans_try_k( &test_mc[i], sample, kmin+i, testnum );
    if( option[ZM_KMEANS_VERBOSE].flag ) printf( "score = %g, ", score );
    if( !kmeans_check_silhouette( &test_mc[i], score ) ){
      if( option[ZM_KMEANS_VERBOSE].flag ) printf( "omitted.\n" );
      continue;
    }
    if( ( evenness = zVecMClusterEvenness( &test_mc[i] ) ) < evenness_min ){
      evenness_min = evenness;
      ik = i;
    }
    if( option[ZM_KMEANS_VERBOSE].flag ) printf( "evenness = %g\n", evenness );
  }
  zVecMClusterMove( &test_mc[ik], mc );

 TERMINATE:
  for( i=0; i<testnum; i++ )
    zVecMClusterDestroy( &test_mc[i] );
  free( test_mc );
  return kmin + ik;
}

int main(int argc, char *argv[])
{
  int k, kmin, kmax, n;

  if( !zm_kmeans_init( argc, argv+1 ) ) return EXIT_FAILURE;
  kmin = atoi( option[ZM_KMEANS_KMIN].arg );
  kmax = atoi( option[ZM_KMEANS_KMAX].arg );
  n    = atoi( option[ZM_KMEANS_TRY_N].arg );
  zRandInit();
  k = kmeans_try( &mc, &sample, kmin, kmax, n );
  if( option[ZM_KMEANS_VERBOSE].flag )
    printf( "determined number of clusters = %d\n", k );
  zVecMClusterValuePrintFile( &mc, option[ZM_KMEANS_CLUSTER_FILENAME].arg );
  zVecMClusterSilhouettePrintFile( &mc, option[ZM_KMEANS_SILHOUETTE_FILENAME].arg );
  zVecMClusterDestroy( &mc );
  zVecListDestroy( &sample );
  return EXIT_SUCCESS;
}
