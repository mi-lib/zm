#include <zm/zm_vec.h>
#include <sys/time.h>

int deltatime(struct timeval *tv1, struct timeval *tv2)
{
  return (int)( (tv2->tv_sec-tv1->tv_sec)*1000000+tv2->tv_usec-tv1->tv_usec );
}

#define N   1000
#define DIM 3

int main(int argc, char *argv[])
{
  zVecList list; /* for comparison */
  zVecTree tree, *node;
  zVec v, nn;
  int i, n, dim;
  struct timeval tv1, tv2;
  double dmin1, dmin2;

  n = argc > 1 ? atoi( argv[1] ) : N;
  dim = argc > 2 ? atoi( argv[2] ) : DIM;
  zRandInit();

  zListInit( &list );
  zVecTreeInit( &tree, dim );
  v = zVecAlloc( dim );
  for( i=0; i<n; i++ ){
    zVecRandUniform( v, -10, 10 );
    zVecListInsertHead( &list, v, true );
    zVecTreeAdd( &tree, v );
  }
  zVecRandUniform( v, -10, 10 );

  gettimeofday( &tv1, NULL );
  dmin1 = zVecTreeNN( &tree, v, &node );
  gettimeofday( &tv2, NULL );
  eprintf( "kd-tree: %g - ", dmin1 ); zVecFWrite( stderr, node->v );
  printf( "%d T=%d ", n, deltatime(&tv1,&tv2) );

  /* for comparison */
  gettimeofday( &tv1, NULL );
  nn = zVecListNN( &list, v, &dmin2 );
  gettimeofday( &tv2, NULL );
  eprintf( "naive  : %g - ", dmin2 ); zVecFWrite( stderr, nn );
  printf( "Tn=%d\n", deltatime(&tv1,&tv2) );

  zVecListDestroy( &list, true );
  zVecTreeDestroy( &tree );
  return 0;
}
