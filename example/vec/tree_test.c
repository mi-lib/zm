#include <zm/zm_vec.h>
#ifdef __WINDOWS__
#include <time.h>
#else
#include <sys/time.h>
#endif /* __ WINDOWS__ */

unsigned long deltatime(struct timeval *tv1, struct timeval *tv2)
{
  return (unsigned long)( (tv2->tv_sec-tv1->tv_sec)*1000000+tv2->tv_usec-tv1->tv_usec );
}

#ifdef __WINDOWS__
/* FILETIME of 1stJan1970 */
static const unsigned __int64 Epoch = ((unsigned __int64)116444736000000000i64);
int gettimeofday(struct timeval* tp, struct timezone* dummy)
{
    FILETIME file_time;
    SYSTEMTIME system_time;
    ULARGE_INTEGER ularge;
    GetSystemTime(&system_time);
    SystemTimeToFileTime(&system_time, &file_time);
    ularge.LowPart = file_time.dwLowDateTime;
    ularge.HighPart = file_time.dwHighDateTime;
    tp->tv_sec = (long)((ularge.QuadPart - Epoch) / 10000000i64);
    tp->tv_usec = (long)(system_time.wMilliseconds * 1000);
    return 0;
}
#endif

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
    zVecListInsertHead( &list, v );
    zVecTreeAdd( &tree, v );
  }
  zVecRandUniform( v, -10, 10 );

  gettimeofday( &tv1, NULL );
  dmin1 = zVecTreeNN( &tree, v, &node );
  gettimeofday( &tv2, NULL );
  eprintf( "kd-tree: %g - ", dmin1 ); zVecFPrint( stderr, node->v );
  printf( "%d T=%d ", n, deltatime(&tv1,&tv2) );

  /* for comparison */
  gettimeofday( &tv1, NULL );
  nn = zVecListNN( &list, v, &dmin2 );
  gettimeofday( &tv2, NULL );
  eprintf( "naive  : %g - ", dmin2 ); zVecFPrint( stderr, nn );
  printf( "Tn=%d\n", deltatime(&tv1,&tv2) );

  zVecListDestroy( &list );
  return 0;
}
