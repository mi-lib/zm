#include <zm/zm_rrt.h>

bool testchk(zVec v, void *util)
{
  /* wall 1 */
  if( zVecElem(v,0) < 6 && zVecElem(v,1) > 4 && zVecElem(v,1) < 5 )
    return true;
  if( zVecElem(v,0) > 4 && zVecElem(v,0) < 5 && zVecElem(v,1) > -2 && zVecElem(v,1) < 5 )
    return true;
  /* wall 3 */
  if( zVecElem(v,0) >-6 && zVecElem(v,1) >-6 && zVecElem(v,1) <-4 )
    return true;
  if( zVecElem(v,0) >-5 && zVecElem(v,0) <-4 && zVecElem(v,1) >-6 && zVecElem(v,1) < 2 )
    return true;
  return false;
}

bool testgoal(zVec v, void *util)
{
  return zVecElem(v,0) > -10 && zVecElem(v,0) < -9 &&
         zVecElem(v,1) >   9 && zVecElem(v,1) < 10 ? true : false;
}

void output_wall(FILE *fp)
{
  /* wall 1 */
  fprintf( fp, "-10 5\n" );
  fprintf( fp, "6 5\n" );
  fprintf( fp, "6 4\n" );
  fprintf( fp, "5 4\n" );
  fprintf( fp, "5 -2\n" );
  fprintf( fp, "4 -2\n" );
  fprintf( fp, "4 4\n" );
  fprintf( fp, "-10 4\n\n" );
  /* wall 2 */
  fprintf( fp, "10 -6\n" );
  fprintf( fp, "-6 -6\n" );
  fprintf( fp, "-6 -4\n" );
  fprintf( fp, "-5 -4\n" );
  fprintf( fp, "-5  2\n" );
  fprintf( fp, "-4  2\n" );
  fprintf( fp, "-4 -4\n" );
  fprintf( fp, "10 -4\n" );
}

int main(int argc, char *argv[])
{
  zRRT rrt;
  zRRTListCell *rc;
  zVecList path;
  zVec min, max, start;
  double cost;
  FILE *fp;

  zRandInit();
  min = zVecCreateList( 2,-10.0,-10.0 );
  max = zVecCreateList( 2, 10.0, 10.0 );
  start = zVecCreateList( 2, 9.0, -9.0 );
  zRRTInit( &rrt, min, max, 0.5, NULL, NULL, testchk, testgoal );
  zRRTFindPathOpt( &rrt, start, 0, NULL, &path, &cost );
  printf( "cost: %g\n", cost );

  fp = fopen( "a", "w" );
  zListForEach( &rrt.slist, rc ){
    if( !rc->data.parent ) continue;
    zVecDataFPrint( fp, rc->data.parent->v );
    zVecDataFPrint( fp, rc->data.v );
    fprintf( fp, "\n" );
  }
  fclose( fp );

  fp = fopen( "c", "w" );
  zVecListFPrint( fp, &path );
  fclose( fp );

  fp = fopen( "e", "w" );
  output_wall( fp );
  fclose( fp );

  zRRTDestroy( &rrt );
  zVecListDestroy( &path );
  zVecFreeAtOnce( 3, min, max, start );
  return 0;
}
