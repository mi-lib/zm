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
  return false; /* dummy */
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
  zVec min, max, start, goal;
  FILE *fp;

  zRandInit();
  min = zVecCreateList( 2,-10.0,-10.0 );
  max = zVecCreateList( 2, 10.0, 10.0 );
  start = zVecCreateList( 2, 9.0, -9.0 );
  goal = zVecCreateList( 2,-9.0, 9.0 );
  zRRTInit( &rrt, min, max, 1e-1, NULL, NULL, testchk );
  zRRTConnect( &rrt, start, goal, 0, NULL, &path );

  fp = fopen( "a", "w" );
  zListForEach( &rrt.slist, rc ){
    if( !rc->data.parent ) continue;
    zVecDataFWrite( fp, rc->data.parent->v );
    zVecDataFWrite( fp, rc->data.v );
    fprintf( fp, "\n" );
  }
  fclose( fp );

  fp = fopen( "b", "w" );
  zListForEach( &rrt.glist, rc ){
    if( !rc->data.parent ) continue;
    zVecDataFWrite( fp, rc->data.parent->v );
    zVecDataFWrite( fp, rc->data.v );
    fprintf( fp, "\n" );
  }
  fclose( fp );

  fp = fopen( "c", "w" );
  zVecListFWrite( fp, &path );
  fclose( fp );

  zRRTPathShortcut( &rrt, NULL, &path );
  fp = fopen( "d", "w" );
  zVecListFWrite( fp, &path );
  fclose( fp );

  fp = fopen( "e", "w" );
  output_wall( fp );
  fclose( fp );

  zRRTDestroy( &rrt );
  zRRTPathDestroy( &path );
  return 0;
}
