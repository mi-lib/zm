#include <zm/zm_rrt.h>

typedef double box_t[4];

void box_create(box_t box, double x1, double y1, double x2, double y2)
{
  box[0] = zMin( x1, x2 );
  box[1] = zMin( y1, y2 );
  box[2] = zMax( x1, x2 );
  box[3] = zMax( y1, y2 );
}

bool testchk_box(zVec v, box_t box)
{
  return zVecElem(v,0) >= box[0] && zVecElem(v,0) <= box[2] &&
         zVecElem(v,1) >= box[1] && zVecElem(v,1) <= box[3] ?
    true : false;
}

bool testchk(zVec v, void *util)
{
  bool ret[4];

  ret[0] = testchk_box( v, ((box_t*)util)[0] );
  ret[1] = testchk_box( v, ((box_t*)util)[1] );
  ret[2] = testchk_box( v, ((box_t*)util)[2] );
  ret[3] = testchk_box( v, ((box_t*)util)[3] );
  return ret[0] || ret[1] || ret[2] || ret[3];
}

void output_box(FILE *fp, box_t box[])
{
  int i;

  for( i=0; i<4; i++ ){
    fprintf( fp, "%g %g\n", box[i][0], box[i][1] );
    fprintf( fp, "%g %g\n", box[i][2], box[i][1] );
    fprintf( fp, "%g %g\n", box[i][2], box[i][3] );
    fprintf( fp, "%g %g\n", box[i][0], box[i][3] );
    fprintf( fp, "%g %g\n\n", box[i][0], box[i][1] );
  }
}

#define W 10.0
#define HW (0.5 *W)
#define QW (0.25*W)
#define MW (0.1 *W)

int main(int argc, char *argv[])
{
  zRRT rrt;
  zRRTListCell *rc;
  box_t box[4];
  zVec min, max, start, goal;
  FILE *fp;

  zRandInit();
  box_create( box[0],    QW-zRandF(MW,HW),    QW-zRandF(MW,HW),    QW+zRandF(MW,HW),    QW+zRandF(MW,HW) );
  box_create( box[1], HW+QW-zRandF(MW,HW),    QW-zRandF(MW,HW), HW+QW+zRandF(MW,HW),    QW+zRandF(MW,HW) );
  box_create( box[2],    QW-zRandF(MW,HW), HW+QW-zRandF(MW,HW),    QW+zRandF(MW,HW), HW+QW+zRandF(MW,HW) );
  box_create( box[3], HW+QW-zRandF(MW,HW), HW+QW-zRandF(MW,HW), HW+QW+zRandF(MW,HW), HW+QW+zRandF(MW,HW) );

  fp = fopen( "a", "w" );
  output_box( fp, box );
  fclose( fp );

  start = zVecAlloc( 2 );
  goal = zVecAlloc( 2 );
  min = zVecCreateList( 2, 0.0, 0.0 );
  max = zVecCreateList( 2, W, W );
  zVecRand( start, min, max );

  zRRTInit( &rrt, min, max, 0.5, NULL, NULL, testchk, NULL );
  zRRTEscape( &rrt, start, 0, box, goal );

  fp = fopen( "s", "w" );
  fprintf( fp, "%g %g\n", zVecElem(start,0), zVecElem(start,1) );
  fclose( fp );
  fp = fopen( "p", "w" );
  zListForEach( &rrt.slist, rc ){
    if( !rc->data.parent ) continue;
    zVecValueFPrint( fp, rc->data.parent->v );
    zVecValueFPrint( fp, rc->data.v );
    fprintf( fp, "\n" );
  }
  fclose( fp );
  fp = fopen( "g", "w" );
  fprintf( fp, "%g %g\n", zVecElem(goal,0), zVecElem(goal,1) );
  fclose( fp );

  zRRTDestroy( &rrt );
  zVecFreeAtOnce( 4, min, max, start, goal );
  return 0;
}
