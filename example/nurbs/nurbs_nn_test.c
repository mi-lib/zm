#include <zm/zm_nurbs.h>

/* nearest neighbor on NURBS */

double _zNURBSVecNNDiv(zNURBS *nurbs, zVec v, zVec nn, int div)
{
  double s1, s2, s1old, s2old, si;
  double d, dmin1, dmin2;
  zVec vs;
  register int i;

  if( !( vs = zVecAlloc( zVecSizeNC(v) ) ) )
    return zNURBSKnot0(nurbs); /* dummy */
  s1 = zNURBSKnot0(nurbs);
  s2 = zNURBSKnotE(nurbs);
  dmin1 = dmin2 = HUGE_VAL;
  do{
    s1old = s1;
    s2old = s2;
    for( i=0; i<=div; i++ ){
      si = (s2old-s1old)*i/div + s1old;
      zNURBSVec( nurbs, si, vs );
      if( ( d = zVecDist( v, vs ) ) <= dmin1 ){
        dmin2 = dmin1; s2 = s1;
        dmin1 = d;     s1 = si;
      } else
      if( d <= dmin2 ){
        dmin2 = d;
        s2 = si;
      }
    }
  } while( !zIsTiny( s1 - s2 ) && !zIsTiny( dmin1 - dmin2 ) );
  zNURBSVec( nurbs, ( si = s1 ), nn );
  zVecFree( vs );
  return si;
}

double _zNURBSVecNNDivRef(zNURBS *nurbs, zVec v, double sr, zVec nn, int div)
{
  double s1, s2, s1old, s2old, si;
  double d, dmin1, dmin2;
  zVec vs;
  register int i;

  if( !( vs = zVecAlloc( zVecSizeNC(v) ) ) )
    return zNURBSKnot0(nurbs); /* dummy */
  s1 = zNURBSKnot0(nurbs);
  s2 = zNURBSKnotE(nurbs);
  dmin1 = dmin2 = HUGE_VAL;
  do{
    s1old = s1;
    s2old = s2;
    for( i=0; i<=div; i++ ){
      si = (s2old-s1old)*i/div + s1old;
      zNURBSVec( nurbs, si, vs );
      if( ( d = zVecDist( v, vs ) + zSqr(si-sr) ) <= dmin1 ){
        dmin2 = dmin1; s2 = s1;
        dmin1 = d;     s1 = si;
      } else
      if( d <= dmin2 ){
        dmin2 = d;
        s2 = si;
      }
    }
  } while( !zIsTiny( s1 - s2 ) && !zIsTiny( dmin1 - dmin2 ) );
  zNURBSVec( nurbs, ( si = s1 ), nn );
  zVecFree( vs );
  return si;
}



#define STEP 100
#define N     50

void test(zNURBS *nurbs)
{
  register int i;
  FILE *fp[5];
  zVec v, n, nn;
  double sr;

  v = zVecAlloc( zVecSizeNC( zNURBSCP(nurbs,0) ) );
  n = zVecAlloc( zVecSizeNC( zNURBSCP(nurbs,0) ) );
  nn = zVecAlloc( zVecSizeNC( zNURBSCP(nurbs,0) ) );

  fp[0] = fopen( "p", "w" );
  for( i=0; i<=STEP; i++ ){
    if( zNURBSVec( nurbs, (double)i/STEP, v ) )
      zVecDataFWrite( fp[0], v );
  }
  fclose( fp[0] );
  fp[0] = fopen( "nn", "w" );
  fp[1] = fopen( "nn1", "w" );
  fp[2] = fopen( "nn2", "w" );
  fp[3] = fopen( "nn3", "w" );
  fp[4] = fopen( "nnr", "w" );
  for( i=0; i<N; i++ ){
#if 0
    zVecSetElemList( v, zRandF(-1,7), zRandF(-3,5) );
    zNURBSVecNN( nurbs, v, nn );
    zVecDataFWrite( fp[0], v );
    zVecDataFWrite( fp[0], nn );
    fprintf( fp[0], "\n" );
#else
    sr = zNURBSKnotE(nurbs)*i/N + zNURBSKnot0(nurbs);
    zNURBSVecDiff( nurbs, sr, n, 2 );
    zNURBSVec( nurbs, sr, nn );
    zVecCat( nn, zRandF(-0.01,0.01), n, v );
    zVecDataFWrite( fp[0], v );
    zVecDataFWrite( fp[0], nn );
    fprintf( fp[0], "\n" );
#endif

    _zNURBSVecNNDiv( nurbs, v, nn, 10 );
    zVecDataFWrite( fp[1], v );
    zVecDataFWrite( fp[1], nn );
    fprintf( fp[1], "\n" );
    _zNURBSVecNNDiv( nurbs, v, nn, 20 );
    zVecDataFWrite( fp[2], v );
    zVecDataFWrite( fp[2], nn );
    fprintf( fp[2], "\n" );
    _zNURBSVecNNDiv( nurbs, v, nn, 30 );
    zVecDataFWrite( fp[3], v );
    zVecDataFWrite( fp[3], nn );
    fprintf( fp[3], "\n" );
    _zNURBSVecNNDivRef( nurbs, v, sr, nn, 30 );
    zVecDataFWrite( fp[4], v );
    zVecDataFWrite( fp[4], nn );
    fprintf( fp[4], "\n" );
  }
  fclose( fp[0] );
  fclose( fp[1] );
  fclose( fp[2] );
  fclose( fp[3] );
  fclose( fp[4] );

  zVecFree( v );
  zVecFree( n );
  zVecFree( nn );
}

void output_src(zSeq *seq)
{
  FILE *fp;
  zSeqListCell *cp;

  fp = fopen( "src", "w" );
  zListForEach( seq, cp )
    fprintf( fp, "%f %f\n", zVecElem(cp->data.v,0), zVecElem(cp->data.v,1) );
  fclose( fp );
}

#define DIM 3

int main(int argc, char *argv[])
{
  zNURBS nurbs;
  zSeq seq;
  zVec v;
  int num, i;
  /* example data array */
  double xp[] = { 2.0, 3.0, 5.0, 4.0, 5.0, 7.0 };
  double yp[] = { 3.0,-1.0,-2.0, 0.0, 4.0, 1.5 };

  zRandInit();
  /* creation of x-values and y-values vector */
  num = sizeof(xp) / sizeof(double);
  zListInit( &seq );
  for( i=0; i<num; i++ ){
    v = zVecCreateList( 2, xp[i], yp[i] );
    zSeqEnqueue( &seq, v, 1.0 );
  }
  output_src( &seq );

  /* creation of spline interpolator */
  if( zNURBSCreate( &nurbs, &seq, DIM ) ){
    zNURBSKnotNormalize( &nurbs );
    test( &nurbs );
    zNURBSDestroy( &nurbs );
  }
  zSeqFree( &seq );
  return 0;  
}
