#include <zm/zm_nle.h>

zVec zSSSolveSteffensen(zVec (* f)(zVec,zVec,void*), zVec x, void *util, int iter)
{
  register int i, j;
  zMat dx, ddx;
  zVec f0, f1;
  double l;

  dx = zMatAllocSqr( zVecSizeNC(x) );
  ddx = zMatAllocSqr( zVecSizeNC(x) );
  f0 = zVecAlloc( zVecSizeNC(x) );
  f1 = zVecAlloc( zVecSizeNC(x) );

  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    zVecCopyNC( x, f0 );
zVecPrint(x);
    for( j=0; ; j++ ){
      f( f0, f1, util );
printf("f1: "); zVecPrint(f1); getchar();
      if( j == zVecSizeNC(x) ) break;
      zRawVecSub( zVecBufNC(f1), zVecBufNC(f0), zMatRowBufNC(dx,j), zVecSizeNC(x) );
      if( zRawVecIsTiny( zMatRowBufNC(dx,j), zVecSizeNC(x) ) ){
        zVecCopyNC( f1, x );
        goto TERMINATE;
      }
      zSwap( zVec, f0, f1 );
    }
    if( zVecIsTiny( zVecSubNCDRC( f1, f0 ) ) ){
      zVecCopyNC( f1, x );
      goto TERMINATE;
    }
    l = zVecSqrNorm( f1 );
    for( j=0; j<zMatRowSizeNC(dx)-1; j++ ){
      zRawVecSub( zMatRowBufNC(dx,j+1), zMatRowBufNC(dx,j), zVecBufNC(f0), zVecSizeNC(f0) );
      zMatPutCol( ddx, j, f0 );
      zMatElemNC(ddx,j,j) += 1.0;
    }
    zRawVecSubDRC( zVecBufNC(f1), zMatRowBufNC(dx,j), zVecSizeNC(f1) );
    zMatPutCol( ddx, j, f1 );
    zMatElemNC(ddx,j,j) += 1.0;

    zMatGetCol( dx, 0, f1 );
zMatPrint(dx); zMatPrint(ddx); getchar();
    if( !zLESolveGauss( ddx, f1, f0 ) ){
      zVecCopyNC( f1, x );
      continue;
    }
    zMulMatTVec( dx, f0, f1 );
    zVecSubNCDRC( x, f1 );
  }
  ZITERWARN( iter );
 TERMINATE:
  zMatFree( dx );
  zMatFree( ddx );
  zVecFree( f0 );
  zVecFree( f1 );
  return x;
}


zVec f(zVec x, zVec y, void *util)
{
  zVecSetElem( y, 0, sqrt(zVecElem(x,0)) );
  zVecSetElem( y, 1, cos(zVecElem(x,0)) );
  zVecSetElem( y, 2, 2.0/(zSqr(zVecElem(x,0))+1.0)+1 );
  zVecSetElem( y, 3, 0.5*(zVecElem(x,0)+3) );
  zVecSetElem( y, 4, (1-zSqr(zVecElem(x,0)))/(1+sqrt(1+zSqr(zVecElem(x,0)))) );
  zVecSetElem( y, 5, 1.0/(1.0+zVecElem(x,0)) );
  return y;
}

int main(void)
{
  zVec x;

  x = zVecAlloc( 6 );

  zVecSetAll( x, 0.1 );
  zSSSolve( f, x, NULL, 0 );
  zVecPrint( x );
  printf( ">> assertion <<\n" );
  f( x, x, NULL );
  zVecPrint( x );

  zVecSetAll( x, 0.1 );
  zSSSolveSteffensen( f, x, NULL, 0 );
  zVecPrint( x );
  printf( ">> assertion <<\n" );
  f( x, x, NULL );
  zVecPrint( x );

  zVecFree( x );
  return 0;
}
