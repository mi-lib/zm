#include <zm/zm_opt.h>

static void _zLPTableauFPrint(FILE *fp, zLPTableau *tab)
{
  fprintf( fp, "A: " ); zMatFPrint( fp, tab->a );
  fprintf( fp, "b: " ); zVecFPrint( fp, tab->b );
  fprintf( fp, "c: " ); zVecFPrint( fp, tab->c );
  fprintf( fp, "d: = %f\n", tab->d );
  fprintf( fp, "(Ib): " ); zIndexFPrint( fp, tab->ib );
  fprintf( fp, "(In): " ); zIndexFPrint( fp, tab->in );
  fprintf( fp, "(Ir): " ); zIndexFPrint( fp, tab->ir );
}
#define _zLPTableauPrint( tab )  _zLPTableauFPrint( stdout, tab )

typedef struct{
  zMat qinv;
  zVec qinvc;
  zVec xtmp;
  zIndex ia;
  /* workspace */
  zMat _m;
  zVec _v;
  zVec _lambda;
} zQPASM;

bool _zQPASMInitBase(zQPASM *asm, zMat a, zVec b, zVec x)
{
  zLPTableau tab;
  register int i, k;
  int n;

  tab.a = zMatAlloc( zMatRowSizeNC(a), 2 * ( zMatColSizeNC(a) + zMatRowSizeNC(a) ) );
  tab.b = zVecAlloc( zMatRowSizeNC(a) );
  tab.c = zVecAlloc( zMatColSizeNC(tab.a) );
  tab.ib = zIndexCreate( zMatRowSizeNC(a) );
  tab.in = zIndexCreate( zMatColSizeNC(tab.a) - zMatRowSizeNC(a) );
  tab.ir = zIndexCreate( zMatRowSizeNC(a) );
  if( !tab.a || !tab.b ||!tab.c || !tab.ib || !tab.in || !tab.ir ){
    zLPTableauDestroy( &tab );
    return false;
  }
  n = 2 * zMatColSizeNC(a);
  for( i=0; i<zVecSizeNC(b); i++ ){
    if( zVecElemNC(b,i) >= 0 ){
      zRawVecCopy( zMatRowBufNC(a,i), zMatRowBufNC(tab.a,i), zMatColSizeNC(a) );
      zRawVecRev( zMatRowBufNC(a,i), zMatRowBufNC(tab.a,i)+zMatColSize(a), zMatColSizeNC(a) );
      zMatSetElemNC( tab.a, i, n+i,-1.0 );
      zVecSetElemNC( tab.b, i, zVecElemNC(b,i) );
    } else{
      zRawVecRev( zMatRowBufNC(a,i), zMatRowBufNC(tab.a,i), zMatColSizeNC(a) );
      zRawVecCopy( zMatRowBufNC(a,i), zMatRowBufNC(tab.a,i)+zMatColSize(a), zMatColSizeNC(a) );
      zMatSetElemNC( tab.a, i, n+i, 1.0 );
      zVecSetElemNC( tab.b, i,-zVecElemNC(b,i) );
    }
    zMatSetElemNC( tab.a, i, zArraySize(tab.in)+i, 1.0 );
  }
  for( i=zArraySize(tab.in); i<zVecSizeNC(tab.c); i++ )
    zVecSetElemNC( tab.c, i, 1.0 );
  tab.d = 0;
  zIndexOrder( tab.ib, zArraySize(tab.in) );
  zIndexOrder( tab.in, 0 );
  zIndexOrder( tab.ir, 0 );
  if( !zLPTableauSimplex( &tab ) || !zIsTiny(tab.d) ){
    ZRUNWARN( ZM_ERR_OPT_UNSOLVE );
    goto TERMINATE;
  }
  zLPTableauFindBase( &tab );

  zVecZero( x );
  for( i=0; i<zArraySize(tab.ib); i++ ){
    if( zIndexElemNC(tab.ib,i) < zMatColSizeNC(a) )
      zVecSetElemNC( x, zIndexElemNC(tab.ib,i), zVecElemNC(tab.b,i) );
    else
    if( zIndexElemNC(tab.ib,i) < n )
      zVecSetElemNC( x, zIndexElemNC(tab.ib,i)-zMatColSizeNC(a),-zVecElemNC(tab.b,i) );
    else
    if( ( k = zIndexElemNC(tab.ib,i) - n ) < zMatRowSizeNC(a) &&
        zVecElemNC(tab.b,i) > 0 ){
      zIndexRemove( asm->ia, k );
    }
  }
_zLPTableauPrint( &tab );
zVecPrint( x );
zIndexPrint( asm->ia );

 TERMINATE:
  zLPTableauDestroy( &tab );
  return true;
}

static zQPASM *_zQPASMInit(zQPASM *asm, zMat q, zVec c, zMat a, zVec b, zVec x)
{
  asm->qinv = zMatAllocSqr( zVecSizeNC(c) );
  asm->qinvc = zVecAlloc( zVecSizeNC(c) );
  asm->xtmp = zVecAlloc( zVecSizeNC(c) );
  asm->ia = zIndexAlloc( zMatRowSizeNC(a) );
  asm->_m = zMatAllocSqr( zMatRowSizeNC(a) );
  asm->_v = zVecAlloc( zMatRowSizeNC(a) );
  asm->_lambda = zVecAlloc( zMatRowSizeNC(a) );

  if( !asm->qinv || !asm->qinvc || !asm->xtmp || !asm->ia ||
      !asm->_m || !asm->_v || !asm->_lambda ) return NULL;
  zMatInv( q, asm->qinv );
  zMulMatVec( asm->qinv, c, asm->qinvc );

  _zQPASMInitBase( asm, a, b, x );
  return asm;
}

static bool _zQPASMSolveEq(zQPASM *asm, zMat a, zVec b)
{
  register int i;

  zMatSetSize( asm->_m, zArraySize(asm->ia), zArraySize(asm->ia) );
  zVecSetSize( asm->_v, zArraySize(asm->ia) );
  zVecSetSize( asm->_lambda, zArraySize(asm->ia) );
  for( i=0; i<zArraySize(asm->ia); i++ ){

  }
  return true;
}

static void _zQPASMDestroy(zQPASM *asm)
{
  zMatFree( asm->qinv );
  zVecFree( asm->qinvc );
  zVecFree( asm->xtmp );
  zIndexFree( asm->ia );
  zMatFree( asm->_m );
  zVecFree( asm->_v );
  zVecFree( asm->_lambda );
}


bool _zQPSolveASM(zMat q, zVec c, zMat a, zVec b, zVec ans, double *cost)
{
  zQPASM asm;

  _zQPASMInit( &asm, q, c, a, b, ans );
  _zQPASMSolveEq( &asm, a, b );
  _zQPASMDestroy( &asm );
  return true;
}

int main(void)
{
  double qarray[] = {
    1, 0, 0, 1,
  };
  double carray[] = {
    1, -6,
  };
#if 0
  double aarray[] = {
    1,-1,
   -2,-1,
   -1,-3,
   -1, 1,
    3, 4,
  };
  double barray[] = {
   -4, -9, -12, -3, -5,
  };
  int n = 2, m = 5;
#else
  double aarray[] = {
    1.0, 1.0,
   -1.0,-4.0,
   -5.0, 6.0,
    5.0,-2.0,
  };
  double barray[] = {
    10.0, -60.0, -25.0, -6.0,
  };
  int n = 2, m = 4;
#endif
  zMat q, a;
  zVec c, b, x;
  double cost;

  q = zMatCloneArray( qarray, n, n );
  c = zVecCloneArray( carray, n );
  a = zMatCloneArray( aarray, m, n );
  b = zVecCloneArray( barray, m );
  x = zVecAlloc( n );

  _zQPSolveASM( q, c, a, b, x, &cost );

  zMatFree( q );
  zMatFree( a );
  zVecFree( c );
  zVecFree( b );
  zVecFree( x );
  return 0;
}
