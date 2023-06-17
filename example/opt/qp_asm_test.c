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
  bool *is_active;
  int n_active;
  /* workspace */
  zMat _m;
  zVec _v;
  zVec _lambda;
} zQPASM;

bool _zQPASMInitBase(zQPASM *qpasm, zMat a, zVec b, zVec x)
{
  zLPTableau tab;
  int i, k, n;

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
    qpasm->is_active[i] = false;
  }
  for( i=zArraySize(tab.in); i<zVecSizeNC(tab.c); i++ )
    zVecSetElemNC( tab.c, i, 1.0 );
  tab.d = 0;
  zIndexOrder( tab.ib, zArraySize(tab.in) );
  zIndexOrder( tab.in, 0 );
  zIndexOrder( tab.ir, 0 );
  _zLPTableauPrint( &tab );
  if( !zLPTableauSimplex( &tab ) || !zIsTiny(tab.d) ){
    ZRUNWARN( ZM_ERR_OPT_UNSOLVE );
    goto TERMINATE;
  }
  zLPTableauFindBase( &tab );

  zVecZero( x );
  for( i=0; i<zArraySize(tab.ib); i++ ){
    if( zIndexElemNC(tab.ib,i) < zMatColSizeNC(a) ){
      zVecSetElemNC( x, zIndexElemNC(tab.ib,i), zVecElemNC(tab.b,i) );
      qpasm->is_active[zIndexElemNC(tab.ir,i)] = true;
    }
    else
    if( zIndexElemNC(tab.ib,i) < n ){
      zVecSetElemNC( x, zIndexElemNC(tab.ib,i)-zMatColSizeNC(a),-zVecElemNC(tab.b,i) );
      qpasm->is_active[zIndexElemNC(tab.ir,i)] = true;
    }
    else
    if( ( k = zIndexElemNC(tab.ib,i) - n ) < zMatRowSizeNC(a) &&
        zVecElemNC(tab.b,i) > 0 ){
    }
  }
  _zLPTableauPrint( &tab );
  zVecPrint( x );
  for( i=0; i<zMatRowSizeNC(a); i++ ){
    printf( "constraint #%d: %s\n", i, qpasm->is_active[i] ? "active" : "inactive" );
  }

 TERMINATE:
  zLPTableauDestroy( &tab );
  return true;
}

static zQPASM *_zQPASMInit(zQPASM *qpasm, zMat q, zVec c, zMat a, zVec b, zVec x)
{
  qpasm->qinv = zMatAllocSqr( zVecSizeNC(c) );
  qpasm->qinvc = zVecAlloc( zVecSizeNC(c) );
  qpasm->xtmp = zVecAlloc( zVecSizeNC(c) );
  qpasm->is_active = zAlloc( bool, zMatRowSizeNC(a) );
  qpasm->_m = zMatAllocSqr( zMatRowSizeNC(a) );
  qpasm->_v = zVecAlloc( zMatRowSizeNC(a) );
  qpasm->_lambda = zVecAlloc( zMatRowSizeNC(a) );

  if( !qpasm->qinv || !qpasm->qinvc || !qpasm->xtmp || !qpasm->is_active ||
      !qpasm->_m || !qpasm->_v || !qpasm->_lambda ) return NULL;
  qpasm->n_active = 0;
  zMatInv( q, qpasm->qinv );
  zMulMatVec( qpasm->qinv, c, qpasm->qinvc );

  _zQPASMInitBase( qpasm, a, b, x );
  return qpasm;
}

static bool _zQPASMSolveEq(zQPASM *qpasm, zMat a, zVec b)
{
  int i;

  zMatSetSize( qpasm->_m, qpasm->n_active, qpasm->n_active );
  zVecSetSize( qpasm->_v, qpasm->n_active );
  zVecSetSize( qpasm->_lambda, qpasm->n_active );
  for( i=0; i<qpasm->n_active; i++ ){

  }
  return true;
}

static void _zQPASMDestroy(zQPASM *qpasm)
{
  zMatFree( qpasm->qinv );
  zVecFree( qpasm->qinvc );
  zVecFree( qpasm->xtmp );
  zFree( qpasm->is_active );
  zMatFree( qpasm->_m );
  zVecFree( qpasm->_v );
  zVecFree( qpasm->_lambda );
}


bool _zQPSolveASM(zMat q, zVec c, zMat a, zVec b, zVec ans, double *cost)
{
  zQPASM qpasm;

  _zQPASMInit( &qpasm, q, c, a, b, ans );
  _zQPASMSolveEq( &qpasm, a, b );
  _zQPASMDestroy( &qpasm );
  return true;
}

#define TEST 0

int main(void)
{
  double qarray[] = {
    1, 0, 0, 1,
  };
  double carray[] = {
    1, -6,
  };
#if TEST == 0
  double aarray[] = {
   -1,-3,
   -1, 1,
    3, 4,
    1,-1,
   -2,-1,
    3, 4,
  };
  double barray[] = {
   -12, -3, -5, -4, -9, -5,
  };
  int n = 2, m = 6;
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
