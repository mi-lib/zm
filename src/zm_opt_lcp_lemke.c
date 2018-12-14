/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_lcp_lemke - optimization tools:
 * Lemke's linear complementarity pivot method.
 */

#include <zm/zm_opt.h>

/* Lemke tableau */
typedef struct{
  zMat m;
  zVec q;
  zIndex ib, in;
  int act;
} _zLemke;

static _zLemke *_zLemkeCreate(_zLemke *tab, zMat m, zVec q);
static void _zLemkeDestroy(_zLemke *tab);
static bool _zLemkeInit(_zLemke *tab);
static void _zLemkeSweep(_zLemke *tab, int p);
static bool _zLemkeSwap(_zLemke *tab, int p);
static int _zLemkePivot(_zLemke *tab);
static bool _zLemkeIter(_zLemke *tab);
static void _zLemkeAnswer(_zLemke *tab, zVec w, zVec z);

#ifdef DEBUG
void _zLemkeWrite(_zLemke *tab);
#endif /* DEBUG */

/* (static)
 * _zLemkeCreate
 * - create Lemke tableau and lexicon.
 */
_zLemke *_zLemkeCreate(_zLemke *tab, zMat m, zVec q)
{
  register int i, j, n;

  n = zMatRowSizeNC(m);
  tab->m = zMatAlloc( n, 2*n+1 );
  tab->q = zVecClone( q );
  tab->ib = zIndexCreate( n );
  tab->in = zIndexCreate( n );
  if( !tab->m || !tab->q || !tab->ib || !tab->in ){
    _zLemkeDestroy( tab );
    return NULL;
  }
  tab->act = 2*n;
  zIndexOrder( tab->ib, n );
  zIndexOrder( tab->in, 0 );
  for( i=0; i<n; i++ ){
    zMatSetElem( tab->m, i, i+n, 1.0 );
    zMatSetElem( tab->m, i, tab->act, -1.0 );
  }
  for( i=0; i<n; i++ )
    for( j=0; j<n; j++ )
      zMatSetElem( tab->m, i, j, -zMatElem(m,i,j) );
  return tab;
}

/* (static)
 * _zLemkeDestroy
 * - destroy Lemke tableau and lexicon.
 */
void _zLemkeDestroy(_zLemke *tab)
{
  zMatFree( tab->m );
  zVecFree( tab->q );
  zIndexFree( tab->ib );
  zIndexFree( tab->in );
}

/* (static)
 * _zLemkeInit
 * - initialize Lemke tableau and lexicon.
 */
bool _zLemkeInit(_zLemke *tab)
{
  int p;
  double qmin;

  if( ( qmin = zVecMin( tab->q, &p ) ) >= 0 ) return true;
  _zLemkeSweep( tab, p );
  _zLemkeSwap( tab, p );
  return false;
}

/* (static)
 * _zLemkeSweep
 * - sweep-out the active column of Lemke tableau.
 */
void _zLemkeSweep(_zLemke *tab, int p)
{
  register int i, j;
  double ap;

  /* normalize pivot row */
  ap = zMatElem( tab->m, p, tab->act );
  for( j=0; j<zArrayNum(tab->in); j++ )
    zMatElem( tab->m, p, zIndexElem(tab->in,j) ) /= ap;
  zMatElem( tab->m, p, zIndexElem(tab->ib,p) ) /= ap;
  zMatSetElem( tab->m, p, tab->act, 1.0 );
  zVecElem(tab->q,p) /= ap;
  /* sweep-out rest row */
  for( i=0; i<zVecSizeNC(tab->q); i++ ){
    if( i == p ) continue;
    ap = zMatElem( tab->m, i, tab->act );
    for( j=0; j<zArrayNum(tab->in); j++ )
      zMatElem( tab->m, i, zIndexElem(tab->in,j) )
        -= zMatElem( tab->m, p, zIndexElem(tab->in,j) ) * ap;
    zMatElem( tab->m, i, zIndexElem(tab->ib,p) )
      =- zMatElem( tab->m, p, zIndexElem(tab->ib,p) ) * ap;
    zMatSetElem( tab->m, i, tab->act, 0.0 );
    zVecElem(tab->q,i) -= zVecElem(tab->q,p) * ap;
  }
}

/* (static)
 * _zLemkeSwap
 * - swap Lemke lexicon and active variable index.
 */
bool _zLemkeSwap(_zLemke *tab, int p)
{
  register int i, ib;

  ib = zIndexElem(tab->ib,p);
  zIndexSetElem( tab->ib, p, tab->act );
  /* choose complementary pivot */
  tab->act = ( ib + zArrayNum(tab->ib) ) % ( 2*zArrayNum(tab->ib) );
  for( i=0; i<zArrayNum(tab->in); i++ )
    if( zIndexElem(tab->in,i) == tab->act ){
      zIndexSetElem( tab->in, i, ib );
      break;
    }
  return ib == zMatColSizeNC(tab->m) - 1 ? true : false;
}

/* (static)
 * _zLemkePivot
 * - pivot Lemke tableau.
 */
int _zLemkePivot(_zLemke *tab)
{ /* next pivot */
  register int i, np;
  double a, p, p_min;

  for( p_min=HUGE_VAL, np=-1, i=0; i<zVecSizeNC(tab->q); i++ ){
    if( ( a = zMatElem(tab->m,i,tab->act) ) < zTOL )
      continue;
    p = zVecElem(tab->q,i) / a;
    if( p < p_min ){
      p_min = p;
      np = i;
    }
  }
  return np;
}

/* (static)
 * _zLemkeIter
 * - iterate Lemke's pivoting procedure.
 */
bool _zLemkeIter(_zLemke *tab)
{
  int p;

  do{
    if( ( p = _zLemkePivot( tab ) ) < 0 ){
      ZRUNWARN( ZM_ERR_OPT_UNSOLVE );
      return false;
    }
    _zLemkeSweep( tab, p );
  } while( !_zLemkeSwap( tab, p ) );
  return true;
}

/* (static)
 * _zLemkeAnswer
 * - complementary vectors.
 */
void _zLemkeAnswer(_zLemke *tab, zVec w, zVec z)
{
  register int i, idx;

  if( w ){
    zVecClear( w );
    for( i=0; i<zArrayNum(tab->ib); i++ ){
      idx = zIndexElem(tab->ib,i) - zVecSizeNC(w);
      if( idx >= 0 && idx < zVecSizeNC(w) )
        zVecSetElem( w, idx, zVecElem(tab->q,i) );
    }
  }
  zVecClear( z );
  for( i=0; i<zArrayNum(tab->ib); i++ )
    if( zIndexElem(tab->ib,i) < zVecSizeNC(z) )
      zVecSetElem( z, zIndexElem(tab->ib,i), zVecElem(tab->q,i) );
}

#ifdef DEBUG
/* (static)
 * _zLemkeWrite
 * - output tableau contents (for debug).
 */
void _zLemkeWrite(_zLemke *tab)
{ /* for debug. */
  printf( "tableau: " ); zMatWrite( tab->m );
  printf( "answer: " ); zVecWrite( tab->q );
  printf( "base lex.: " ); zIndexWrite( tab->ib );
  printf( "nonbase lex.: " ); zIndexWrite( tab->in );
  printf( "active var.: %d\n", tab->act );
}
#endif /* DEBUG */

/* zLCPSolveLemke
 * - solve linear complementarity problem by Lemke's method.
 */
bool zLCPSolveLemke(zMat m, zVec q, zVec w, zVec z)
{
  _zLemke tab;
  bool ret = true;

  if( !zMatIsSqr(m) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return false;
  }
  if( !zMatRowVecSizeIsEqual(m,q) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return false;
  }
  if( ( w && ( !zVecSizeIsEqual(w,q) || !zVecSizeIsEqual(w,z) ) ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return false;
  }
  if( !_zLemkeCreate( &tab, m, q ) ){
    ZALLOCERROR();
    return false;
  }
  if( ( ret = ( _zLemkeInit( &tab ) || _zLemkeIter( &tab ) ) ) )
    _zLemkeAnswer( &tab, w, z );
  _zLemkeDestroy( &tab );
  return ret;
}
