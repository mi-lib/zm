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

#ifdef DEBUG
/* print out tableau contents (for debug). */
static void _zLemkePrint(_zLemke *tab)
{ /* for debug. */
  printf( "tableau: " ); zMatPrint( tab->m );
  printf( "answer: " ); zVecPrint( tab->q );
  printf( "base lex.: " ); zIndexPrint( tab->ib );
  printf( "nonbase lex.: " ); zIndexPrint( tab->in );
  printf( "active var.: %d\n", tab->act );
}
#endif /* DEBUG */

/* destroy Lemke tableau and lexicon. */
static void _zLemkeDestroy(_zLemke *tab)
{
  zMatFree( tab->m );
  zVecFree( tab->q );
  zIndexFree( tab->ib );
  zIndexFree( tab->in );
}

/* create Lemke tableau and lexicon. */
static _zLemke *_zLemkeCreate(_zLemke *tab, zMat m, zVec q)
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
    zMatSetElemNC( tab->m, i, i+n, 1.0 );
    zMatSetElemNC( tab->m, i, tab->act, -1.0 );
  }
  for( i=0; i<n; i++ )
    for( j=0; j<n; j++ )
      zMatSetElemNC( tab->m, i, j, -zMatElemNC(m,i,j) );
  return tab;
}

/* sweep-out the active column of Lemke tableau. */
static void _zLemkeSweep(_zLemke *tab, int p)
{
  register int i, j;
  double ap;

  /* normalize pivot row */
  ap = zMatElemNC( tab->m, p, tab->act );
  for( j=0; j<zArraySize(tab->in); j++ )
    zMatElemNC( tab->m, p, zIndexElemNC(tab->in,j) ) /= ap;
  zMatElemNC( tab->m, p, zIndexElemNC(tab->ib,p) ) /= ap;
  zMatSetElemNC( tab->m, p, tab->act, 1.0 );
  zVecElemNC(tab->q,p) /= ap;
  /* sweep-out rest row */
  for( i=0; i<zVecSizeNC(tab->q); i++ ){
    if( i == p ) continue;
    ap = zMatElemNC( tab->m, i, tab->act );
    for( j=0; j<zArraySize(tab->in); j++ )
      zMatElemNC( tab->m, i, zIndexElemNC(tab->in,j) )
        -= zMatElemNC( tab->m, p, zIndexElemNC(tab->in,j) ) * ap;
    zMatElemNC( tab->m, i, zIndexElemNC(tab->ib,p) )
      =- zMatElemNC( tab->m, p, zIndexElemNC(tab->ib,p) ) * ap;
    zMatSetElemNC( tab->m, i, tab->act, 0.0 );
    zVecElemNC(tab->q,i) -= zVecElemNC(tab->q,p) * ap;
  }
}

/* swap Lemke lexicon and active variable index. */
static bool _zLemkeSwap(_zLemke *tab, int p)
{
  register int i, ib;

  ib = zIndexElemNC(tab->ib,p);
  zIndexSetElemNC( tab->ib, p, tab->act );
  /* choose complementary pivot */
  tab->act = ( ib + zArraySize(tab->ib) ) % ( 2*zArraySize(tab->ib) );
  for( i=0; i<zArraySize(tab->in); i++ )
    if( zIndexElemNC(tab->in,i) == tab->act ){
      zIndexSetElemNC( tab->in, i, ib );
      break;
    }
  return ib == zMatColSizeNC(tab->m) - 1 ? true : false;
}

/* find next pivot in Lemke tableau. */
static int _zLemkePivot(_zLemke *tab)
{
  register int i, np;
  double a, p, p_min;

  for( p_min=HUGE_VAL, np=-1, i=0; i<zVecSizeNC(tab->q); i++ ){
    if( ( a = zMatElemNC(tab->m,i,tab->act) ) < zTOL )
      continue;
    p = zVecElemNC(tab->q,i) / a;
    if( p < p_min ){
      p_min = p;
      np = i;
    }
  }
  return np;
}

/* iterate Lemke's pivoting procedure. */
static bool _zLemkeIter(_zLemke *tab)
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

/* complementary vectors. */
static void _zLemkeAnswer(_zLemke *tab, zVec w, zVec z)
{
  register int i, idx;

  if( w ){
    zVecZero( w );
    for( i=0; i<zArraySize(tab->ib); i++ ){
      idx = zIndexElemNC(tab->ib,i) - zVecSizeNC(w);
      if( idx >= 0 && idx < zVecSizeNC(w) )
        zVecSetElemNC( w, idx, zVecElemNC(tab->q,i) );
    }
  }
  zVecZero( z );
  for( i=0; i<zArraySize(tab->ib); i++ )
    if( zIndexElemNC(tab->ib,i) < zVecSizeNC(z) )
      zVecSetElemNC( z, zIndexElemNC(tab->ib,i), zVecElemNC(tab->q,i) );
}

/* initialize Lemke tableau and lexicon. */
static bool _zLemkeInit(_zLemke *tab)
{
  int p;
  double qmin;

  if( ( qmin = zVecMin( tab->q, &p ) ) >= 0 ) return true;
  _zLemkeSweep( tab, p );
  _zLemkeSwap( tab, p );
  return false;
}

/* solve linear complementarity problem by Lemke's method. */
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
