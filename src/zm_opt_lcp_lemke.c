/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_lcp_lemke - optimization tools: Lemke's linear complementarity pivot method.
 */

#include <zm/zm_opt.h>

/* Lemke tableau */
typedef struct{
  zMat tab;
  zVec p;
  zIndex ib, in;
  int act;
} _zLemke;

#ifdef DEBUG
/* print out tableau contents (for debug). */
static void _zLemkePrint(const _zLemke *lemke)
{ /* for debug. */
  printf( "tableau: " );      zMatPrint( lemke->tab );
  printf( "answer: " );       zVecPrint( lemke->p );
  printf( "base lex.: " );    zIndexPrint( lemke->ib );
  printf( "nonbase lex.: " ); zIndexPrint( lemke->in );
  printf( "active var.: %d\n", lemke->act );
}
#endif /* DEBUG */

/* destroy Lemke tableau and lexicon. */
static void _zLemkeDestroy(_zLemke *lemke)
{
  zMatFree( lemke->tab );
  zVecFree( lemke->p );
  zIndexFree( lemke->ib );
  zIndexFree( lemke->in );
}

/* create Lemke tableau and lexicon. */
static _zLemke *_zLemkeCreate(_zLemke *lemke, const zMat m, const zVec q)
{
  int i, j, n;

  n = zMatRowSizeNC(m);
  lemke->tab = zMatAlloc( n, 2*n+1 );
  lemke->p = zVecClone( q );
  lemke->ib = zIndexCreate( n );
  lemke->in = zIndexCreate( n );
  if( !lemke->tab || !lemke->p || !lemke->ib || !lemke->in ){
    _zLemkeDestroy( lemke );
    return NULL;
  }
  lemke->act = 2*n;
  zIndexOrder( lemke->ib, n );
  zIndexOrder( lemke->in, 0 );
  for( i=0; i<n; i++ ){
    for( j=0; j<n; j++ )
      zMatSetElemNC( lemke->tab, i, j, -zMatElemNC(m,i,j) );
    zMatSetElemNC( lemke->tab, i, i+n, 1.0 );
    zMatSetElemNC( lemke->tab, i, lemke->act, -1.0 );
  }
  return lemke;
}

/* sweep-out the active column of Lemke tableau. */
static void _zLemkeSweep(_zLemke *lemke, int p)
{
  int i, j;
  double ap;

  /* normalize pivot row */
  ap = zMatElemNC( lemke->tab, p, lemke->act );
  for( j=0; j<zIndexSizeNC(lemke->in); j++ )
    zMatElemNC( lemke->tab, p, zIndexElemNC(lemke->in,j) ) /= ap;
  zMatElemNC( lemke->tab, p, zIndexElemNC(lemke->ib,p) ) /= ap;
  zMatSetElemNC( lemke->tab, p, lemke->act, 1.0 );
  zVecElemNC(lemke->p,p) /= ap;
  /* sweep-out rest row */
  for( i=0; i<zVecSizeNC(lemke->p); i++ ){
    if( i == p ) continue;
    ap = zMatElemNC( lemke->tab, i, lemke->act );
    for( j=0; j<zIndexSizeNC(lemke->in); j++ )
      zMatElemNC( lemke->tab, i, zIndexElemNC(lemke->in,j) )
        -= zMatElemNC( lemke->tab, p, zIndexElemNC(lemke->in,j) ) * ap;
    zMatElemNC( lemke->tab, i, zIndexElemNC(lemke->ib,p) )
      =- zMatElemNC( lemke->tab, p, zIndexElemNC(lemke->ib,p) ) * ap;
    zMatSetElemNC( lemke->tab, i, lemke->act, 0.0 );
    zVecElemNC(lemke->p,i) -= zVecElemNC(lemke->p,p) * ap;
  }
}

/* swap Lemke lexicon and active variable index. */
static bool _zLemkeSwap(_zLemke *lemke, int p)
{
  int i, ib;

  ib = zIndexElemNC(lemke->ib,p);
  zIndexSetElemNC( lemke->ib, p, lemke->act );
  /* choose complementary pivot */
  lemke->act = ( ib + zIndexSizeNC(lemke->ib) ) % ( 2*zIndexSizeNC(lemke->ib) );
  for( i=0; i<zIndexSizeNC(lemke->in); i++ )
    if( zIndexElemNC(lemke->in,i) == lemke->act ){
      zIndexSetElemNC( lemke->in, i, ib );
      break;
    }
  return ib == zMatColSizeNC(lemke->tab) - 1 ? true : false;
}

/* find next pivot in Lemke tableau. */
static int _zLemkePivot(_zLemke *lemke)
{
  int i, np;
  double a, p, p_min;

  for( p_min=HUGE_VAL, np=-1, i=0; i<zVecSizeNC(lemke->p); i++ ){
    if( ( a = zMatElemNC(lemke->tab,i,lemke->act) ) < zTOL )
      continue;
    p = zVecElemNC(lemke->p,i) / a;
    if( p < p_min ){
      p_min = p;
      np = i;
    }
  }
  return np;
}

/* iterate Lemke's pivoting procedure. */
static bool _zLemkeIter(_zLemke *lemke)
{
  int p;

  do{
    if( ( p = _zLemkePivot( lemke ) ) < 0 ){
      ZRUNWARN( ZM_ERR_OPT_UNSOLVABLE );
      return false;
    }
    _zLemkeSweep( lemke, p );
  } while( !_zLemkeSwap( lemke, p ) );
  return true;
}

/* complementary vectors. */
static void _zLemkeAnswer(_zLemke *lemke, zVec w, zVec z)
{
  int i, idx;

  if( w ){
    zVecZero( w );
    for( i=0; i<zIndexSizeNC(lemke->ib); i++ ){
      idx = zIndexElemNC(lemke->ib,i) - zVecSizeNC(w);
      if( idx >= 0 && idx < zVecSizeNC(w) )
        zVecSetElemNC( w, idx, zVecElemNC(lemke->p,i) );
    }
  }
  zVecZero( z );
  for( i=0; i<zIndexSizeNC(lemke->ib); i++ )
    if( zIndexElemNC(lemke->ib,i) < zVecSizeNC(z) )
      zVecSetElemNC( z, zIndexElemNC(lemke->ib,i), zVecElemNC(lemke->p,i) );
}

/* initialize Lemke tableau and lexicon. */
static bool _zLemkeInit(_zLemke *lemke)
{
  int p;
  double qmin;

  if( ( qmin = zVecElemMin( lemke->p, &p ) ) >= 0 ) return true;
  _zLemkeSweep( lemke, p );
  _zLemkeSwap( lemke, p );
  return false;
}

/* solve linear complementarity problem by Lemke's method. */
bool zLCPSolveLemke(const zMat m, const zVec p, zVec w, zVec z)
{
  _zLemke lemke;
  bool ret = true;

  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return false;
  }
  if( !zMatRowVecSizeEqual( m, p ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return false;
  }
  if( ( w && ( !zVecSizeEqual( w, p ) || !zVecSizeEqual( w, z ) ) ) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return false;
  }
  if( !_zLemkeCreate( &lemke, m, p ) ){
    ZALLOCERROR();
    return false;
  }
  if( ( ret = ( _zLemkeInit( &lemke ) || _zLemkeIter( &lemke ) ) ) )
    _zLemkeAnswer( &lemke, w, z );
  _zLemkeDestroy( &lemke );
  return ret;
}
