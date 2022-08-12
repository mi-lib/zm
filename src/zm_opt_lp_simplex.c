/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_lp_simplex - optimization tools:
 * linear programming by simplex method.
 */

#include <zm/zm_opt.h>

/* simplex method tableau */
typedef struct{
  zMat a;
  zVec b, c;
  double d;
  zIndex ib, in;
  zIndex ir;
} _zLPTableau;

#ifdef DEBUG
/* print out tableau contents to a file (for debug). */
static void _zLPTableauFPrint(FILE *fp, _zLPTableau *tab)
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
#endif /* DEBUG */

/* create initial simplex tableau with slack variables. */
static bool _zLPTableauCreate(_zLPTableau *tab, zMat a, zVec b)
{
  register int i;

  tab->a = zMatAlloc( zMatRowSizeNC(a), zMatColSizeNC(a)+zMatRowSizeNC(a) );
  tab->b = zVecAlloc( zMatRowSizeNC(a) );
  tab->c = zVecAlloc( zMatColSizeNC(tab->a) );
  tab->ib = zIndexCreate( zMatRowSizeNC(a) );
  tab->in = zIndexCreate( zMatColSizeNC(a) );
  tab->ir = zIndexCreate( zMatRowSizeNC(a) );
  if( !tab->a || !tab->b ||!tab->c || !tab->ib || !tab->in || !tab->ir )
    return false;
  for( i=0; i<zVecSizeNC(b); i++ ){
    if( zVecElemNC(b,i) >= 0 ){
      zRawVecCopy( zMatRowBufNC(a,i), zMatRowBufNC(tab->a,i), zMatColSizeNC(a) );
      zVecSetElemNC( tab->b, i, zVecElemNC(b,i) );
    } else{
      zRawVecRev( zMatRowBufNC(a,i), zMatRowBufNC(tab->a,i), zMatColSizeNC(a) );
      zVecSetElemNC( tab->b, i, -zVecElemNC(b,i) );
    }
    zMatSetElemNC( tab->a, i, zMatColSizeNC(a)+i, 1.0 );
  }
  for( i=zMatColSizeNC(a); i<zVecSizeNC(tab->c); i++ )
    zVecSetElemNC( tab->c, i, 1.0 );
  tab->d = 0;
  zIndexOrder( tab->ib, zMatColSizeNC(a) );
  zIndexOrder( tab->in, 0 );
  zIndexOrder( tab->ir, 0 );
  return true;
}

/* destroy simplex tableau. */
static void _zLPTableauDestroy(_zLPTableau *tab)
{
  zMatFree( tab->a );
  zVecFree( tab->b );
  zVecFree( tab->c );
  zIndexFree( tab->ib );
  zIndexFree( tab->in );
  zIndexFree( tab->ir );
}

/* sweep-out one coefficient of cost function corresponding to old bases. */
static void _zLPTableauSweepC1(_zLPTableau *tab, int p)
{
  register int i;
  double cp;

  cp = zVecElemNC(tab->c,zIndexElemNC(tab->ib,p));
  for( i=0; i<zArraySize(tab->in); i++ )
    zVecElemNC(tab->c,zIndexElemNC(tab->in,i))
      -= zMatElemNC(tab->a,zIndexElemNC(tab->ir,p),zIndexElemNC(tab->in,i)) * cp;
  tab->d += zVecElemNC(tab->b,zIndexElemNC(tab->ir,p)) * cp;
  zVecSetElemNC( tab->c, zIndexElemNC(tab->ib,p), 0 );
}

/* sweep-out coefficients of cost function corresponding to bases. */
static void _zLPTableauSweepC(_zLPTableau *tab)
{
  register int i;

  for( i=0; i<zArraySize(tab->ib); i++ )
    _zLPTableauSweepC1( tab, i );
}

/* find next axis column of tableau to be sweeped-out. */
static int _zLPTableauFindNA(_zLPTableau *tab)
{ /* next axis for normal case */
  register int i, na;
  double c_min;

  c_min = zVecElemNC( tab->c, zIndexElemNC(tab->in,(na=0)) );
  for( i=1; i<zArraySize(tab->in); i++ )
    if( zVecElemNC(tab->c,zIndexElemNC(tab->in,i)) < c_min )
      c_min = zVecElemNC( tab->c, zIndexElemNC(tab->in,(na=i)) );
  return c_min < -zTOL ? na : -1;
}

/* find next axis column of tableau in degenerate case. */
static int _zLPTableauFindNA_deg(_zLPTableau *tab)
{ /* next axis for degenerated case */
  register int i, na = -1, idx_min = INT_MAX;

  for( i=0; i<zArraySize(tab->in); i++ )
    if( zVecElemNC(tab->c,zIndexElemNC(tab->in,i)) < -zTOL ){
      if( zIndexElemNC(tab->in,i) < idx_min ){
        idx_min = zIndexElemNC(tab->in,i);
        na = i;
      }
    }
  return na;
}

/* find next pivot in axis column to be sweeped-out. */
static int _zLPTableauFindNP(_zLPTableau *tab, int *na)
{ /* next pivot */
  register int i, np;
  double a, p, p_min;
  bool f_try = true;

 RETRY:
  for( p_min=HUGE_VAL, np=-1, i=0; i<zArraySize(tab->ir); i++ ){
    /* a should be comared with 0 in theory. But zTOL instead of 0 does work. */
    if( ( a = zMatElemNC(tab->a,zIndexElemNC(tab->ir,i),zIndexElemNC(tab->in,*na)) ) < zTOL )
      continue;
    p = zVecElemNC(tab->b,zIndexElemNC(tab->ir,i)) / a;
    if( zIsTiny( p ) && f_try ){ /* degenerated case */
      *na = _zLPTableauFindNA_deg( tab );
      f_try = false;
      goto RETRY;
    }
    if( p < p_min ){
      p_min = p;
      np = i;
    }
  }
  return np;
}

/* sweep-out tabeau matrix of constraint equation. */
static void _zLPTableauSweepA(_zLPTableau *tab, int np, int na)
{
  register int i, j, r;
  double ap;

  /* normalize pivot row */
  ap = zMatElemNC(tab->a,zIndexElemNC(tab->ir,np),zIndexElemNC(tab->in,na));
  for( j=0; j<zArraySize(tab->in); j++ )
    zMatElemNC(tab->a,zIndexElemNC(tab->ir,np),zIndexElemNC(tab->in,j)) /= ap;
  zMatElemNC(tab->a,zIndexElemNC(tab->ir,np),zIndexElemNC(tab->ib,np)) /= ap;
  zVecElemNC(tab->b,zIndexElemNC(tab->ir,np)) /= ap;
  /* sweep-out rest row */
  for( i=0; i<zArraySize(tab->ir); i++ ){
    if( i == np ) continue;
    r = zIndexElemNC(tab->ir,i);
    ap = zMatElemNC(tab->a,r,zIndexElemNC(tab->in,na));
    for( j=0; j<zArraySize(tab->in); j++ )
      zMatElemNC(tab->a,r,zIndexElemNC(tab->in,j))
        -= zMatElemNC(tab->a,zIndexElemNC(tab->ir,np),zIndexElemNC(tab->in,j)) * ap;
    zMatElemNC(tab->a,r,zIndexElemNC(tab->ib,np))
      -= zMatElemNC(tab->a,zIndexElemNC(tab->ir,np),zIndexElemNC(tab->ib,np)) * ap;
    zVecElemNC(tab->b,r) -= zVecElemNC(tab->b,zIndexElemNC(tab->ir,np)) * ap;
  }
}

/* swap base/non-base pivot. */
static void _zLPTableauSwapPivot(_zLPTableau *tab, int np, int na)
{
  int tmp;

  tmp = zIndexElemNC( tab->in, na );
  zIndexSetElemNC( tab->in, na, zIndexElemNC(tab->ib,np) );
  zIndexSetElemNC( tab->ib, np, tmp );
}

/* simplex method for initialized tableau.
 *      1  0 a_1(m+1) ... a_1n     b_0
 *   A=  .    .             .   b=  .
 *        .   .             .       .
 *      0  1 a_m(m+1) ... a_mn     b_m
 *
 *   C= [c_1 ...          c_n]  d=  0
 *
 *   Ib= [1 ... m], In= [m+1 ... n+m]
 */
static bool _zLPTableauSimplex(_zLPTableau *tab)
{
  int na, np;
  register int i = 0;

  _zLPTableauSweepC( tab );
  while( ( na = _zLPTableauFindNA( tab ) ) >= 0 ){
    if( i++ > zMatColSizeNC(tab->a) ) /* probably degenerated case */
      return false;
    if( ( np = _zLPTableauFindNP( tab, &na ) ) < 0 ) /* non-convex case */
      return false;
    _zLPTableauSweepA( tab, np, na );
    _zLPTableauSwapPivot( tab, np, na );
    _zLPTableauSweepC1( tab, np );
  }
  return true;
}

/* reset tableau to the second stage. */
static bool _zLPTableauReset(_zLPTableau *tab, zVec c)
{
  register int i, j, n;

  tab->d = 0; /* precautionary touch-up */
  n = zVecSizeNC(tab->c) - zVecSizeNC(tab->b);
  for( i=0; i<zArraySize(tab->ib); i++ ){
    if( zIndexElemNC(tab->ib,i) < n ) continue;
    for( j=0; j<zArraySize(tab->in); j++ ){
      if( zIndexElemNC(tab->in,j) < n && !zIsTiny( zMatElemNC(tab->a,i,zIndexElemNC(tab->in,j)) ) ){
        _zLPTableauSweepA( tab, i, j );
        _zLPTableauSwapPivot( tab, i, j );
        goto NEXT;
      }
    }
    zIndexRemove( tab->ir, i );
    zIndexRemove( tab->ib, i );
    i--;
   NEXT: ;
  }
  /* rearrange nonfeasible bases */
  for( i=0; i<zArraySize(tab->in); i++ )
    if( zIndexElemNC(tab->in,i) >= n ){
      zIndexRemove( tab->in, i );
      i--;
    }
  zVecSetSize( tab->c, zVecSizeNC(c) );
  zVecCopyNC( c, tab->c );
  return true;
}

/* arrange answer vector. */
static void _zLPTableauAns(_zLPTableau *tab, zVec ans)
{
  register int i;

  zVecZero( ans );
  for( i=0; i<zArraySize(tab->ib); i++ )
    zVecSetElemNC( ans, zIndexElemNC(tab->ib,i), zVecElemNC(tab->b,i) );
}

/* dual-phase simplex method for linear programming. */
bool zLPSolveSimplex(zMat a, zVec b, zVec c, zVec ans, double *cost)
{
  _zLPTableau tab;
  bool ret = false;

  if( !zMatColVecSizeIsEqual(a,ans) || !zMatRowVecSizeIsEqual(a,b) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return false;
  }
  if( !zVecSizeIsEqual(c,ans) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return false;
  }
  if( !_zLPTableauCreate( &tab, a, b ) ){
    ret = false;
    goto TERMINATE;
  }
  /* first phase: feasible base */
  if( !_zLPTableauSimplex( &tab ) || !zIsTiny(tab.d) ){
    ZRUNWARN( ZM_ERR_OPT_UNSOLVE );
    goto TERMINATE;
  }
  _zLPTableauReset( &tab, c );
  /* second phase: body problem */
  if( !_zLPTableauSimplex( &tab ) ){
    ZRUNERROR( ZM_ERR_OPT_INF );
    goto TERMINATE;
  }
  if( cost ) *cost = tab.d;
  _zLPTableauAns( &tab, ans );
  ret = true;
 TERMINATE:
  _zLPTableauDestroy( &tab );
  return ret;
}

/* find a feasible base under Ax=b and x>=0 based on simplex method. */
bool zLPFeasibleBase(zMat a, zVec b, zVec base)
{
  _zLPTableau tab;
  bool ret = false;

  if( !zMatRowVecSizeIsEqual(a,b) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return false;
  }
  if( !_zLPTableauCreate( &tab, a, b ) ){
    ZALLOCERROR();
    return false;
  }
  if( !_zLPTableauSimplex( &tab ) || !zIsTiny(tab.d) ){
    ZRUNWARN( ZM_ERR_OPT_UNSOLVE );
  } else{
    _zLPTableauAns( &tab, base );
    ret = true;
  }
  _zLPTableauDestroy( &tab );
  return ret;
}
