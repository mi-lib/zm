/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_lp_simplex - optimization tools:
 * linear programming by simplex method.
 */

#include <zm/zm_opt.h>

#ifdef DEBUG
/* print out tableau contents to a file (for debug). */
static void _zLPTableauFPrint(FILE *fp, const zLPTableau *tab)
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
bool zLPTableauCreate(zLPTableau *tab, const zMat a, const zVec b)
{
  int i;

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
void zLPTableauDestroy(zLPTableau *tab)
{
  zMatFree( tab->a );
  zVecFree( tab->b );
  zVecFree( tab->c );
  zIndexFree( tab->ib );
  zIndexFree( tab->in );
  zIndexFree( tab->ir );
}

/* sweep-out one coefficient of cost function corresponding to old bases. */
static void _zLPTableauSweepC1(zLPTableau *tab, int p)
{
  int i;
  double cp;

  cp = zVecElemNC(tab->c,zIndexElemNC(tab->ib,p));
  for( i=0; i<zIndexSizeNC(tab->in); i++ )
    zVecElemNC(tab->c,zIndexElemNC(tab->in,i))
      -= zMatElemNC(tab->a,zIndexElemNC(tab->ir,p),zIndexElemNC(tab->in,i)) * cp;
  tab->d += zVecElemNC(tab->b,zIndexElemNC(tab->ir,p)) * cp;
  zVecSetElemNC( tab->c, zIndexElemNC(tab->ib,p), 0 );
}

/* sweep-out coefficients of cost function corresponding to bases. */
static void _zLPTableauSweepC(zLPTableau *tab)
{
  int i;

  for( i=0; i<zIndexSizeNC(tab->ib); i++ )
    _zLPTableauSweepC1( tab, i );
}

/* find next axis column of tableau to be sweeped-out. */
static int _zLPTableauFindNA(zLPTableau *tab, double ctol)
{ /* next axis for normal case */
  int i, na;
  double c_min;

  c_min = zVecElemNC( tab->c, zIndexElemNC(tab->in,(na=0)) );
  for( i=1; i<zIndexSizeNC(tab->in); i++ )
    if( zVecElemNC(tab->c,zIndexElemNC(tab->in,i)) < c_min )
      c_min = zVecElemNC( tab->c, zIndexElemNC(tab->in,(na=i)) );
  return c_min < -ctol ? na : -1;
}

/* find next axis column of tableau in degenerate case. */
static int _zLPTableauFindNA_deg(zLPTableau *tab, double ctol)
{ /* next axis for degenerated case */
  int i, na = -1, idx_min = INT_MAX;

  for( i=0; i<zIndexSizeNC(tab->in); i++ )
    if( zVecElemNC(tab->c,zIndexElemNC(tab->in,i)) < -ctol ){
      if( zIndexElemNC(tab->in,i) < idx_min ){
        idx_min = zIndexElemNC(tab->in,i);
        na = i;
      }
    }
  return na;
}

/* find next pivot in axis column to be sweeped-out.
 * if the step is degenerated (zero ratio), the axis column is chosen again by Bland's rule, and the pivot
 * is chosen among the rows with the zero ratio as the basic variable with the smallest index leaves the
 * bases, which prevents the iteration from cycling. */
static int _zLPTableauFindNP(zLPTableau *tab, int *na, double ctol)
{ /* next pivot */
  int i, np;
  double a, p, p_min, atol;
  bool bland = false;

 RETRY:
  /* a should be compared with 0 in theory, but an element which is small relative to the maximum
   * absolute value in the column is not chosen as the pivot, since it might be a round-off error. */
  for( atol=0, i=0; i<zIndexSizeNC(tab->ir); i++ )
    if( ( a = fabs( zMatElemNC(tab->a,zIndexElemNC(tab->ir,i),zIndexElemNC(tab->in,*na)) ) ) > atol ) atol = a;
  atol *= ZM_LP_SIMPLEX_PIVOT_TOL;
  for( p_min=HUGE_VAL, np=-1, i=0; i<zIndexSizeNC(tab->ir); i++ ){
    if( ( a = zMatElemNC(tab->a,zIndexElemNC(tab->ir,i),zIndexElemNC(tab->in,*na)) ) <= atol )
      continue;
    p = zVecElemNC(tab->b,zIndexElemNC(tab->ir,i)) / a;
    if( zIsTiny( p ) ){ /* degenerated case */
      if( !bland ){
        *na = _zLPTableauFindNA_deg( tab, ctol );
        bland = true;
        goto RETRY;
      }
      if( !zIsTiny( p_min ) || zIndexElemNC(tab->ib,i) < zIndexElemNC(tab->ib,np) ){
        p_min = p;
        np = i;
      }
      continue;
    }
    if( p < p_min ){
      p_min = p;
      np = i;
    }
  }
  return np;
}

/* sweep-out tabeau matrix of constraint equation. */
static void _zLPTableauSweepA(zLPTableau *tab, int np, int na)
{
  int i, j, r;
  double ap;

  /* normalize pivot row */
  ap = zMatElemNC(tab->a,zIndexElemNC(tab->ir,np),zIndexElemNC(tab->in,na));
  for( j=0; j<zIndexSizeNC(tab->in); j++ )
    zMatElemNC(tab->a,zIndexElemNC(tab->ir,np),zIndexElemNC(tab->in,j)) /= ap;
  zMatElemNC(tab->a,zIndexElemNC(tab->ir,np),zIndexElemNC(tab->ib,np)) /= ap;
  zVecElemNC(tab->b,zIndexElemNC(tab->ir,np)) /= ap;
  /* sweep-out rest row */
  for( i=0; i<zIndexSizeNC(tab->ir); i++ ){
    if( i == np ) continue;
    r = zIndexElemNC(tab->ir,i);
    ap = zMatElemNC(tab->a,r,zIndexElemNC(tab->in,na));
    for( j=0; j<zIndexSizeNC(tab->in); j++ )
      zMatElemNC(tab->a,r,zIndexElemNC(tab->in,j))
        -= zMatElemNC(tab->a,zIndexElemNC(tab->ir,np),zIndexElemNC(tab->in,j)) * ap;
    zMatElemNC(tab->a,r,zIndexElemNC(tab->ib,np))
      -= zMatElemNC(tab->a,zIndexElemNC(tab->ir,np),zIndexElemNC(tab->ib,np)) * ap;
    zVecElemNC(tab->b,r) -= zVecElemNC(tab->b,zIndexElemNC(tab->ir,np)) * ap;
  }
}

/* swap base/non-base pivot. */
static void _zLPTableauSwapPivot(zLPTableau *tab, int np, int na)
{
  int tmp;

  tmp = zIndexElemNC( tab->in, na );
  zIndexSetElemNC( tab->in, na, zIndexElemNC(tab->ib,np) );
  zIndexSetElemNC( tab->ib, np, tmp );
}

/* index from which artificial variables are numbered in the first stage. */
#define _zLPTableauArtIndex(tab) ( zVecSizeNC((tab)->c) - zVecSizeNC((tab)->b) )

/* number of artificial variables in the bases. */
static int _zLPTableauNumArt(const zLPTableau *tab)
{
  int i, n = 0;

  for( i=0; i<zIndexSizeNC(tab->ib); i++ )
    if( zIndexElemNC(tab->ib,i) >= _zLPTableauArtIndex(tab) ) n++;
  return n;
}

/* sum of values of artificial variables in the bases, which is the cost of the first stage.
 * It is computed directly from the bases instead of the cost accumulated through the
 * iteration, which suffers from round-off errors. */
static double _zLPTableauSumArtVal(const zLPTableau *tab)
{
  int i;
  double sum = 0;

  for( i=0; i<zIndexSizeNC(tab->ib); i++ )
    if( zIndexElemNC(tab->ib,i) >= _zLPTableauArtIndex(tab) )
      sum += zVecElemNC(tab->b,zIndexElemNC(tab->ir,i));
  return sum;
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
static bool _zLPTableauSimplex(zLPTableau *tab, double tol)
{ /* if tol >= 0 (the first stage), the iteration terminates as soon as all artificial variables
   * leave the bases or the cost becomes smaller than tol in the absolute value. */
  int iter = 0, na, np, n_art = 0;
  double ctol; /* tolerance for reduced costs relative to the given cost coefficients */

  ctol = ZM_LP_SIMPLEX_TOL * ( 1 + _zVecElemAbsMax( tab->c, NULL ) );
  if( tol >= 0 ) n_art = _zLPTableauNumArt( tab );
  _zLPTableauSweepC( tab );
  while( !( tol >= 0 && ( n_art == 0 || zIsTol( tab->d, tol ) ) ) &&
         ( na = _zLPTableauFindNA( tab, ctol ) ) >= 0 ){
    /* the iteration never cycles by Bland's rule in theory, but it is regarded as trapped if the number
     * of iterations exceeds five times the number of columns (the number of constraints plus the number
     * of variables), which would be due to round-off errors or a hard degenerated problem. */
    if( ++iter > 5 * zMatColSizeNC(tab->a) ){
      ZRUNWARN( ZM_ERR_OPT_TRAPPED );
      return false;
    }
    if( ( np = _zLPTableauFindNP( tab, &na, ctol ) ) < 0 ){ /* unbounded case */
      ZRUNERROR( ZM_ERR_OPT_INFINITESOLUTION );
      return false;
    }
    _zLPTableauSweepA( tab, np, na );
    if( tol >= 0 ){ /* count artificial variables in the bases */
      if( zIndexElemNC(tab->ib,np) >= _zLPTableauArtIndex(tab) ) n_art--;
      if( zIndexElemNC(tab->in,na) >= _zLPTableauArtIndex(tab) ) n_art++;
    }
    _zLPTableauSwapPivot( tab, np, na );
    _zLPTableauSweepC1( tab, np );
  }
  return true;
}
bool zLPTableauSimplex(zLPTableau *tab)
{
  return _zLPTableauSimplex( tab, -1 ); /* never terminates by the cost */
}

/* arrange the bases for the second stage after the first stage. */
static void _zLPTableauArrangeBase(zLPTableau *tab, double amax)
{
  int i, j, jmax, n, r;
  double a, binv;

  tab->d = 0;
  n = _zLPTableauArtIndex(tab);
  /* artificial variables remaining in the bases are replaced by non-artificial variables, or removed
   * together with the constraints if the constraints are redundant. */
  for( i=0; i<zIndexSizeNC(tab->ib); i++ ){
    if( zIndexElemNC(tab->ib,i) < n ) continue;
    r = zIndexElemNC(tab->ir,i);
    for( binv=0, j=n; j<zMatColSizeNC(tab->a); j++ )
      if( ( a = fabs( zMatElemNC(tab->a,r,j) ) ) > binv ) binv = a;
    for( jmax=-1, j=0; j<zIndexSizeNC(tab->in); j++ )
      if( zIndexElemNC(tab->in,j) < n &&
          ( a = fabs( zMatElemNC(tab->a,r,zIndexElemNC(tab->in,j)) ) ) > ZM_LP_SIMPLEX_PIVOT_TOL * binv * amax ){
        if( jmax < 0 || a > fabs( zMatElemNC(tab->a,r,zIndexElemNC(tab->in,jmax)) ) ) jmax = j;
      }
    if( jmax >= 0 ){
      /* the artificial variable, which vanishes within the tolerance, is replaced at zero level
       * so as not to make the other variables negative. */
      zVecSetElemNC( tab->b, r, 0 );
      _zLPTableauSweepA( tab, i, jmax );
      _zLPTableauSwapPivot( tab, i, jmax );
    } else{ /* redundant constraint */
      zIndexRemove( tab->ir, i );
      zIndexRemove( tab->ib, i );
      i--;
    }
  }
  /* remove artificial variables from non-bases */
  for( i=0; i<zIndexSizeNC(tab->in); i++ )
    if( zIndexElemNC(tab->in,i) >= n ){
      zIndexRemove( tab->in, i );
      i--;
    }
}

/* the first stage of simplex method to find a feasible base. */
bool zLPTableauFindFeasibleBase(zLPTableau *tab)
{
  int i, j;
  double scale = 1.0; /* initial sum of artificial variables (+1) */
  double amax = 0, a; /* maximum absolute value of the initial coefficients of non-artificial variables */

  for( i=0; i<zVecSizeNC(tab->b); i++ ){
    scale += fabs( zVecElemNC(tab->b,i) );
    for( j=0; j<_zLPTableauArtIndex(tab); j++ )
      if( ( a = fabs( zMatElemNC(tab->a,i,j) ) ) > amax ) amax = a;
  }
  /* the iteration terminates as soon as the sum of artificial variables vanishes, since
   * the first stage has a degenerate optimum, at which the round-off errors of reduced
   * costs would let the iteration go astray. */
  if( !_zLPTableauSimplex( tab, zTOL * scale ) ) return false;
  if( !zIsTol( _zLPTableauSumArtVal( tab ), zTOL * scale ) ){ /* infeasible case */
    ZRUNWARN( ZM_ERR_OPT_UNSOLVABLE );
    return false;
  }
  _zLPTableauArrangeBase( tab, amax );
  return true;
}

static void _zLTableauSetC(zLPTableau *tab, const zVec c)
{
  zVecSetSize( tab->c, zVecSizeNC(c) );
  zVecCopyNC( c, tab->c );
}

/* arrange answer vector. */
static void _zLPTableauAns(zLPTableau *tab, zVec ans)
{
  int i;

  zVecZero( ans );
  for( i=0; i<zIndexSizeNC(tab->ib); i++ )
    zVecSetElemNC( ans, zIndexElemNC(tab->ib,i), zVecElemNC(tab->b,zIndexElemNC(tab->ir,i)) );
}

/* dual-phase simplex method for linear programming. */
bool zLPSolveSimplex(const zMat a, const zVec b, const zVec c, zVec ans, double *cost)
{
  zLPTableau tab;
  bool ret = false;

  if( !zMatColVecSizeEqual(a,ans) || !zMatRowVecSizeEqual(a,b) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return false;
  }
  if( !zVecSizeEqual(c,ans) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return false;
  }
  if( !zLPTableauCreate( &tab, a, b ) ){
    ret = false;
    goto TERMINATE;
  }
  /* first phase: initial feasible base */
  if( !zLPTableauFindFeasibleBase( &tab ) ) goto TERMINATE;
  /* second phase: body problem */
  _zLTableauSetC( &tab, c );
  if( !zLPTableauSimplex( &tab ) ) goto TERMINATE;
  if( cost ) *cost = tab.d;
  _zLPTableauAns( &tab, ans );
  ret = true;
 TERMINATE:
  zLPTableauDestroy( &tab );
  return ret;
}

/* find a feasible base under Ax=b and x>=0 based on simplex method. */
bool zLPFeasibleBase(const zMat a, const zVec b, zVec base)
{
  zLPTableau tab;
  bool ret = false;

  if( !zMatRowVecSizeEqual(a,b) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return false;
  }
  if( !zLPTableauCreate( &tab, a, b ) ){
    ZALLOCERROR();
    return false;
  }
  if( zLPTableauFindFeasibleBase( &tab ) ){
    _zLPTableauAns( &tab, base );
    ret = true;
  }
  zLPTableauDestroy( &tab );
  return ret;
}
