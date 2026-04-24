/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_lp_megiddodyer - optimization tools: linear programming based on Megiddo-Dyer algorithm.
 */

#include <zm/zm_opt.h>

/* a linear inequality constraint. */

static double _zLPMegiddoDyerConstraintVal(const zLPMegiddoDyerConstraint *constraint, double x)
{
  return constraint->m * x + constraint->y0;
}

#if DEBUG == 1
static double _zLPMegiddoDyerConstraintInvVal(const zLPMegiddoDyerConstraint *constraint, double y)
{
  return zIsTiny( constraint->m ) ? NAN : ( y - constraint->y0 ) / constraint->m;
}
#endif /* DEBUG */

static bool _zLPMegiddoDyerConstraintIntersection(const zLPMegiddoDyerConstraint *c1, const zLPMegiddoDyerConstraint *c2, double *x, double *y)
{
  double den;

  if( zIsTiny( ( den = c1->m - c2->m ) ) ) return false;
  *x = -( c1->y0 - c2->y0 ) / den;
  *y = ( c1->m * c2->y0 - c2->m * c1->y0 ) / den;
  return true;
}

/* list of inequality constraints. */

static double _zLPMegiddoDyerConstraintListVal(const zLPMegiddoDyerConstraintList *list, double x, bool (* cmp)(double,double), double val_init)
{
  zLPMegiddoDyerConstraintListCell *cp;
  double val, tmp;

  val = val_init;
  zListForEach( list, cp )
    if( cmp( ( tmp = _zLPMegiddoDyerConstraintVal( &cp->data, x ) ), val ) ) val = tmp;
  return val;
}

static double _zLPMegiddoDyerConstraintListValGrad(const zLPMegiddoDyerConstraintList *list, double x, bool (* cmp)(double,double), double val_init, double *m_min, double *m_max)
{
  zLPMegiddoDyerConstraintListCell *cp;
  double val, tmp;

  *m_min = HUGE_VAL;
  *m_max =-HUGE_VAL;
  val = val_init;
  zListForEach( list, cp ){
    tmp = _zLPMegiddoDyerConstraintVal( &cp->data, x );
    if( zEqual( tmp, val, zTOL ) ){
      if( cp->data.m > *m_max ) *m_max = cp->data.m;
      if( cp->data.m < *m_min ) *m_min = cp->data.m;
    } else
    if( cmp( tmp, val ) ){
      val = tmp;
      *m_min = *m_max = cp->data.m;
    }
  }
  return val;
}

/* list of pairs of linear inequality constraints. */

static bool _zLPMegiddoDyerConstraintPairListCreate(zLPMegiddoDyerConstraintPairList *pairlist, zLPMegiddoDyerConstraintPtrList *constraintpool, zLPMegiddoDyerConstraintList *constraintlist, bool (* discard_rule)(double,double))
{
  zLPMegiddoDyerConstraintPtrListCell *cpp1, *cpp2;
  zLPMegiddoDyerConstraintPairListCell *pp;
  double x, y;

  while( zListSize(constraintpool) >= 2 ){
    zListDeleteHead( constraintpool, &cpp1 );
    zListDeleteHead( constraintpool, &cpp2 );
    if( !_zLPMegiddoDyerConstraintIntersection( &cpp1->data->data, &cpp2->data->data, &x, &y ) ){
      if( !discard_rule( cpp1->data->data.y0, cpp2->data->data.y0 ) )
        zSwap( zLPMegiddoDyerConstraintPtrListCell*, cpp1, cpp2 ); /* cpp1 -> to discard, cpp2 -> to revert */
      zListPurge( constraintlist, cpp1->data );
      free( cpp1->data );
      free( cpp1 );
      zListInsertHead( constraintpool, cpp2 );
      continue;
    }
    if( !( pp = zAlloc( zLPMegiddoDyerConstraintPairListCell, 1 ) ) ){
      ZALLOCERROR();
      /* push back cells in order to avoid memory-leakage */
      zListInsertHead( constraintpool, cpp1 );
      zListInsertHead( constraintpool, cpp2 );
      return false;
    }
    pp->data.cp1 = cpp1;
    pp->data.cp2 = cpp2;
    pp->data.x = x;
    pp->data.y = y;
    zListInsertHead( pairlist, pp );
  }
  return zListSize(pairlist);
}

ZEDA_DEF_LIST_QUICKSORT( zLPMegiddoDyerConstraintPairList, zLPMegiddoDyerConstraintPairListCell )

static int _zLPMegiddoDyerConstraintPairCompare(zLPMegiddoDyerConstraintPairListCell *c1, zLPMegiddoDyerConstraintPairListCell *c2, void *dummy)
{
  if( c1->data.x > c2->data.x ) return 1;
  if( c1->data.x < c2->data.x ) return -1;
  return 0;
}

static zLPMegiddoDyerConstraintPairListCell *_zLPMegiddoDyerConstraintPairListMedian(zLPMegiddoDyerConstraintPairList *pairlist)
{
  return zLPMegiddoDyerConstraintPairListQuickMedian( pairlist, _zLPMegiddoDyerConstraintPairCompare, NULL );
}

#define ZM_LP_MEGIDDO_DYER_X_POSITIVE_BOUNDED 0x01
#define ZM_LP_MEGIDDO_DYER_X_NEGATIVE_BOUNDED 0x02
#define ZM_LP_MEGIDDO_DYER_DETERMINATE        0x04

#define ZM_LP_MEGIDDO_DYER_X_BOUNDED ( ZM_LP_MEGIDDO_DYER_X_POSITIVE_BOUNDED | ZM_LP_MEGIDDO_DYER_X_NEGATIVE_BOUNDED )

/* initialize properties for Megiddo-Dyer algorithm. */
void zLPMegiddoDyerInit(zLPMegiddoDyer *lp)
{
  zListInit( &lp->upperbound );
  zListInit( &lp->lowerbound );
  zListInit( &lp->upperboundpool );
  zListInit( &lp->lowerboundpool );
  zListInit( &lp->upperboundpair );
  zListInit( &lp->lowerboundpair );
  lp->cost_c1 = lp->cost_c2 = 0;
  lp->x_min =-HUGE_VAL;
  lp->x_max = HUGE_VAL;
  lp->mu_min = lp->ml_min = HUGE_VAL;
  lp->mu_max = lp->ml_max =-HUGE_VAL;
  lp->feasibility_flag = 0;
}

/* destroy properties for Megiddo-Dyer algorithm. */
void zLPMegiddoDyerDestroy(zLPMegiddoDyer *lp)
{
  zListDestroy( zLPMegiddoDyerConstraintListCell, &lp->upperbound );
  zListDestroy( zLPMegiddoDyerConstraintListCell, &lp->lowerbound );
  zListDestroy( zLPMegiddoDyerConstraintPtrListCell, &lp->upperboundpool );
  zListDestroy( zLPMegiddoDyerConstraintPtrListCell, &lp->lowerboundpool );
  zListDestroy( zLPMegiddoDyerConstraintPairListCell, &lp->upperboundpair );
  zListDestroy( zLPMegiddoDyerConstraintPairListCell, &lp->lowerboundpair );
}

/* set coefficients of cost function for Megiddo-Dyer algorithm. */
bool zLPMegiddoDyerSetCostCoefficient(zLPMegiddoDyer *lp, double c1, double c2)
{
  if( zIsTiny( ( lp->cost_c2 = c2 ) ) ){
    ZRUNERROR( ZM_ERR_OPT_LP2_INVALID_COSTFUNC );
    return false;
  }
  lp->cost_c1 = c1;
  return true;
}

static double _zLPMegiddoDyerXYToX2(zLPMegiddoDyer *lp, double x, double y)
{
  return ( y - lp->cost_c1 * x ) / lp->cost_c2;
}

/* add a linear inequality constraint to a linear programming problem to be solved by Megiddo-Dyer algorithm. */
bool zLPMegiddoDyerAddConstraint(zLPMegiddoDyer *lp, double a1, double a2, double b)
{
  zLPMegiddoDyerConstraintListCell *cp;
  zLPMegiddoDyerConstraintPtrListCell *cpp;
  double x_tmp;

  if( zIsTiny( a2 ) ){
    if( zIsTiny( a1 ) ){
      ZRUNWARN( ZM_ERR_OPT_LP2_INVALID_CONSTRAINT );
      return false;
    }
    x_tmp = b / a1;
    if( a1 < 0 && x_tmp > lp->x_min ){
      if( ( lp->x_min = x_tmp ) > lp->x_max ){
        ZRUNERROR( ZM_ERR_OPT_UNSOLVABLE );
        return false;
      }
    } else
    if( a1 > 0 && x_tmp < lp->x_max ){
      if( ( lp->x_max = x_tmp ) < lp->x_min ){
        ZRUNERROR( ZM_ERR_OPT_UNSOLVABLE );
        return false;
      }
    }
    lp->feasibility_flag |= ZM_LP_MEGIDDO_DYER_DETERMINATE;
    return true;
  }
  if( zIsTiny( lp->cost_c2 ) ){
    ZRUNERROR( ZM_ERR_OPT_LP2_INVALID_COSTFUNC );
    return false;
  }
  cp  = zAlloc( zLPMegiddoDyerConstraintListCell, 1 );
  cpp = zAlloc( zLPMegiddoDyerConstraintPtrListCell, 1 );
  if( !cp || !cpp ){
    ZALLOCERROR();
    free( cp );
    free( cpp );
    return false;
  }
  cpp->data = cp;
  a2 /= lp->cost_c2;
  if( !zIsTiny( cp->data.m = lp->cost_c1 - a1 / a2 ) )
    lp->feasibility_flag |= ZM_LP_MEGIDDO_DYER_DETERMINATE;
  cp->data.y0 = b / a2;
  if( a2 > 0 ){
    zListInsertHead( &lp->upperbound, cp );
    zListInsertHead( &lp->upperboundpool, cpp );
    if( cp->data.m > lp->mu_max ) lp->mu_max = cp->data.m;
    if( cp->data.m < lp->mu_min ) lp->mu_min = cp->data.m;
  } else{
    zListInsertHead( &lp->lowerbound, cp );
    zListInsertHead( &lp->lowerboundpool, cpp );
    if( cp->data.m > lp->ml_max ) lp->ml_max = cp->data.m;
    if( cp->data.m < lp->ml_min ) lp->ml_min = cp->data.m;
  }
  /* check if the feasible solution is bounded */
  if( lp->ml_max > 0 ||
      ( !zListIsEmpty( &lp->lowerbound ) && lp->x_max < HUGE_VAL ) ||
      ( !zListIsEmpty( &lp->upperbound ) && !zListIsEmpty( &lp->lowerbound ) && lp->mu_min < lp->ml_max ) )
    lp->feasibility_flag |= ZM_LP_MEGIDDO_DYER_X_POSITIVE_BOUNDED;
  if( lp->ml_min < 0 ||
      ( !zListIsEmpty( &lp->lowerbound ) && lp->x_min >-HUGE_VAL ) ||
      ( !zListIsEmpty( &lp->upperbound ) && !zListIsEmpty( &lp->lowerbound ) && lp->mu_max > lp->ml_min ) )
    lp->feasibility_flag |= ZM_LP_MEGIDDO_DYER_X_NEGATIVE_BOUNDED;
  return true;
}

static double _zLPMegiddoDyerUpperboundVal(const zLPMegiddoDyer *lp, double x)
{
  return _zLPMegiddoDyerConstraintListVal( &lp->upperbound, x, zIsLess, HUGE_VAL );
}
static double _zLPMegiddoDyerLowerboundVal(const zLPMegiddoDyer *lp, double x )
{
  return _zLPMegiddoDyerConstraintListVal( &lp->lowerbound, x, zIsGreater,-HUGE_VAL );
}
static double _zLPMegiddoDyerUpperboundValGrad(const zLPMegiddoDyer *lp, double x, double *m_min, double *m_max)
{
  return _zLPMegiddoDyerConstraintListValGrad( &lp->upperbound, x, zIsLess, HUGE_VAL, m_min, m_max );
}
static double _zLPMegiddoDyerLowerboundValGrad(const zLPMegiddoDyer *lp, double x, double *m_min, double *m_max)
{
  return _zLPMegiddoDyerConstraintListValGrad( &lp->lowerbound, x, zIsGreater,-HUGE_VAL, m_min, m_max );
}

static bool _zLPMegiddoDyerCreateUpperboundPairList(zLPMegiddoDyer *lp)
{
  return _zLPMegiddoDyerConstraintPairListCreate( &lp->upperboundpair, &lp->upperboundpool, &lp->upperbound, zIsGreater );
}
static bool _zLPMegiddoDyerCreateLowerboundPairList(zLPMegiddoDyer *lp)
{
  return _zLPMegiddoDyerConstraintPairListCreate( &lp->lowerboundpair, &lp->lowerboundpool, &lp->lowerbound, zIsLess );
}

#if DEBUG == 1 /* -v-v-v- for debug -v-v-v- */

static void _zLPMegiddoDyerOutputBoundaryLine(FILE *fp, zLPMegiddoDyerConstraint *constraint, double x_left, double x_right, double y_bottom, double y_top)
{
  double x1, y1, x2, y2;

  y1 = _zLPMegiddoDyerConstraintVal( constraint, ( x1 = x_left ) );
  y2 = _zLPMegiddoDyerConstraintVal( constraint, ( x2 = x_right ) );
  if( y1 < y_bottom && y2 > y_bottom ) x1 = _zLPMegiddoDyerConstraintInvVal( constraint, ( y1 = y_bottom ) );
  if( y1 > y_top    && y2 < y_top    ) x1 = _zLPMegiddoDyerConstraintInvVal( constraint, ( y1 = y_top    ) );
  if( y1 < y_top    && y2 > y_top    ) x2 = _zLPMegiddoDyerConstraintInvVal( constraint, ( y2 = y_top    ) );
  if( y1 > y_bottom && y2 < y_bottom ) x2 = _zLPMegiddoDyerConstraintInvVal( constraint, ( y2 = y_bottom ) );
  fprintf( fp, "%g %g\n%g %g\n\n", x1, y1, x2, y2 );
}

static void _zLPMegiddoDyerOutputBoundarySet(zLPMegiddoDyerConstraintList *list, double x_left, double x_right, double y_bottom, double y_top, const char *filename)
{
  zLPMegiddoDyerConstraintListCell *cp;
  FILE *fp;

  if( zListIsEmpty( list ) ) return;
  fp = fopen( filename, "w" );
  zListForEach( list, cp )
    _zLPMegiddoDyerOutputBoundaryLine( fp, &cp->data, x_left, x_right, y_bottom, y_top );
  fclose( fp );
}

static void _zLPMegiddoDyerOutputBoundaryFunc(zLPMegiddoDyerConstraintList *list, double x_left, double x_right, int div, bool (*cmp)(double,double), double val_init, const char *filename)
{
  FILE *fp;
  double x;
  int i;

  if( zListIsEmpty( list ) ) return;
  fp = fopen( filename, "w" );
  for( i=0; i<div; i++ ){
    x = x_left + ( x_right - x_left ) * i / div;
    fprintf( fp, "%g %g\n", x, _zLPMegiddoDyerConstraintListVal( list, x, cmp, val_init ) );
  }
  fclose( fp );
}

static void _zLPMegiddoDyerOutputBoundary(zLPMegiddoDyer *lp, double x_left, double x_right, double y_bottom, double y_top)
{
  const int div = 100;

  _zLPMegiddoDyerOutputBoundarySet( &lp->upperbound, x_left, x_right, y_bottom, y_top, "us" );
  _zLPMegiddoDyerOutputBoundaryFunc( &lp->upperbound, x_left, x_right, div, zIsLess, HUGE_VAL, "ub" );
  _zLPMegiddoDyerOutputBoundarySet( &lp->lowerbound, x_left, x_right, y_bottom, y_top, "ls" );
  _zLPMegiddoDyerOutputBoundaryFunc( &lp->lowerbound, x_left, x_right, div, zIsGreater,-HUGE_VAL, "lb" );
}

static void _zLPMegiddoDyerOutputBoundaryPair(zLPMegiddoDyerConstraintPairList *pairlist, double x_left, double x_right, double y_bottom, double y_top, const char *filename)
{
  zLPMegiddoDyerConstraintPairListCell *cp;
  FILE *fp;

  if( zListIsEmpty( pairlist ) ) return;
  fp = fopen( filename, "w" );
  zListForEach( pairlist, cp ){
    _zLPMegiddoDyerOutputBoundaryLine( fp, &cp->data.cp1->data->data, x_left, x_right, y_bottom, y_top );
    _zLPMegiddoDyerOutputBoundaryLine( fp, &cp->data.cp2->data->data, x_left, x_right, y_bottom, y_top );
  }
  fclose( fp );
}

static void _zLPMegiddoDyerOutputIntersection(zLPMegiddoDyerConstraintPairList *pairlist, const char *filename)
{
  zLPMegiddoDyerConstraintPairListCell *cp;
  FILE *fp;

  if( zListIsEmpty( pairlist ) ) return;
  fp = fopen( filename, "w" );
  zListForEach( pairlist, cp )
    fprintf( fp, "%g %g\n", cp->data.x, cp->data.y );
  fclose( fp );
}

static void _zLPMegiddoDyerOutputPoint(double x, double y, const char *filename)
{
  FILE *fp;

  fp = fopen( filename, "w" );
  fprintf( fp, "%g %g\n", x, y );
  fclose( fp );
}

static void _zLPMegiddoDyerOutputValGrad(double x, double y, double m_min, double m_max, double x_left, double x_right, double y_bottom, double y_top, const char *filename)
{
  FILE *fp;
  zLPMegiddoDyerConstraint c;

  fp = fopen( filename, "w" );
  c.m  = m_min;
  c.y0 = y - m_min * x;
  _zLPMegiddoDyerOutputBoundaryLine( fp, &c, x_left, x_right, y_bottom, y_top );
  c.m  = m_max;
  c.y0 = y - m_max * x;
  _zLPMegiddoDyerOutputBoundaryLine( fp, &c, x_left, x_right, y_bottom, y_top );
  fclose( fp );
}

#endif /* DEBUG */ /* -^-^-^- for debug -^-^-^- */

static bool _zLPMegiddoDyerXYIsFeasible(const zLPMegiddoDyer *lp, double x, double y)
{
  return x >= lp->x_min - zTOL && x <= lp->x_max + zTOL &&
         y <= _zLPMegiddoDyerUpperboundVal( lp, x ) + zTOL &&
         y >= _zLPMegiddoDyerLowerboundVal( lp, x ) - zTOL;
}

static bool _zLPMegiddoDyerTestAndUpdateDirect(const zLPMegiddoDyer *lp, double *x, double *y, double test_x, double test_y)
{
  if( !_zLPMegiddoDyerXYIsFeasible( lp, test_x, test_y ) || test_y > *y ) return false;
  *x = test_x;
  *y = test_y;
  return true;
}

static bool _zLPMegiddoDyerSolveDirect(const zLPMegiddoDyer *lp, double *x, double *y)
{
  zLPMegiddoDyerConstraintListCell *cp1, *cp2;
  double _x, _y;
  bool has_solution = false;

  *x = NAN;
  *y = HUGE_VAL;
  /* upperbound vs lowerbound */
  zListForEach( &lp->upperbound, cp1 ){
    zListForEach( &lp->lowerbound, cp2 ){
      if( !_zLPMegiddoDyerConstraintIntersection( &cp1->data, &cp2->data, &_x, &_y ) ) continue;
      if( _zLPMegiddoDyerTestAndUpdateDirect( lp, x, y, _x, _y ) ){
        has_solution = true;
      }
    }
  }
  /* lowerbound vs lowerbound */
  zListForEach( &lp->lowerbound, cp1 ){
    for( cp2=zListCellNext(cp1); cp2!=zListRoot(&lp->lowerbound); cp2=zListCellNext(cp2) ){
      if( !_zLPMegiddoDyerConstraintIntersection( &cp1->data, &cp2->data, &_x, &_y ) ) continue;
      if( _zLPMegiddoDyerTestAndUpdateDirect( lp, x, y, _x, _y ) ){
        has_solution = true;
      }
    }
  }
  /* minimum x-boundary */
  if( ( _x = lp->x_min ) != -HUGE_VAL ){
    _y = _zLPMegiddoDyerLowerboundVal( lp, _x );
    if( _zLPMegiddoDyerTestAndUpdateDirect( lp, x, y, _x, _y ) ){
      has_solution = true;
    }
  }
  /* maximum x-boundary */
  if( ( _x = lp->x_max ) != HUGE_VAL ){
    _y = _zLPMegiddoDyerLowerboundVal( lp, _x );
    if( _zLPMegiddoDyerTestAndUpdateDirect( lp, x, y, _x, _y ) ){
      has_solution = true;
    }
  }
  return has_solution;
}

typedef enum{
  ZM_LP_MEGIDDO_DYER_INVALID = -1,
  ZM_LP_MEGIDDO_DYER_UPDATE_XMIN,
  ZM_LP_MEGIDDO_DYER_UPDATE_XMAX,
  ZM_LP_MEGIDDO_DYER_UNSOLVABLE,
  ZM_LP_MEGIDDO_DYER_FOUND_OPTIMUM,
} zLPMegiddoDyerUpdateCode;

static zLPMegiddoDyerUpdateCode _zLPMegiddoDyerUpdateXMin(zLPMegiddoDyer *lp, double xm)
{
  lp->x_min = xm;
  return ZM_LP_MEGIDDO_DYER_UPDATE_XMIN;
}

static zLPMegiddoDyerUpdateCode _zLPMegiddoDyerUpdateXMax(zLPMegiddoDyer *lp, double xm)
{
  lp->x_max = xm;
  return ZM_LP_MEGIDDO_DYER_UPDATE_XMAX;
}

static zLPMegiddoDyerUpdateCode _zLPMegiddoDyerFoundOptimum(zLPMegiddoDyer *lp, double *x, double *y, double xo, double yo)
{
  *x = xo;
  *y = yo;
  return ZM_LP_MEGIDDO_DYER_FOUND_OPTIMUM;
}

static zLPMegiddoDyerUpdateCode _zLPMegiddoDyerUpdateXRegion(zLPMegiddoDyer *lp, double *x, double *y)
{
  zLPMegiddoDyerConstraintPairListCell *median;
  double yu, yl, mu_min, mu_max, ml_min, ml_max;

  _zLPMegiddoDyerCreateUpperboundPairList( lp );
  _zLPMegiddoDyerCreateLowerboundPairList( lp );
  if( !( median = _zLPMegiddoDyerConstraintPairListMedian( &lp->lowerboundpair ) ) ){
    if( !( median = _zLPMegiddoDyerConstraintPairListMedian( &lp->upperboundpair ) ) ){
      ZRUNWARN( ZM_ERR_OPT_LP2_NO_CONSTRAINT_PAIR );
      return ZM_LP_MEGIDDO_DYER_INVALID;
    }
  }
  yu = _zLPMegiddoDyerUpperboundValGrad( lp, median->data.x, &mu_min, &mu_max );
  yl = _zLPMegiddoDyerLowerboundValGrad( lp, median->data.x, &ml_min, &ml_max );
  if( yl > yu ){
    if( ml_max < mu_min ){
      return _zLPMegiddoDyerUpdateXMin( lp, median->data.x );
    } else
    if( ml_min > mu_max ){
      return _zLPMegiddoDyerUpdateXMax( lp, median->data.x );
    }
    return ZM_LP_MEGIDDO_DYER_UNSOLVABLE;
  } else
  if( yl < yu ){
    if( ml_max < 0 ){
      return _zLPMegiddoDyerUpdateXMin( lp, median->data.x );
    } else
    if( ml_min > 0 ){
      return _zLPMegiddoDyerUpdateXMax( lp, median->data.x );
    }
    return _zLPMegiddoDyerFoundOptimum( lp, x, y, median->data.x, yl );
  }
  if( ml_max < 0 && ml_max <= mu_min ){
    return _zLPMegiddoDyerUpdateXMin( lp, median->data.x );
  } else
  if( ml_min > 0 && ml_min >= mu_max ){
    return _zLPMegiddoDyerUpdateXMax( lp, median->data.x );
  }
  return _zLPMegiddoDyerFoundOptimum( lp, x, y, median->data.x, yl );
}

static void _zLPMegiddoDyerConstraintPairDiscard(zLPMegiddoDyerConstraintList *list, zLPMegiddoDyerConstraintPtrList *pool, zLPMegiddoDyerConstraintPtrListCell *to_discard, zLPMegiddoDyerConstraintPtrListCell *to_revert)
{
  zListPurge( list, to_discard->data );
  free( to_discard->data );
  free( to_discard );
  zListInsertHead( pool, to_revert );
}

static void _zLPMegiddoDyerDiscardConstraint(zLPMegiddoDyer *lp, double xm, bool (* cmp_x)(double,double), bool (* cmp_grad1)(double,double), bool (* cmp_grad2)(double,double))
{
  zLPMegiddoDyerConstraintPairListCell *cp, *cp_tmp;

  zListForEach( &lp->upperboundpair, cp ){
    if( cmp_x( cp->data.x, xm ) ) continue;
    if( cmp_grad1( cp->data.cp1->data->data.m, cp->data.cp2->data->data.m ) )
     _zLPMegiddoDyerConstraintPairDiscard( &lp->upperbound, &lp->upperboundpool, cp->data.cp1, cp->data.cp2 );
    else
     _zLPMegiddoDyerConstraintPairDiscard( &lp->upperbound, &lp->upperboundpool, cp->data.cp2, cp->data.cp1 );
    cp = zListCellPrev(cp);
    zListDeleteNext( &lp->upperboundpair, cp, &cp_tmp );
    free( cp_tmp );
  }
  zListForEach( &lp->lowerboundpair, cp ){
    if( cmp_x( cp->data.x, xm ) ) continue;
    if( cmp_grad2( cp->data.cp1->data->data.m, cp->data.cp2->data->data.m ) )
     _zLPMegiddoDyerConstraintPairDiscard( &lp->lowerbound, &lp->lowerboundpool, cp->data.cp1, cp->data.cp2 );
    else
     _zLPMegiddoDyerConstraintPairDiscard( &lp->lowerbound, &lp->lowerboundpool, cp->data.cp2, cp->data.cp1 );
    cp = zListCellPrev(cp);
    zListDeleteNext( &lp->lowerboundpair, cp, &cp_tmp );
    free( cp_tmp );
  }
}

static void _zLPMegiddoDyerDiscardConstraintSmallerX(zLPMegiddoDyer *lp)
{
  _zLPMegiddoDyerDiscardConstraint( lp, lp->x_min, zIsGreater, zIsGreater, zIsLess );
}

static void _zLPMegiddoDyerDiscardConstraintGreaterX(zLPMegiddoDyer *lp)
{
  _zLPMegiddoDyerDiscardConstraint( lp, lp->x_max, zIsLess, zIsLess, zIsGreater );
}

static bool _zLPMegiddoDyerCheckFeasibility(zLPMegiddoDyer *lp)
{
  if( ( lp->feasibility_flag & ZM_LP_MEGIDDO_DYER_X_BOUNDED ) != ZM_LP_MEGIDDO_DYER_X_BOUNDED ){
    ZRUNERROR( ZM_ERR_OPT_INFINITESOLUTION );
    return false;
  }
  if( !( lp->feasibility_flag & ZM_LP_MEGIDDO_DYER_DETERMINATE ) ){
    ZRUNERROR( ZM_ERR_OPT_NONUNIQUE );
    return false;
  }
  return true;
}

/* solve a linear programming problem by Megiddo-Dyer algorithm. */
bool zLPMegiddoDyerSolve(zLPMegiddoDyer *lp, double *x1, double *x2, double *cost)
{
  double _cost;

  *x1 = *x2 = NAN;
  if( !_zLPMegiddoDyerCheckFeasibility( lp ) ) return false;
  if( !cost ) cost = &_cost;
  while( zListSize(&lp->upperbound) + zListSize( &lp->lowerbound ) > 4 ){
    switch( _zLPMegiddoDyerUpdateXRegion( lp, x1, cost ) ){
    case ZM_LP_MEGIDDO_DYER_FOUND_OPTIMUM:
      goto ZM_LP_MEGIDDO_DYER_FOUND_SOLUTION;
    case ZM_LP_MEGIDDO_DYER_UPDATE_XMIN:
      _zLPMegiddoDyerDiscardConstraintSmallerX( lp );
      break;
    case ZM_LP_MEGIDDO_DYER_UPDATE_XMAX:
      _zLPMegiddoDyerDiscardConstraintGreaterX( lp );
      break;
    case ZM_LP_MEGIDDO_DYER_UNSOLVABLE:
      ZRUNERROR( ZM_ERR_OPT_UNSOLVABLE );
      return false;
    default: /* ZM_LP_MEGIDDO_DYER_INVALID */
      goto ZM_LP_MEGIDDO_DYER_SOLVE_DIRECTLY;
    }
  }
 ZM_LP_MEGIDDO_DYER_SOLVE_DIRECTLY:
  if( !_zLPMegiddoDyerSolveDirect( lp, x1, cost ) ) return false;
 ZM_LP_MEGIDDO_DYER_FOUND_SOLUTION:
  *x2 = _zLPMegiddoDyerXYToX2( lp, *x1, *cost );
  return true;
}

/* solve a linear programming problem directly by checking all intersections of boundaries of linear inequality constraints. */
bool zLPMegiddoDyerSolveDirect(zLPMegiddoDyer *lp, double *x1, double *x2, double *cost)
{
  double _cost;

  *x1 = *x2 = NAN;
  if( !_zLPMegiddoDyerCheckFeasibility( lp ) ) return false;
  if( !cost ) cost = &_cost;
  if( !_zLPMegiddoDyerSolveDirect( lp, x1, cost ) ){
    ZRUNERROR( ZM_ERR_OPT_UNSOLVABLE );
    return false;
  }
  *x2 = _zLPMegiddoDyerXYToX2( lp, *x1, *cost );
  return true;
}

#if DEBUG == 1
/* read a linear programming problem from a file for Megiddo-Dyer algorithm. */
zLPMegiddoDyer *zLPMegiddoDyerReadFile(zLPMegiddoDyer *lp, const char *filename)
{
  FILE *fp;
  int nc, i;
  double a1, a2, b;
  zLPMegiddoDyer *retval = NULL;

  if( !( fp = fopen( filename, "rb" ) ) ){
    ZOPENERROR( filename );
    return NULL;
  }
  zLPMegiddoDyerInit( lp );
  if( fread( &lp->cost_c1, sizeof(double), 1, fp ) != 1 ) goto ZM_LP_MEGIDDO_DYER_READ_FILE_TERMINATE;
  if( fread( &lp->cost_c2, sizeof(double), 1, fp ) != 1 ) goto ZM_LP_MEGIDDO_DYER_READ_FILE_TERMINATE;
  if( fread( &nc, sizeof(int), 1, fp ) != 1 ) goto ZM_LP_MEGIDDO_DYER_READ_FILE_TERMINATE;
  for( i=0; i<nc; i++ ){
    if( fread( &a1, sizeof(double), 1, fp ) != 1 ) goto ZM_LP_MEGIDDO_DYER_READ_FILE_TERMINATE;
    if( fread( &a2, sizeof(double), 1, fp ) != 1 ) goto ZM_LP_MEGIDDO_DYER_READ_FILE_TERMINATE;
    if( fread( &b,  sizeof(double), 1, fp ) != 1 ) goto ZM_LP_MEGIDDO_DYER_READ_FILE_TERMINATE;
    zLPMegiddoDyerAddConstraint( lp, a1, a2, b );
  }
  retval = lp;
 ZM_LP_MEGIDDO_DYER_READ_FILE_TERMINATE:
  fclose( fp );
  return retval;
}

/* write a linear programming problem to a file for Megiddo-Dyer algorithm. */
zLPMegiddoDyer *zLPMegiddoDyerWriteFile(zLPMegiddoDyer *lp, const char *filename)
{
  FILE *fp;
  int nc;
  double a1, a2, b;
  zLPMegiddoDyerConstraintListCell *cp;

  if( !( fp = fopen( filename, "wb" ) ) ){
    ZOPENERROR( filename );
    return NULL;
  }
  fwrite( &lp->cost_c1, sizeof(double), 1, fp );
  fwrite( &lp->cost_c2, sizeof(double), 1, fp );
  nc = zListSize(&lp->upperbound) + zListSize(&lp->lowerbound);
  fwrite( &nc, sizeof(int), 1, fp );
  zListForEach( &lp->upperbound, cp ){
    a1 =-cp->data.m;
    a2 = 1.0;
    b  = cp->data.y0;
    fwrite( &a1, sizeof(double), 1, fp );
    fwrite( &a2, sizeof(double), 1, fp );
    fwrite( &b,  sizeof(double), 1, fp );
  }
  zListForEach( &lp->lowerbound, cp ){
    a1 = cp->data.m;
    a2 =-1.0;
    b  =-cp->data.y0;
    fwrite( &a1, sizeof(double), 1, fp );
    fwrite( &a2, sizeof(double), 1, fp );
    fwrite( &b,  sizeof(double), 1, fp );
  }
  fclose( fp );
  return lp;
}

/* write a linear programming problem to a file in an ASCII format for Megiddo-Dyer algorithm. */
zLPMegiddoDyer *zLPMegiddoDyerWriteFileASCII(zLPMegiddoDyer *lp, const char *filename)
{
  FILE *fp;
  zLPMegiddoDyerConstraintListCell *cp;

  if( !( fp = fopen( filename, "w" ) ) ){
    ZOPENERROR( filename );
    return NULL;
  }
  fprintf( fp, "%g %g\n", lp->cost_c1, lp->cost_c2 );
  fprintf( fp, "%d\n", zListSize(&lp->upperbound) + zListSize(&lp->lowerbound) );
  zListForEach( &lp->upperbound, cp ){
    fprintf( fp, "%g 1.0 %g\n", -cp->data.m, cp->data.y0 );
  }
  zListForEach( &lp->upperbound, cp ){
    fprintf( fp, "%g -1.0 %g\n", cp->data.m,-cp->data.y0 );
  }
  fclose( fp );
  return lp;
}
#endif /* DEBUG */
