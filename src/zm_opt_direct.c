/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_direct - optimization tools: DIRECT method.
 */
#include <zm/zm_opt.h>

/* hyper-rectangle region */
typedef struct{
  zVec x;
  double val;
  int *level;
} zOptDIRECTRect;

/* free memory for a hyper-rectangle region */
static zOptDIRECTRect *_zOptDIRECTRectFree(zOptDIRECTRect *rect)
{
  zVecFree( rect->x );
  free( rect->level );
  free( rect );
  return NULL;
}

/* allocate memory for a hyper-rectangle region */
static zOptDIRECTRect *_zOptDIRECTRectAlloc(int dim)
{
  zOptDIRECTRect *rect;

  if( !( rect = zAlloc( zOptDIRECTRect, 1 ) ) ) return NULL;
  rect->x = zVecAlloc( dim );
  rect->level = zAlloc( int, dim );
  if( !rect->x || !rect->level ) return _zOptDIRECTRectFree( rect );
  return rect;
}

/* clone a hyper-rectangle region */
static zOptDIRECTRect *_zOptDIRECTRectClone(zOptDIRECTRect *src)
{
  zOptDIRECTRect *dest;

  if( !( dest = _zOptDIRECTRectAlloc( zVecSizeNC(src->x) ) ) ) return NULL;
  zVecCopyNC( src->x, dest->x );
  dest->val = src->val;
  memcpy( dest->level, src->level, sizeof(int)*zVecSizeNC(src->x) );
  return dest;
}

/* evaluate the center value of a hyper-rectangle region */
static double _zOptDIRECTRectEval(zOptDIRECTRect *rect, double (* f)(const zVec,void*), void *util, const zVec min, const zVec max, zVec v)
{
  return rect->val = f( zVecScale( rect->x, min, max, v ), util );
}

#ifdef DEBUG
/* print out properties of a hyper-rectangle region */
static void _zOptDIRECTRectFPrint(FILE *fp, zOptDIRECTRect *rect)
{
  int i;

  fprintf( fp, "center: " );
  zVecFPrint( fp, rect->x );
  fprintf( fp, "val = %.10g\n", rect->val );
  fprintf( fp, "level:" );
  for( i=0; i<zVecSizeNC(rect->x); i++ )
    fprintf( fp, " %d", rect->level[i] );
  fprintf( fp, "\n" );
}
#endif

/* list of hyper-rectangle regions */
zListClass(zOptDIRECTRectList, zOptDIRECTRectListCell, zOptDIRECTRect*);

/* insert a hyper-rectangle region into a list in an ascending order based on the evaluation */
static bool _zOptDIRECTRectListInsert(zOptDIRECTRectList *list, zOptDIRECTRect *rect)
{
  zOptDIRECTRectListCell *cp, *ncp;

  if( !( ncp = zAlloc( zOptDIRECTRectListCell, 1 ) ) ) return false;
  ncp->data = rect;
  zListForEach( list, cp )
    if( cp->data->val > rect->val ) break;
  zListInsertPrev( list, cp, ncp );
  return true;
}

/* destroy a list of hyper-rectangle regions */
static void _zOptDIRECTRectListDestroy(zOptDIRECTRectList *list)
{
  zOptDIRECTRectListCell *cp;

  while( !zListIsEmpty( list ) ){
    zListDeleteHead( list, &cp );
    _zOptDIRECTRectFree( cp->data );
    free( cp );
  }
}

#ifdef DEBUG
/* print out a list of hyper-rectangle regions */
static void _zOptDIRECTRectListFPrint(FILE *fp, zOptDIRECTRectList *list)
{
  zOptDIRECTRectListCell *cp;
  int i = 0;

  zListForEach( list, cp ){
    fprintf( fp, "<#%d>\n", i++ );
    _zOptDIRECTRectFPrint( fp, cp->data );
  }
}

/* print out the centers of hyper-rectangle regions in a list */
static void _zOptDIRECTRectListCenterFPrint(FILE *fp, zOptDIRECTRectList *list)
{
  zOptDIRECTRectListCell *cp;
  int i;

  zListForEach( list, cp ){
    for( i=0; i<zVecSizeNC(cp->data->x); i++ )
      fprintf( fp, "%.10g ", zVecElemNC(cp->data->x,i) );
    fprintf( fp, "%.10g\n", cp->data->val );
  }
}
#endif

/* a set of hyper-rectangle regions grouped based on levels */
typedef struct{
  int level_max;
  zOptDIRECTRectList *list;
  double min_val;
  int min_level;
} zOptDIRECTRectSet;

/* allocate a set of hyper-rectangle regions */
static zOptDIRECTRectSet *_zOptDIRECTRectSetAlloc(zOptDIRECTRectSet *set, int level_max)
{
  int i;

  if( !( set->list = zAlloc( zOptDIRECTRectList, level_max+1 ) ) ) return NULL;
  for( i=0; i<=level_max; i++ )
    zListInit( &set->list[i] );
  set->level_max = level_max;
  set->min_level = 0;
  return set;
}

/* destroy a set of hyper-rectangle regions */
static void _zOptDIRECTRectSetDestroy(zOptDIRECTRectSet *set)
{
  int i;

  for( i=0; i<=set->level_max; i++ )
    _zOptDIRECTRectListDestroy( &set->list[i] );
  zFree( set->list );
  set->level_max = 0;
}

/* insert a hyper-rectangle region at a level to a set */
static bool _zOptDIRECTRectSetInsert(zOptDIRECTRectSet *set, int level, zOptDIRECTRect *rect)
{
  return _zOptDIRECTRectListInsert( &set->list[level], rect ) ? true : false;
}

/* pick up a cell of a hyper-rectangle region and delete it from a set */
static zOptDIRECTRectListCell *_zOptDIRECTRectSetDeleteCell(zOptDIRECTRectSet *set, int level)
{
  zOptDIRECTRectListCell *cell;

  zListDeleteTail( &set->list[level], &cell );
  return cell;
}

/* initialize a set of hyper-rectangle regions */
static zOptDIRECTRectSet *_zOptDIRECTRectSetInit(zOptDIRECTRectSet *set, double (* f)(const zVec,void*), void *util, int dim, const zVec min, const zVec max, zVec x)
{
  zOptDIRECTRect *rect;

  if( !( rect = _zOptDIRECTRectAlloc( dim ) ) ) return NULL;
  zVecSetAll( rect->x, 0.5 );
  _zOptDIRECTRectEval( rect, f, util, min, max, x );
  memset( rect->level, 0, sizeof(int)*dim );
  _zOptDIRECTRectSetInsert( set, 0, rect );
  return set;
}

/* find the minimum value in a set of hyper-rectangle regions */
static bool _zOptDIRECTRectSetFindMin(zOptDIRECTRectSet *set)
{
  int i = 0;

  while( zListIsEmpty( &set->list[i] ) ) i++;
  if( i > set->level_max ){
    ZRUNERROR( ZM_ERR_FATAL );
    return false;
  }
  set->min_val = zListTail( &set->list[( set->min_level = i )] )->data->val;
  for( ; i<=set->level_max; i++ ){
    if( zListIsEmpty( &set->list[i] ) ) continue;
    if( set->min_val > zListTail( &set->list[i] )->data->val )
      set->min_val = zListTail( &set->list[( set->min_level = i )] )->data->val;
  }
  return true;
}

/* a dividing pattern of a hyper-rectangle region */
typedef struct{
  bool is_active;
  zOptDIRECTRect *rect_pos;
  zOptDIRECTRect *rect_neg;
  double w;
} zOptDIRECTRectPattern;

/* allocate a dividing pattern of a hyper-rectangle region */
static zOptDIRECTRectPattern *_zOptDIRECTRectPatternAlloc(int dim)
{
  zOptDIRECTRectPattern *pat;
  int i;

  if( !( pat = zAlloc( zOptDIRECTRectPattern, dim ) ) ) return NULL;
  for( i=0; i<dim; i++ ){
    pat[i].is_active = false;
    pat[i].rect_pos = pat[i].rect_neg = NULL;
    pat[i].w = HUGE_VAL;
  }
  return pat;
}

/* destroy a dividing pattern of a hyper-rectangle region */
static void _zOptDIRECTRectPatternDestroy(zOptDIRECTRectPattern *pat, int dim)
{
  int i;

  for( i=0; i<dim; i++ ){
    pat[i].is_active = false;
    if( pat[i].rect_pos ) _zOptDIRECTRectFree( pat[i].rect_pos );
    if( pat[i].rect_neg ) _zOptDIRECTRectFree( pat[i].rect_neg );
  }
  free( pat );
}

/* create a dividing pattern of a hyper-rectangle region */
static bool _zOptDIRECTRectPatternCreate(zOptDIRECTRectPattern *pat, int i, int level, zOptDIRECTRect *rect, double (* f)(const zVec,void*), void *util, const zVec min, const zVec max, zVec x)
{
  double d;

  pat->is_active = true;
  pat->rect_pos = _zOptDIRECTRectClone( rect );
  pat->rect_neg = _zOptDIRECTRectClone( rect );
  if( !pat->rect_pos || !pat->rect_neg ){
    pat->rect_pos = pat->rect_neg = NULL;
    return false;
  }
  d = pow( 3.0, -(level+1) );
  zVecElemNC(pat->rect_pos->x,i) += d;
  _zOptDIRECTRectEval( pat->rect_pos, f, util, min, max, x );
  zVecElemNC(pat->rect_neg->x,i) -= d;
  _zOptDIRECTRectEval( pat->rect_neg, f, util, min, max, x );
  pat->w = zMin( pat->rect_pos->val, pat->rect_neg->val );
  return true;
}

/* divide a hyper-rectangle region */
static bool _zOptDIRECTRectSetRectDivide(zOptDIRECTRectSet *set, int level, zOptDIRECTRect *rect, double (* f)(const zVec,void*), void *util, const zVec min, const zVec max, zVec x)
{
  zOptDIRECTRectPattern *pat;
  zIntList index;
  zIntListCell *cp, *ncp;
  bool ret = true;
  int i;

  zListInit( &index );
  if( !( pat = _zOptDIRECTRectPatternAlloc( zVecSizeNC(rect->x) ) ) ) return false;
  for( i=0; i<zVecSizeNC(rect->x); i++ ){
    if( rect->level[i] != level ) continue;
    /* find directions of longest edges */
    if( !_zOptDIRECTRectPatternCreate( &pat[i], i, level, rect, f, util, min, max, x ) ){
      ret = false;
      goto TERMINATE;
    }
    zListForEach( &index, cp )
      if( pat[cp->data].w > pat[i].w ) break;
    if( !( ncp = zAlloc( zIntListCell, 1 ) ) ){
      ret = false;
      goto TERMINATE;
    }
    ncp->data = i;
    zListInsertPrev( &index, cp, ncp );
  }
  /* divide the region in directions of longest edges */
  zListForEach( &index, cp ){
    rect->level[cp->data] = level + 1;
    memcpy( pat[cp->data].rect_pos->level, rect->level, sizeof(int)*zVecSizeNC(rect->x) );
    memcpy( pat[cp->data].rect_neg->level, rect->level, sizeof(int)*zVecSizeNC(rect->x) );
    if( cp == zListHead(&index) ) level++;
    if( !_zOptDIRECTRectSetInsert( set, level, pat[cp->data].rect_pos ) ||
        !_zOptDIRECTRectSetInsert( set, level, pat[cp->data].rect_neg ) ){
      ret = false;
      goto TERMINATE;
    }
    pat[cp->data].rect_pos = pat[cp->data].rect_neg = NULL;
  }
  if( !_zOptDIRECTRectSetInsert( set, level, rect ) ) ret = false;
 TERMINATE:
  zListDestroy( zIntListCell, &index );
  _zOptDIRECTRectPatternDestroy( pat, zVecSizeNC(rect->x) );
  return ret;
}

#ifdef DEBUG
/* print out a set of hyper-rectangle regions */
static void _zOptDIRECTRectSetFPrint(FILE *fp, zOptDIRECTRectSet *set)
{
  int i = 0;

  fprintf( fp, "*** rectangle set ***\n" );
  for( i=0; i<=set->level_max; i++ ){
    fprintf( fp, "[level=%d]\n", i );
    _zOptDIRECTRectListFPrint( fp, &set->list[i] );
  }
}

/* print out centers of hyper-rectangle regions in a set */
static void _zOptDIRECTRectSetCenterFPrint(FILE *fp, zOptDIRECTRectSet *set)
{
  int i;

  for( i=0; i<=set->level_max; i++ )
    _zOptDIRECTRectListCenterFPrint( fp, &set->list[i] );
}
#endif

/* a set of potentially optimal hyper-rectangle regions */
typedef struct{
  int level_max;
  zOptDIRECTRectListCell **cell;
} zOptDIRECTRectPOR;

/* allocate memory for a set of potentially optimal hyper-rectangle regions */
static zOptDIRECTRectPOR *_zOptDIRECTRectPORAlloc(zOptDIRECTRectPOR *por, int level_max)
{
  int i;

  if( !( por->cell = zAlloc( zOptDIRECTRectListCell*, level_max+1 ) ) ) return NULL;
  por->level_max = level_max;
  for( i=0; i<=level_max; i++ )
    por->cell[i] = NULL;
  return por;
}

/* free memory of a set of potentially optimal hyper-rectangle regions */
static void _zOptDIRECTRectPORFree(zOptDIRECTRectPOR *por)
{
  zFree( por->cell );
  por->level_max = 0;
}

/* set a cell of a set of potentially optimal hyper-rectangle regions */
static void _zOptDIRECTRectPORSetCell(zOptDIRECTRectPOR *por, int level, zOptDIRECTRectListCell *cell)
{
  por->cell[level] = cell;
}

/* identify potentially optimal hyper-rectangle regions */
static bool _zOptDIRECTRectPORIdent(zOptDIRECTRectPOR *por, zOptDIRECTRectSet *set)
{
  zOptDIRECTRectListCell *cell;
  double gmax, g;
  int i = 0, j, k;

  if( !_zOptDIRECTRectSetFindMin( set ) ) return false;
  if( set->min_level == set->level_max ) return true;
  while( zListIsEmpty( &set->list[i] ) ) i++;
  while( 1 ){
    cell = _zOptDIRECTRectSetDeleteCell( set, i );
    _zOptDIRECTRectPORSetCell( por, i, cell );
    if( i == set->min_level ) break;
    gmax = 0;
    j = i + 1;
    for( k=j; k<=set->min_level; k++ ){
      if( zListIsEmpty( &set->list[k] ) ) continue;
      g = pow(3.0,i) * ( cell->data->val - zListTail(&set->list[k])->data->val ) / ( 1 - pow(3.0,i-k) );
      if( g >= gmax ){
        gmax = g;
        j = k;
      }
    }
    i = j;
  }
  return true;
}

/* divide potentially optimal hyper-rectangle regions */
static bool _zOptDIRECTRectPORSetDivide(zOptDIRECTRectPOR *por, zOptDIRECTRectSet *set, double (* f)(const zVec,void*), void *util, const zVec min, const zVec max, zVec x)
{
  int i;

  for( i=0; i<por->level_max; i++ ){
    if( por->cell[i] )
      if( !_zOptDIRECTRectSetRectDivide( set, i, por->cell[i]->data, f, util, min, max, x ) ) return false;
    por->cell[i] = NULL;
  }
  return true;
}

#ifdef DEBUG
/* print out a set of potentially optimal hyper-rectangle regions */
static void _zOptDIRECTRectPORFPrint(FILE *fp, zOptDIRECTRectPOR *por)
{
  int i;

  fprintf( fp, "*** potentially optimal rectangle set ***\n" );
  for( i=0; i<=por->level_max; i++ ){
    fprintf( fp, "[level=%d]\n", i );
    if( por->cell[i] )
      _zOptDIRECTRectFPrint( fp, por->cell[i]->data );
  }
}
#endif

/* solve an optimization problem by DIRECT method. */
int zOptSolveDIRECT(double (* f)(const zVec,void*), void *util, const zVec min, const zVec max, int iter, double tol, zVec ans, double *eval)
{
  zOptDIRECTRectSet set;
  zOptDIRECTRectPOR por;
  int level_max;
  int i = 0;

  if( !zVecSizeEqual( min, ans ) || !zVecSizeEqual( max, ans ) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return -1;
  }
  level_max = floor( -zLog( 3, tol ) );
  if( !_zOptDIRECTRectSetAlloc( &set, level_max ) ) return -1;
  if( !_zOptDIRECTRectPORAlloc( &por, level_max ) ) goto TERMINATE2;
  if( !_zOptDIRECTRectSetInit( &set, f, util, zVecSizeNC(ans), min, max, ans ) )
    goto TERMINATE2;

  if( iter == 0 ) iter = ZOPT_DIRECT_MAX_ITER_NUM;
  for( i=0; i<iter; i++ ){
    if( !_zOptDIRECTRectPORIdent( &por, &set ) ) goto TERMINATE1;
    if( set.min_level == set.level_max ) break;
    _zOptDIRECTRectPORSetDivide( &por, &set, f, util, min, max, ans );
  }
  zVecScale( zListTail( &set.list[set.min_level] )->data->x, min, max, ans );
  *eval = zListTail( &set.list[set.min_level] )->data->val;
 TERMINATE1:
  _zOptDIRECTRectPORFree( &por );
 TERMINATE2:
  _zOptDIRECTRectSetDestroy( &set );
  return i;
}
