/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_pex - polynomial expression class.
 */

#include <zm/zm_pex.h>

/* ********************************************************** */
/* CLASS: zPex
 * polynomial expression class
 * ********************************************************** */

/* zPexRgl
 * - regulate polynomial.
 */
zPex zPexRgl(zPex *p)
{
  int size;
  zPex newp;

  for( size=zPexDim(*p); zPexCoeff(*p,size)==0; size-- );
  if( size != zPexDim(*p) ){
    zPexSetDim( *p, size );
    if( !( newp = zPexAlloc( size ) ) ){
      ZALLOCERROR();
      return NULL;
    }
    zPexCopy( *p, newp );
    zPexFree( *p );
    *p = newp;
  }
  return *p;
}

/* zPexAddDRC
 * - direct addision of polynomial expression.
 */
zPex zPexAddDRC(zPex p1, zPex p2)
{
  register int i, dim;

  if( zPexDim(p1) < ( dim = zPexDim(p2) ) ){
    ZRUNERROR( ZM_ERR_PEX_DIMMIS );
    return NULL;
  }
  for( i=0; i<=dim; i++ )
    zPexCoeff(p1,i) += zPexCoeff(p2,i);
  return p1;
}

/* zPexSubDRC
 * - direct substraction of polynomial expression.
 */
zPex zPexSubDRC(zPex p1, zPex p2)
{
  register int i, dim;

  if( zPexDim(p1) < ( dim = zPexDim(p2) ) ){
    ZRUNERROR( ZM_ERR_PEX_DIMMIS );
    return NULL;
  }
  for( i=0; i<=dim; i++ )
    zPexCoeff(p1,i) -= zPexCoeff(p2,i);
  return p1;
}

/* zPexAdd
 * - addition of polynomial expression.
 */
zPex zPexAdd(zPex p1, zPex p2)
{
  zPex p;

  if( !( p = zPexAlloc( zMax( zPexDim(p1), zPexDim(p2) ) ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  zPexAddDRC( p, p1 );
  zPexAddDRC( p, p2 );
  zPexRgl( &p );
  return p;
}

/* zPexSub
 * - substraction of polynomial expression.
 */
zPex zPexSub(zPex p1, zPex p2)
{
  zPex p;

  if( !( p = zPexAlloc( zMax( zPexDim(p1), zPexDim(p2) ) ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  zPexAddDRC( p, p1 );
  zPexSubDRC( p, p2 );
  zPexRgl( &p );
  return p;
}

/* zPexMul
 * - multiplication of polynomial expression.
 */
zPex zPexMul(zPex p1, zPex p2)
{
  register int i, j, dim1, dim2;
  zPex p;

  dim1 = zPexDim( p1 );
  dim2 = zPexDim( p2 );
  if( !( p = zPexAlloc( dim1 + dim2 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  for( i=0; i<=dim2; i++ )
    for( j=0; j<=dim1; j++ )
      zPexCoeff(p,i+j) += zPexCoeff(p2,i)*zPexCoeff(p1,j);
  zPexRgl( &p );
  return p;
}

static void _zPexDivDRC(zPex p, zPex f, zPex q, zPex *r);

/* (static)
 * _zPexDivDRC
 * - division of polynomial expression.
 */
void _zPexDivDRC(zPex p, zPex f, zPex q, zPex *r)
{
  register int i, n, m;
  double a;

  n = zPexDim( p );
  m = zPexDim( f );
  if( n < m ){
    zPexCopy( p, *r );
    zPexRgl( r );
    return;
  }
  a = zPexCoeff( p, n ) / zPexCoeff( f, m );
  zPexSetCoeff( q, n-m, a );
  zPexSetDim( q, n-m-1 );
  for( i=1; i<=m; i++ )
    zPexSetCoeff( p, n-i, zPexCoeff(p,n-i) - a*zPexCoeff(f,m-i) );
  zPexSetDim( p, n-1 );
  _zPexDivDRC( p, f, q, r );
  zPexSetDim( p, n );
  zPexSetDim( q, n-m );
  return;
}

/* zPexDiv
 * - division of polynomial expression.
 */
bool zPexDiv(zPex p, zPex f, zPex *q, zPex *r)
{
  zPex pcp;
  int dim;
  bool ret = true;

  dim = zPexDim(p) - zPexDim(f);
  if( !( pcp = zPexAlloc( zPexDim(p) ) ) ){
    ZALLOCERROR();
    return false;
  }
  zPexCopy( p, pcp );
  *q = zPexAlloc( dim );
  *r = zPexAlloc( zPexDim(f) - 1 );
  if( *q && *r )
    _zPexDivDRC( pcp, f, *q, r );
  else{
    ZALLOCERROR();
    ret = false;
  }
  zPexFree( pcp );
  return ret;
}

static zPex _zPexExp(zVec factor);
/* (static)
 * _zPexExp
 * - expand factors into polynomial expression.
 *   (internal function which may modify 'factor')
 */
zPex _zPexExp(zVec factor)
{
  zPex p1, p2, p = NULL;
  uint size, hsize;

  size = zVecSize( factor );
  if( size == 1 ){
    if( !( p = zPexAlloc( 1 ) ) ){
      ZALLOCERROR();
      return NULL;
    }
    zPexSetCoeff( p, 0,-zVecElem( factor, 0 ) );
    zPexSetCoeff( p, 1, 1 );
    return p;
  }
  hsize = size / 2;
  zVecSetSize( factor, hsize );
  p1 = _zPexExp( factor );
  zVecSetSize( factor, size-hsize );
  memcpy( zVecBuf(factor), zVecBuf(factor)+hsize, sizeof(double)*zVecSizeNC(factor) );
  p2 = _zPexExp( factor );
  if( !p1 || !p2 ) goto TERMINATE;
  p = zPexMul( p1, p2 );

 TERMINATE:
  zPexFree( p1 );
  zPexFree( p2 );
  return p;
}

/* zPexExp
 * - expand factors into polynomial expression.
 */
zPex zPexExp(zVec factor)
{
  zVec f;
  zPex result;

  if( !( f = zVecClone( factor ) ) ) return NULL;
  result = _zPexExp( f );
  zVecFree( f );
  return result;
}

/* zPexModulo
 * - modulo of a primary expression.
 */
zPex zPexModulo(zPex p1, double a, zPex p2)
{
  register int i, j;

  if( !zPexCopy( p1, p2 ) ){
    ZRUNERROR( ZM_ERR_PEX_DIMMIS );
    return NULL;
  }
  for( i=zPexDim(p2)-1; i>=0; i-- )
    for( j=i; j<zPexDim(p2); j++ )
      zPexCoeff(p2,j) -= a * zPexCoeff(p2,j+1);
  return p2;
}

/* zPexDif
 * - differentiation of polynomial expression.
 */
zPex zPexDif(zPex p)
{
  register int i, dim;
  zPex q;

  if( !( q = zPexAlloc( ( dim = zPexDim(p) - 1 ) ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  for( i=0; i<=dim; i++ )
    zPexSetCoeff( q, i, zPexCoeff(p,i+1)*(i+1) );
  return q;
}

/* zPexIntg
 * - integration of polynomial expression.
 */
zPex zPexIntg(zPex p)
{
  register int i, dim;
  zPex q;

  if( !( q = zPexAlloc( ( dim=zPexDim(p) ) + 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  for( i=0; i<=dim; i++ )
    zPexSetCoeff( q, i+1, zPexCoeff(p,i)/(i+1) );
  return q;
}

/* zPexVal
 * - evaluate polynomial expression for one argument.
 */
double zPexVal(zPex p, double arg)
{
  register int i;
  double result;

  i = zPexDim( p );
  result = zPexCoeff( p, i );
  for( ; i>0; i-- )
    result = zPexCoeff( p, i-1 ) + result * arg;
  return result;
}

/* zPexCVal
 * - evaluate polynomial expression for a complex number argument.
 */
zComplex *zPexCVal(zPex p, zComplex *arg, zComplex *c)
{
  register int i;
  zComplex tmp;

  i = zPexDim(p);
  zComplexCreate( c, zPexCoeff(p,i), 0 );
  for( ; i>0; i-- ){
    zComplexCMul( c, arg, &tmp );
    zComplexCreate( c, zPexCoeff(p,i-1)+tmp.re, tmp.im );
  }
  return c;
}

/* zPexDifVal
 * - evaluate differential value of polynomial expression.
 */
double zPexDifVal(zPex p, int dim, double arg)
{
  register int i;
  double c, result;

  i = zPexDim( p );
  c = zPermut( i, dim );
  result = c * zPexCoeff( p, i );
  for( ; i>dim; i-- ){
    c *= (double)( i - dim ) / i;
    result = c * zPexCoeff( p, i-1 ) + result * arg;
  }
  return result;
}

/* zPexFRead
 * - read a polynomial expression from file.
 */
zPex zPexFRead(FILE *fp)
{
  register int i, dim;
  zPex p;

  dim = zFInt( fp );
  if( !( p = zPexAlloc( dim ) ) ) return NULL;

  for( i=0; i<=dim; i++ )
    zPexSetCoeff( p, i, zFDouble( fp ) );
  return p;
}

/* zPexFWrite
 * - write a polynomial expression to file.
 */
void zPexFWrite(FILE *fp, zPex p)
{
  register int i;

  if( !p ) return;
  fprintf( fp, "%d (", zPexDim(p) );
  for( i=0; i<=zPexDim(p); i++ )
    fprintf( fp, " %.10g", zPexCoeff(p,i) );
  fprintf( fp, " )\n" );
}

/* zPexFExpr
 * - expression of polynomial expression.
 */
void zPexFExpr(FILE *fp, zPex p, char c)
{
  register int i, dim;

  if( ( dim = zPexDim(p) ) < 0 )
    fprintf( fp, "0" );
  else
    for( i=dim; i>=0; i-- ){
      if( i<dim && zPexCoeff(p,i)>0 )
        fprintf( fp, " + " );
      if( zPexCoeff(p,i) != 0 ){
        if( zPexCoeff(p,i) != 1 || i == 0 )
          fprintf( fp, "%.10g ", zPexCoeff(p,i) );
        if( i > 0 ) fprintf( fp, "%c", c );
        if( i > 1 ) fprintf( fp, "^%d", i );
      }
    }
  fprintf( fp, "\n" );
}
