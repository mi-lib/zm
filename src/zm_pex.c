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

/* regulate polynomial. */
zPex zPexRgl(zPex *p)
{
  int size;
  zPex newp;

  for( size=zPexDim(*p); zPexCoeff(*p,size)==0; size-- );
  if( size != zPexDim(*p) ){
    zPexSetDim( *p, size );
    if( !( newp = zPexAlloc( size ) ) ) return NULL;
    zPexCopy( *p, newp );
    zPexFree( *p );
    *p = newp;
  }
  return *p;
}

/* add a polynomial expression directly to another. */
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

/* subtract a polynomial expression directly from another. */
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

/* add two polynomial expressions. */
zPex zPexAdd(zPex p1, zPex p2)
{
  zPex p;

  if( !( p = zPexAlloc( zMax( zPexDim(p1), zPexDim(p2) ) ) ) ) return NULL;
  zPexAddDRC( p, p1 );
  zPexAddDRC( p, p2 );
  zPexRgl( &p );
  return p;
}

/* substract a polynomial expression from another. */
zPex zPexSub(zPex p1, zPex p2)
{
  zPex p;

  if( !( p = zPexAlloc( zMax( zPexDim(p1), zPexDim(p2) ) ) ) ) return NULL;
  zPexAddDRC( p, p1 );
  zPexSubDRC( p, p2 );
  zPexRgl( &p );
  return p;
}

/* multiply a polynomial expression by another. */
zPex zPexMul(zPex p1, zPex p2)
{
  register int i, j, dim1, dim2;
  zPex p;

  dim1 = zPexDim( p1 );
  dim2 = zPexDim( p2 );
  if( !( p = zPexAlloc( dim1 + dim2 ) ) ) return NULL;
  for( i=0; i<=dim2; i++ )
    for( j=0; j<=dim1; j++ )
      zPexCoeff(p,i+j) += zPexCoeff(p2,i)*zPexCoeff(p1,j);
  zPexRgl( &p );
  return p;
}

/* divide a polynomial expression directly by another. */
static void _zPexDivDRC(zPex p, zPex f, zPex q, zPex *r)
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

/* divide a polynomial expression by another. */
bool zPexDiv(zPex p, zPex f, zPex *q, zPex *r)
{
  zPex pcp;
  int dim;
  bool ret = true;

  dim = zPexDim(p) - zPexDim(f);
  if( !( pcp = zPexAlloc( zPexDim(p) ) ) ) return false;
  zPexCopy( p, pcp );
  *q = zPexAlloc( dim );
  *r = zPexAlloc( zPexDim(f) - 1 );
  if( *q && *r )
    _zPexDivDRC( pcp, f, *q, r );
  else
    ret = false;
  zPexFree( pcp );
  return ret;
}

/* expand factors into polynomial expression.
 * (internal function to modify \a factor) */
static zPex _zPexExp(zVec factor)
{
  zPex p1, p2, p = NULL;
  uint size, hsize;

  size = zVecSizeNC( factor );
  if( size == 1 ){
    if( !( p = zPexAlloc( 1 ) ) ) return NULL;
    zPexSetCoeff( p, 0,-zVecElemNC(factor,0) );
    zPexSetCoeff( p, 1, 1 );
    return p;
  }
  hsize = size / 2;
  zVecSetSizeNC( factor, hsize );
  p1 = _zPexExp( factor );
  zVecSetSizeNC( factor, size-hsize );
  memcpy( zVecBufNC(factor), zVecBufNC(factor)+hsize, sizeof(double)*zVecSizeNC(factor) );
  p2 = _zPexExp( factor );
  if( !p1 || !p2 ) goto TERMINATE;
  p = zPexMul( p1, p2 );

 TERMINATE:
  zPexFree( p1 );
  zPexFree( p2 );
  return p;
}

/* expand factors into a polynomial expression. */
zPex zPexExp(zVec factor)
{
  zVec f;
  zPex result;

  if( !( f = zVecClone( factor ) ) ) return NULL;
  result = _zPexExp( f );
  zVecFree( f );
  return result;
}

/* expand imaginary factors into polynomial expression.
 * conjugate complex numbers are supposed to be paired as adjacencies.
 */
static zPex _zPexExpIm(zCVec ifactor)
{
  zPex p1, p2, p = NULL;
  uint size, hsize, hhsize;

  hsize = ( size = zVecSizeNC( ifactor ) ) / 2;
  if( hsize == 1 ){
    if( !( p = zPexAlloc( 2 ) ) ) return NULL;
    zPexSetCoeff( p, 0, zComplexSqrAbs(zCVecElemNC(ifactor,0)) );
    zPexSetCoeff( p, 1,-2*zCVecElemNC(ifactor,0)->re );
    zPexSetCoeff( p, 2, 1 );
    return p;
  }
  hhsize = hsize / 2;
  zVecSetSizeNC( ifactor, ( hsize = hhsize * 2 ) );
  p1 = _zPexExpIm( ifactor );
  zVecSetSizeNC( ifactor, size-hsize );
  memcpy( zCVecBufNC(ifactor), zCVecBufNC(ifactor)+hsize, sizeof(zComplex)*zCVecSizeNC(ifactor) );
  p2 = _zPexExpIm( ifactor );
  if( !p1 || !p2 ) goto TERMINATE;
  p = zPexMul( p1, p2 );

 TERMINATE:
  zPexFree( p1 );
  zPexFree( p2 );
  return p;
}

/* expand complex factors into a polynomial expression. */
zPex zPexCExp(zCVec factor)
{
  zVec rfactor;
  zCVec ifactor;
  zPex pr, pc, p = NULL;

  if( !zCVecToReIm( factor, &rfactor, &ifactor ) ) return NULL;
  pr = rfactor ? _zPexExp( rfactor ) : NULL;
  pc = ifactor ? _zPexExpIm( ifactor ) : NULL;
  if( pr && pc )
    p = zPexMul( pr, pc );
  else if( pr )
    p = zPexClone( pr );
  else
    p = zPexClone( pc );
  zPexFree( pr );
  zPexFree( pc );
  zVecFree( rfactor );
  zCVecFree( ifactor );
  return p;
}

/* modulo of a primary expression. */
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

/* differentiate a polynomial expression. */
zPex zPexDif(zPex p)
{
  register int i, dim;
  zPex q;

  if( !( q = zPexAlloc( ( dim = zPexDim(p) - 1 ) ) ) ) return NULL;
  for( i=0; i<=dim; i++ )
    zPexSetCoeff( q, i, zPexCoeff(p,i+1)*(i+1) );
  return q;
}

/* integrate a polynomial expression. */
zPex zPexIntg(zPex p)
{
  register int i, dim;
  zPex q;

  if( !( q = zPexAlloc( ( dim=zPexDim(p) ) + 1 ) ) ) return NULL;
  for( i=0; i<=dim; i++ )
    zPexSetCoeff( q, i+1, zPexCoeff(p,i)/(i+1) );
  return q;
}

/* evaluate a polynomial expression given an argument. */
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

/* evaluate a polynomial expression given a complex number argument. */
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

/* evaluate the differential value of a polynomial expression. */
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

/* scan a polynomial expression from a ZTK processor. */
zPex zPexFromZTK(ZTK *ztk)
{
  register int i, dim;
  zPex p;

  if( !( p = zPexAlloc( ( dim = ZTKInt(ztk) ) ) ) ) return NULL;
  for( i=0; i<=dim; i++ )
    zPexSetCoeff( p, i, ZTKDouble(ztk) );
  return p;
}

/* scan a polynomial expression from a file. */
zPex zPexFScan(FILE *fp)
{
  register int i, dim;
  zPex p;

  if( !( p = zPexAlloc( ( dim = zFInt(fp) ) ) ) ) return NULL;
  for( i=0; i<=dim; i++ )
    zPexSetCoeff( p, i, zFDouble(fp) );
  return p;
}

/* print a polynomial expression to a file. */
void zPexFPrint(FILE *fp, zPex p)
{
  register int i;

  if( !p ) return;
  fprintf( fp, "%d (", zPexDim(p) );
  for( i=0; i<=zPexDim(p); i++ )
    fprintf( fp, " %.10g", zPexCoeff(p,i) );
  fprintf( fp, " )\n" );
}

/* present a polynomial expression. */
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
