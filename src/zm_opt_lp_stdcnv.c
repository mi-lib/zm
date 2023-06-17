/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_lp_stdcnv - optimization tools: problem convertor
 * to the standard form of linear programming.
 */

#include <zm/zm_opt.h>

/* convert linear programming problem from inequality constraint form
 * to standard equation form.
 */
bool zLPIneq2Std(zMat a, zVec c, zVec x, zMat *as, zVec *cs, zVec *xs)
{
  int i;

  *as = zMatAlloc( zMatRowSizeNC(a), zMatColSizeNC(a)+zMatRowSizeNC(a) );
  *cs = zVecAlloc( zVecSizeNC(c)+zMatRowSizeNC(a) );
  *xs = zVecAlloc( zVecSizeNC(x)+zMatRowSizeNC(a) );
  if( !*as || !*cs || !*xs ){
    zMatFree( *as );
    zVecFree( *cs );
    zVecFree( *xs );
    return false;
  }
  zMatPut( *as, 0, 0, a );
  for( i=0; i<zMatRowSizeNC(a); i++ )
    zMatSetElemNC( *as, i, zMatColSizeNC(a)+i, 1.0 );
  zVecPut( *cs, 0, c );
  return true;
}

/* convert linear programming problem from inequality constraint form
 * with unbounded variables to standard equation form with positive variables.
 */
bool zLPUnb2Std(zMat a, zVec c, zVec x, zMat *as, zVec *cs, zVec *xs)
{
  int i, j;

  *as = zMatAlloc( zMatRowSizeNC(a), zMatColSizeNC(a)*2+zMatRowSizeNC(a) );
  *cs = zVecAlloc( zVecSizeNC(c)*2+zMatRowSizeNC(a) );
  *xs = zVecAlloc( zVecSizeNC(x)*2+zMatRowSizeNC(a) );
  if( !*as || !*cs || !*xs ){
    zMatFree( *as );
    zVecFree( *cs );
    zVecFree( *xs );
    return false;
  }
  for( i=0; i<zMatRowSizeNC(a); i++ ){
    for( j=0; j<zMatColSizeNC(a); j++ ){
      zMatSetElemNC( *as, i, j, zMatElem(a,i,j) );
      zMatSetElemNC( *as, i, zMatColSizeNC(a)+j, -zMatElem(a,i,j) );
    }
    zMatSetElemNC( *as, i, zMatColSizeNC(a)*2+i, 1.0 );
  }
  for( i=0; i<zVecSizeNC(c); i++ ){
    zVecSetElemNC( *cs, i, zVecElem(c,i) );
    zVecSetElemNC( *cs, zVecSizeNC(c)+i, -zVecElem(c,i) );
  }
  return true;
}
