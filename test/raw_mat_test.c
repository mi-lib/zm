#include <zm/zm_raw.h>

void set_raw_mat_all(double *m, int colcapacity, int rowsize, int colsize, double val)
{
  int i, j;

  for( i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      m[i*colcapacity+j] = val;
}

void assert_raw_mat_zero(void)
{
  double *m;
  const int rowcapacity = 10;
  const int colcapacity = 15;
  const int rowsize = 7;
  const int colsize =10;
  int i, j;
  bool result1, result2;

  m = zAlloc( double, rowcapacity * colcapacity );
  zRawMatZero( m, colcapacity, rowcapacity, colcapacity );
  set_raw_mat_all( m, colcapacity, rowsize, colsize, 1.0 );
  for( result1=true, i=0; i<rowcapacity; i++ )
    for( j=0; j<colcapacity; j++ )
      if( ( i < rowsize && j < colsize && m[i*colcapacity+j] != 1 ) ||
          ( ( i >= rowsize || j >= colsize ) && m[i*colcapacity+j] != 0 ) ) result1 = false;
  set_raw_mat_all( m, colcapacity, rowcapacity, colcapacity, 1.0 );
  zRawMatZero( m, colcapacity, rowsize, colsize );
  for( result2=true, i=0; i<rowcapacity; i++ )
    for( j=0; j<colcapacity; j++ )
      if( ( i < rowsize && j < colsize && m[i*colcapacity+j] != 0 ) ||
          ( ( i >= rowsize || j >= colsize ) && m[i*colcapacity+j] != 1 ) ) result2 = false;
  free( m );
  zAssert( zRawMatZero, result1 && result2 );
}

void assert_raw_mat_touchup(void)
{
  double *m;
  const int rowcapacity = 10;
  const int colcapacity = 15;
  const int rowsize = 7;
  const int colsize =10;
  int i, j;
  bool result1, result2;

  m = zAlloc( double, rowcapacity * colcapacity );
  set_raw_mat_all( m, colcapacity, rowcapacity, colcapacity, 1.0 );
  set_raw_mat_all( m, colcapacity, rowsize, colsize, 0.01 );
  for( result1=true, i=0; i<rowcapacity; i++ )
    for( j=0; j<colcapacity; j++ )
      if( ( i < rowsize && j < colsize && m[i*colcapacity+j] != 0.01 ) ||
          ( ( i >= rowsize || j >= colsize ) && m[i*colcapacity+j] != 1 ) ) result1 = false;
  zRawMatTouchup( m, colcapacity, rowsize, colsize, 0.1 );
  for( result2=true, i=0; i<rowcapacity; i++ )
    for( j=0; j<colcapacity; j++ )
      if( ( i < rowsize && j < colsize && m[i*colcapacity+j] != 0 ) ||
          ( ( i >= rowsize || j >= colsize ) && m[i*colcapacity+j] != 1 ) ) result2 = false;
  free( m );
  zAssert( zRawMatTouchup, result1 && result2 );
}

void assert_raw_mat_rand(void)
{
  double *m, *mmin, *mmax;
  const int rowcapacity = 8;
  const int colcapacity =10;
  const int rowsize = 5;
  const int colsize = 7;
  const int n = 100;
  int i, j, k;
  bool result1, result2;

  m = zAlloc( double, rowcapacity * colcapacity );
  mmin = zAlloc( double, rowcapacity * colcapacity );
  mmax = zAlloc( double, rowcapacity * colcapacity );
  zRawMatZero( m, colcapacity, rowcapacity, colcapacity );
  for( i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ ){
      mmin[i*colcapacity*j] = i*colcapacity + j;
      mmax[i*colcapacity*j] = mmin[i*colcapacity*j] + 5;
    }
  for( result1=true, k=0; k<n; k++ ){
    zRawMatRand( m, colcapacity, mmin, colcapacity, mmax, colcapacity, rowsize, colsize );
    for( i=0; i<rowcapacity; i++ )
      for( j=0; j<colcapacity; j++ )
        if( i < rowsize && j < colsize &&
            ( m[i*colcapacity+j] < mmin[i*colcapacity+j] ||
              m[i*colcapacity+j] > mmax[i*colcapacity+j] ) ) result1 = false;
  }
  for( result2=true, k=0; k<n; k++ ){
    zRawMatRandUniform( m, colcapacity, -10, 10, rowsize, colsize );
    for( i=0; i<rowcapacity; i++ )
      for( j=0; j<colcapacity; j++ )
        if( i < rowsize && j < colsize &&
            ( m[i*colcapacity+j] <-10 ||
              m[i*colcapacity+j] > 10 ) ) result2 = false;
  }
  free( m );
  free( mmin );
  free( mmax );
  zAssert( zRawMatRand, result1 );
  zAssert( zRawMatRandUniform, result2 );
}

void assert_raw_mat_copy(void)
{
  double *m1, *m2;
  const int rowcapacity1 = 10;
  const int colcapacity1 = 15;
  const int rowcapacity2 =  8;
  const int colcapacity2 =  9;
  int rowsize = 5;
  int colsize = 7;
  bool result;

  m1 = zAlloc( double, rowcapacity1 * colcapacity1 );
  m2 = zAlloc( double, rowcapacity2 * colcapacity2 );
  zRawMatRandUniform( m1, colcapacity1, -10, 10, rowsize, colsize );
  zRawMatCopy( m1, colcapacity1, m2, colcapacity2, rowsize, colsize );
  result = zRawMatMatch( m1, colcapacity1, m2, colcapacity2, rowsize, colsize );
  free( m1 );
  free( m2 );
  zAssert( zRawMatCopy, result );
}

void assert_raw_mat_equal(void)
{
  double *m1, *m2;
  const int rowcapacity1 = 10;
  const int colcapacity1 = 15;
  const int rowcapacity2 =  8;
  const int colcapacity2 =  9;
  int rowsize = 5;
  int colsize = 7;
  int i, j;
  bool result;

  m1 = zAlloc( double, rowcapacity1 * colcapacity1 );
  m2 = zAlloc( double, rowcapacity2 * colcapacity2 );
  zRawMatRandUniform( m1, colcapacity1, -10, 10, rowsize, colsize );
  zRawMatCopy( m1, colcapacity1, m2, colcapacity2, rowsize, colsize );
  for( i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      m2[i*colcapacity2+j] += 0.1*zTOL;
  result = zRawMatEqual( m1, colcapacity1, m2, colcapacity2, rowsize, colsize, zTOL );
  for( i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      m2[i*colcapacity2+j] += 10*zTOL;
  if( zRawMatEqual( m1, colcapacity1, m2, colcapacity2, rowsize, colsize, zTOL ) ) result = false;

  free( m1 );
  free( m2 );
  zAssert( zRawMatEqual, result );
}

void assert_raw_mat_is_tol(void)
{
  double *m;
  const int rowcapacity = 4;
  const int colcapacity = 5;
  const int rowsize = 3;
  const int colsize = 3;
  int i, j;
  bool result;

  m = zAlloc( double, rowcapacity * colcapacity );
  for( i=0; i<rowcapacity; i++ )
    for( j=0; j<colcapacity; j++ )
      m[i*colcapacity+j] = 1.0;
  for( i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      m[i*colcapacity+j] = zTOL*0.1;
  result = zRawMatIsTiny( m, colcapacity, rowsize, colsize );
  for( i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      m[i*colcapacity+j] += zTOL;
  if( zRawMatIsTiny( m, colcapacity, rowsize, colsize ) ) result = false;
  free( m );
  zAssert( zRawMatIsTiny, result );
}

void assert_raw_mat_diag(void)
{
  double *m;
  double diag[] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
  const int rowcapacity = 8;
  const int colcapacity = 10;
  const int rowsize = 5;
  const int colsize = 7;
  int i, j;
  bool result1, result2;

  m = zAlloc( double, rowcapacity * colcapacity );
  zRawMatIdent( m, colcapacity, rowsize, colsize );
  for( result1=true, i=0; i<rowcapacity; i++ )
    for( j=0; j<colcapacity; j++ )
      if( i < rowsize && j < colsize && i == j ){
        if( !zEqual( m[i*colcapacity+j], 1.0, zTOL ) ) result1 = false;
      } else
        if( m[i*colcapacity+j] != 0 ) result1 = false;
  zRawMatDiag( m, colcapacity, diag, rowsize, colsize );
  for( result2=true, i=0; i<rowcapacity; i++ )
    for( j=0; j<colcapacity; j++ )
      if( i < rowsize && j < colsize && i == j ){
        if( !zEqual( m[i*colcapacity+j], i + 1.0, zTOL ) ) result2 = false;
      } else
        if( m[i*colcapacity+j] != 0 ) result2 = false;
  free( m );
  zAssert( zRawMatIdent, result1 );
  zAssert( zRawMatDiag, result2 );
}

void assert_raw_mat_get_put_row_col(void)
{
  double *m, *v, *v_ans;
  const int rowcapacity = 6;
  const int colcapacity = 8;
  const int rowsize = 4;
  const int colsize = 5;
  int i, j;
  bool result1, result2, result3, result4;

  m = zAlloc( double, rowcapacity * colcapacity );
  v = zAlloc( double, zMax( rowcapacity, colcapacity ) );
  v_ans = zAlloc( double, zMax( rowcapacity, colcapacity ) );
  for( i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      m[i*colcapacity+j] = i + j;
  result1 = true;
  for( i=0; i<colsize; i++ ) v_ans[i] = i;
  zRawMatGetRow( m, colcapacity, rowsize, colsize, 0, v );
  if( !zRawVecMatch( v, v_ans, colsize ) ) result1 = false;
  for( i=0; i<colsize; i++ ) v_ans[i] = i + 1;
  zRawMatGetRow( m, colcapacity, rowsize, colsize, 1, v );
  if( !zRawVecMatch( v, v_ans, colsize ) ) result1 = false;
  result2 = true;
  for( i=0; i<rowsize; i++ ) v_ans[i] = i;
  zRawMatGetCol( m, colcapacity, rowsize, colsize, 0, v );
  if( !zRawVecMatch( v, v_ans, rowsize ) ) result2 = false;
  for( i=0; i<rowsize; i++ ) v_ans[i] = i + 1;
  zRawMatGetCol( m, colcapacity, rowsize, colsize, 1, v );
  if( !zRawVecMatch( v, v_ans, rowsize ) ) result2 = false;
  result3 = true;
  for( i=0; i<colsize; i++ ) v[i] = 10.0;
  zRawMatPutRow( m, colcapacity, rowsize, colsize, 0, v );
  for( i=0; i<colsize; i++ )
    if( m[i] != 10 ) result3 = false;
  result4 = true;
  zRawMatPutCol( m, colcapacity, rowsize, colsize, 0, v );
  for( i=0; i<rowsize; i++ )
    if( m[i*colcapacity] != 10 ) result4 = false;

  free( m );
  free( v );
  free( v_ans );
  zAssert( zRawMatGetRow, result1 );
  zAssert( zRawMatGetCol, result2 );
  zAssert( zRawMatPutRow, result3 );
  zAssert( zRawMatPutCol, result4 );
}

void assert_raw_mat_swap_row_col(void)
{
  double *m, *v, *v_ans;
  const int rowcapacity = 6;
  const int colcapacity = 8;
  const int rowsize = 4;
  const int colsize = 5;
  int i, j;
  bool result1, result2;

  m = zAlloc( double, rowcapacity * colcapacity );
  v = zAlloc( double, zMax( rowcapacity, colcapacity ) );
  v_ans = zAlloc( double, zMax( rowcapacity, colcapacity ) );
  for( i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      m[i*colcapacity+j] = i + j;

  result1 = true;
  zRawMatSwapRow( m, colcapacity, rowsize, colsize, 0, 2 );
  for( i=0; i<colsize; i++ )
    if( !zEqual( m[i], i+2, zTOL ) || !zEqual( m[2*colcapacity+i], i, zTOL ) ) result1 = false;
  zRawMatSwapRow( m, colcapacity, rowsize, colsize, 0, 2 );
  for( i=0; i<colsize; i++ )
    if( !zEqual( m[i], i, zTOL ) || !zEqual( m[2*colcapacity+i], i+2, zTOL ) ) result1 = false;
  result2 = true;
  zRawMatSwapCol( m, colcapacity, rowsize, colsize, 0, 2 );
  for( i=0; i<rowsize; i++ )
    if( !zEqual( m[i*colcapacity], i+2, zTOL ) || !zEqual( m[i*colcapacity+2], i, zTOL ) ) result2 = false;
  zRawMatSwapCol( m, colcapacity, rowsize, colsize, 0, 2 );
  for( i=0; i<rowsize; i++ )
    if( !zEqual( m[i*colcapacity], i, zTOL ) || !zEqual( m[i*colcapacity+2], i+2, zTOL ) ) result2 = false;

  free( m );
  free( v );
  free( v_ans );
  zAssert( zRawMatSwapRow, result1 );
  zAssert( zRawMatSwapCol, result2 );
}

void assert_raw_mat_get_put(void)
{
  double *src, *dest, *ans;
  const int rowcapacity1 = 5;
  const int colcapacity1 = 7;
  const int rowsize1 = 3;
  const int colsize1 = 3;
  const int rowcapacity2 = 5;
  const int colcapacity2 = 3;
  const int rowsize2 = 4;
  const int colsize2 = 2;
  int i, j, k, l;
  bool result1, result2, result3, result4;

  src = zAlloc( double, rowcapacity1 * colcapacity1 );
  dest = zAlloc( double, rowcapacity2 * colcapacity2 );
  ans = zAlloc( double, rowcapacity2 * colcapacity2 );
  for( i=0; i<rowsize1; i++ )
    for( j=0; j<colsize1; j++ )
      src[i*colcapacity1+j] = i * 2 + j;

  for( result1=true, i=0; i<rowsize1; i++ )
    for( j=0; j<colsize1; j++ ){
      zRawMatZero( dest, colcapacity2, rowcapacity2, colcapacity2 );
      zRawMatGet( src, colcapacity1, rowsize1, colsize1, i, j, dest, colcapacity2, rowsize2, colsize2 );
      zRawMatZero( ans, colcapacity2, rowcapacity2, colcapacity2 );
      for( k=0; k<rowsize2; k++ )
        for( l=0; l<colsize2; l++ ){
          if( i + k >= rowsize1 ) continue;
          if( j + l >= colsize1 ) continue;
          ans[k*colcapacity2+l] = ( i + k ) * 2 + j + l;
        }
      if( !zRawMatMatch( dest, colcapacity2, ans, colcapacity2, rowsize2, colsize2 ) ) result1 = false;
  }
  for( result2=true, i=0; i<rowsize1; i++ )
    for( j=0; j<colsize1; j++ ){
      zRawMatZero( dest, colcapacity2, rowcapacity2, colcapacity2 );
      zRawMatTGet( src, colcapacity1, rowsize1, colsize1, i, j, dest, colcapacity2, rowsize2, colsize2 );
      zRawMatZero( ans, colcapacity2, rowcapacity2, colcapacity2 );
      for( k=0; k<colsize2; k++ )
        for( l=0; l<rowsize2; l++ ){
          if( i + k >= colsize1 ) continue;
          if( j + l >= rowsize1 ) continue;
          ans[l*colcapacity2+k] = ( i + k ) * 2 + j + l;
        }
      if( !zRawMatMatch( dest, colcapacity2, ans, colcapacity2, rowsize2, colsize2 ) ) result2 = false;
  }
  for( result3=true, i=0; i<rowsize2; i++ )
    for( j=0; j<colsize2; j++ ){
      zRawMatZero( dest, colcapacity2, rowcapacity2, colcapacity2 );
      zRawMatPut( dest, colcapacity2, rowsize2, colsize2, i, j, src, colcapacity1, rowsize1, colsize1 );
      zRawMatZero( ans, colcapacity2, rowcapacity2, colcapacity2 );
      for( k=0; k<rowsize1; k++ )
        for( l=0; l<colsize1; l++ ){
          if( i + k >= rowsize2 ) continue;
          if( j + l >= colsize2 ) continue;
          ans[(i+k)*colcapacity2+(j+l)] = k * 2 + l;
        }
      if( !zRawMatMatch( dest, colcapacity2, ans, colcapacity2, rowsize2, colsize2 ) ) result3 = false;
  }
  for( result4=true, i=0; i<rowsize2; i++ )
    for( j=0; j<colsize2; j++ ){
      zRawMatZero( dest, colcapacity2, rowcapacity2, colcapacity2 );
      zRawMatTPut( dest, colcapacity2, rowsize2, colsize2, i, j, src, colcapacity1, rowsize1, colsize1 );
      zRawMatZero( ans, colcapacity2, rowcapacity2, colcapacity2 );
      for( k=0; k<colsize1; k++ )
        for( l=0; l<rowsize1; l++ ){
          if( i + k >= rowsize2 ) continue;
          if( j + l >= colsize2 ) continue;
          ans[(i+k)*colcapacity2+(j+l)] = k + l * 2;
        }
      if( !zRawMatMatch( dest, colcapacity2, ans, colcapacity2, rowsize2, colsize2 ) ) result4 = false;
  }
  free( src );
  free( dest );
  free( ans );
  zAssert( zRawMatGet, result1 );
  zAssert( zRawMatTGet, result2 );
  zAssert( zRawMatPut, result3 );
  zAssert( zRawMatTPut, result4 );
}

void assert_raw_mat_arith(void)
{
  double *m1, *m2, *m, *mc;
  double c;
  const int rowcapacity1 = 3, colcapacity1 = 4;
  const int rowcapacity2 = 5, colcapacity2 = 5;
  const int rowcapacity3 = 4, colcapacity3 = 4;
  const int rowsize = 2, colsize = 3;
  const int n = 10;
  int i, j, k;
  bool result_add, result_sub, result_rev, result_mul, result_div, result_cat;
  bool result_add_drc, result_sub_drc, result_rev_drc, result_mul_drc, result_div_drc, result_cat_drc;

  m1 = zAlloc( double, rowcapacity1 * colcapacity1 );
  m2 = zAlloc( double, rowcapacity2 * colcapacity2 );
  m  = zAlloc( double, rowcapacity3 * colcapacity3 );
  mc = zAlloc( double, rowcapacity1 * colcapacity1 );
  result_add = result_sub = result_rev = result_mul = result_div = result_cat = true;
  result_add_drc = result_sub_drc = result_rev_drc = result_mul_drc = result_div_drc = result_cat_drc = true;
  for( k=0; k<n; k++ ){
    zRawMatRandUniform( m1, colcapacity1, -10, 10, rowsize, colsize );
    zRawMatRandUniform( m2, colcapacity2, -10, 10, rowsize, colsize );
    c = zRandF( -10, 10 );
    zRawMatAdd( m1, colcapacity1, m2, colcapacity2, m, colcapacity3, rowsize, colsize );
    zRawMatCopy( m1, colcapacity1, mc, colcapacity1, rowsize, colsize );
    zRawMatAddDRC( mc, colcapacity1, m2, colcapacity2, rowsize, colsize );
    for( i=0; i<rowsize; i++ )
      for( j=0; j<colsize; j++ ){
        if( !zIsTiny( m1[i*colcapacity1+j]+m2[i*colcapacity2+j]-m[i*colcapacity3+j] ) ) result_add = false;
        if( !zIsTiny( mc[i*colcapacity1+j]-m[i*colcapacity3+j] ) ) result_add_drc = false;
      }
    zRawMatSub( m1, colcapacity1, m2, colcapacity2, m, colcapacity3, rowsize, colsize );
    zRawMatCopy( m1, colcapacity1, mc, colcapacity1, rowsize, colsize );
    zRawMatSubDRC( mc, colcapacity1, m2, colcapacity2, rowsize, colsize );
    for( i=0; i<rowsize; i++ )
      for( j=0; j<colsize; j++ ){
        if( !zIsTiny( m1[i*colcapacity1+j]-m2[i*colcapacity2+j]-m[i*colcapacity3+j] ) ) result_sub = false;
        if( !zIsTiny( mc[i*colcapacity1+j]-m[i*colcapacity3+j] ) ) result_sub_drc = false;
      }
    zRawMatRev( m1, colcapacity1, m, colcapacity3, rowsize, colsize );
    zRawMatCopy( m1, colcapacity1, mc, colcapacity1, rowsize, colsize );
    zRawMatRevDRC( mc, colcapacity1, rowsize, colsize );
    for( i=0; i<rowsize; i++ )
      for( j=0; j<colsize; j++ ){
        if( !zIsTiny( m1[i*colcapacity1+j]+m[i*colcapacity3+j] ) ) result_rev = false;
        if( !zIsTiny( mc[i*colcapacity1+j]-m[i*colcapacity3+j] ) ) result_rev_drc = false;
      }
    zRawMatMul( m1, colcapacity1, c, m, colcapacity3, rowsize, colsize );
    zRawMatCopy( m1, colcapacity1, mc, colcapacity1, rowsize, colsize );
    zRawMatMulDRC( mc, colcapacity1, c, rowsize, colsize );
    for( i=0; i<rowsize; i++ )
      for( j=0; j<colsize; j++ ){
        if( !zIsTiny( m1[i*colcapacity1+j]*c-m[i*colcapacity3+j] ) ) result_mul = false;
        if( !zIsTiny( mc[i*colcapacity1+j]-m[i*colcapacity3+j] ) ) result_mul_drc = false;
      }
    zRawMatDiv( m1, colcapacity1, c, m, colcapacity3, rowsize, colsize );
    zRawMatCopy( m1, colcapacity1, mc, colcapacity1, rowsize, colsize );
    zRawMatDivDRC( mc, colcapacity1, c, rowsize, colsize );
    for( i=0; i<rowsize; i++ )
      for( j=0; j<colsize; j++ ){
        if( !zIsTiny( m1[i*colcapacity1+j]/c-m[i*colcapacity3+j] ) ) result_div = false;
        if( !zIsTiny( mc[i*colcapacity1+j]-m[i*colcapacity3+j] ) ) result_div_drc = false;
      }
    zRawMatCat( m1, colcapacity1, c, m2, colcapacity2, m, colcapacity3, rowsize, colsize );
    zRawMatCopy( m1, colcapacity1, mc, colcapacity1, rowsize, colsize );
    zRawMatCatDRC( mc, colcapacity1, c, m2, colcapacity2, rowsize, colsize );
    for( i=0; i<rowsize; i++ )
      for( j=0; j<colsize; j++ ){
        if( !zIsTiny( m1[i*colcapacity1+j]+c*m2[i*colcapacity2+j]-m[i*colcapacity3+j] ) ) result_cat = false;
        if( !zIsTiny( mc[i*colcapacity1+j]-m[i*colcapacity3+j] ) ) result_cat_drc = false;
      }
  }

  free( m1 );
  free( m2 );
  free( m );
  free( mc );
  zAssert( zRawMatAdd, result_add );
  zAssert( zRawMatSub, result_sub );
  zAssert( zRawMatRev, result_rev );
  zAssert( zRawMatMul, result_mul );
  zAssert( zRawMatDiv, result_div );
  zAssert( zRawMatCat, result_cat );
  zAssert( zRawMatAddDRC, result_add_drc );
  zAssert( zRawMatSubDRC, result_sub_drc );
  zAssert( zRawMatRevDRC, result_rev_drc );
  zAssert( zRawMatMulDRC, result_mul_drc );
  zAssert( zRawMatDivDRC, result_div_drc );
  zAssert( zRawMatCatDRC, result_cat_drc );
}

void assert_raw_mat_col_arith(void)
{
  double *morg, *m, *v, *colvec, c;
  const int rowcapacity = 4, colcapacity = 5;
  const int rowsize = 3, colsize = 4;
  const int n = 10;
  int i, j, k, col;
  bool result_add, result_sub, result_mul, result_cat, result_ip;

  morg = zAlloc( double, rowcapacity * colcapacity );
  m = zAlloc( double, rowcapacity * colcapacity );
  v = zAlloc( double, rowcapacity );
  colvec = zAlloc( double, rowcapacity );
  result_add = result_sub = result_mul = result_cat = result_ip = true;
  for( k=0; k<n; k++ ){
    zRawMatRandUniform( morg, colcapacity, rowsize, colsize, -10, 10 );
    zRawVecRandUniform( v, rowsize, -10, 10 );
    c = zRandF( -10, 10 );
    zRawMatCopy( morg, colcapacity, m, colcapacity, rowsize, colsize );
    col = zRandI( 0, colsize-1 );
    zRawMatColAddDRC( m, colcapacity, v, rowsize, colsize, col );
    for( i=0; i<rowsize; i++ )
      for( j=0; j<colsize; j++ )
        if( j == col ){
          if( !zIsTiny( morg[i*colcapacity+j]+v[i]-m[i*colcapacity+j] ) ) result_add = false;
        } else{
          if( !zIsTiny( morg[i*colcapacity+j]-m[i*colcapacity+j] ) ) result_add = false;
        }
    zRawMatColSubDRC( m, colcapacity, v, rowsize, colsize, col );
    if( !zRawMatEqual( m, colcapacity, morg, colcapacity, rowsize, colsize, zTOL ) ) result_sub = false;
    zRawMatColMulDRC( m, colcapacity, c, rowsize, colsize, col );
    for( i=0; i<rowsize; i++ )
      for( j=0; j<colsize; j++ )
        if( j == col ){
          if( !zIsTiny( morg[i*colcapacity+j]*c-m[i*colcapacity+j] ) ) result_mul = false;
        } else{
          if( !zIsTiny( morg[i*colcapacity+j]-m[i*colcapacity+j] ) ) result_mul = false;
        }
    zRawMatCopy( morg, colcapacity, m, colcapacity, rowsize, colsize );
    zRawMatColCatDRC( m, colcapacity, c, v, rowsize, colsize, col );
    for( i=0; i<rowsize; i++ )
      for( j=0; j<colsize; j++ )
        if( j == col ){
          if( !zIsTiny( morg[i*colcapacity+j]+c*v[i]-m[i*colcapacity+j] ) ) result_cat = false;
        } else{
          if( !zIsTiny( morg[i*colcapacity+j]-m[i*colcapacity+j] ) ) result_cat = false;
        }
    zRawMatGetCol( m, colcapacity, rowsize, colsize, col, colvec );
    if( !zEqual( zRawMatColInnerProd( m, colcapacity, v, rowsize, colsize, col ), zRawVecInnerProd( colvec, v, rowsize ), zTOL ) ) result_ip = false;
  }
  free( morg );
  free( m );
  free( v );
  free( colvec );
  zAssert( zRawMatColAddDRC, result_add );
  zAssert( zRawMatColSubDRC, result_sub );
  zAssert( zRawMatColMulDRC, result_mul );
  zAssert( zRawMatColCatDRC, result_cat );
  zAssert( zRawMatColInnerProd, result_ip );
}

void assert_raw_mat_norm(void)
{
  double m[] = {
    1.0, 2.0, 3.0, 10.0,
    4.0, 5.0, 6.0, 10.0,
   10.0,10.0,10.0, 10.0,
  };

  zAssert( zRawMatSqrNorm, zEqual( zRawMatSqrNorm( m, 4, 2, 3 ), 1.0*1.0+2.0*2.0+3.0*3.0+4.0*4.0+5.0*5.0+6.0*6.0, zTOL ) );
}

void assert_raw_mat_transpose(void)
{
  double *m, *tm, tr1, tr2;;
  const int rowcapacity = 10, colcapacity = 10;
  int rowsize, colsize;
  const int n = 10;
  int i, j, k;
  bool result_t, result_t_drc, result_tr;

  m = zAlloc( double, rowcapacity * colcapacity );
  tm = zAlloc( double, rowcapacity * colcapacity );
  result_t = result_t_drc = true;
  for( k=0; k<n; k++ ){
    rowsize = zRandI( 1, rowcapacity );
    colsize = zRandI( 1, colcapacity );
    zRawMatZero( m, colcapacity, rowcapacity, colcapacity );
    zRawMatRandUniform( m, colcapacity, -10, 10, rowsize, colsize );
    zRawMatT( m, colcapacity, tm, rowcapacity, rowsize, colsize );
    for( i=0; i<rowsize; i++ )
      for( j=0; j<colsize; j++ )
        if( m[i*colcapacity+j] != tm[j*rowcapacity+i] ) result_t = false;
    zRawMatTDRC( tm, colcapacity, rowcapacity );
    if( !zRawMatMatch( m, colcapacity, tm, colcapacity, rowsize, colsize ) ) result_t_drc = false;
  }
  tr1 = zRawMatTrace( m, colcapacity, rowsize, colsize );
  for( tr2=0, i=0; i<rowsize && i<colsize; i++ )
    tr2 += m[i*colcapacity+i];
  result_tr = zEqual( tr1, tr2, zTOL );

  free( m );
  free( tm );
  zAssert( zRawMatT, result_t );
  zAssert( zRawMatTDRC, result_t_drc );
  zAssert( zRawMatTrace, result_tr );
}

void assert_raw_mul_mat_vec(void)
{
  double m[] = {
    1.0, 2.0, 3.0, 0.0,
    2.0, 3.0, 4.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
  };
  double v[] = {
    3.0, 4.0, 5.0,
  };
  double mv[3];
  double mv_ans[] = {
    26.0, 38.0,
  };
  double mtv_ans[] = {
    11.0, 18.0, 25.0,
  };
  const int colcapacity = 4;
  const int rowsize = 2, colsize = 3;

  zRawMulMatVec( m, colcapacity, v, rowsize, colsize, mv );
  zAssert( zRawMulMatVec, zRawVecEqual( mv, mv_ans, rowsize, zTOL ) );
  zRawMulMatTVec( m, colcapacity, v, rowsize, colsize, mv );
  zAssert( zRawMulMatTVec, zRawVecEqual( mv, mtv_ans, colsize, zTOL ) );
}

void assert_raw_mul_mat_mat(void)
{
  double m1[] = {
   -1.0, 2.0, 3.0, 0.0,
    2.0,-3.0, 4.0, 0.0,
    3.0,-4.0, 5.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
  };
  double m2[] = {
    1.0,-2.0,-3.0, 0.0,
   -2.0,-3.0,-4.0, 0.0,
    3.0, 4.0,-5.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
  };
  double m[15];
  double mm_ans[] = { /* (2x3).(3x2) */
    4.0, 8.0,
   20.0,21.0,
  };
  double mmt_ans[] = { /* (2x3).(2x3)^T */
   -14.0, -16.0,
    -4.0, -11.0,
  };
  double mtm_ans[] = { /* (2x3)^T .(2x3) */
   -5.0, -4.0, -5.0,
    8.0,  5.0,  6.0,
   -5.0,-18.0,-25.0,
  };

  zRawMulMatMat( m1, 4, 2, 3, m2, 4, 3, 2, m, 5 );
  zAssert( zRawMulMatMat, zRawMatEqual( m, 5, mm_ans, 2, 2, 2, zTOL ) );
  zRawMulMatMatT( m1, 4, 2, 3, m2, 4, 2, 3, m, 5 );
  zAssert( zRawMulMatMatT, zRawMatEqual( m, 5, mmt_ans, 2, 2, 2, zTOL ) );
  zRawMulMatTMat( m1, 4, 2, 3, m2, 4, 2, 3, m, 5 );
  zAssert( zRawMulMatTMat, zRawMatEqual( m, 5, mtm_ans, 3, 3, 3, zTOL ) );
}

void assert_raw_mat_dyad(void)
{
  double v1[] = { 1.0, 2.0, 3.0 };
  double v2[] = { 1.0, 3.0, 5.0 };
  double m[20];
  double m_ans[] = {
    1.0, 3.0, 5.0,
    2.0, 6.0,10.0,
    3.0, 9.0,15.0,
  };
  double m_add_ans[] = {
    2.0, 4.0, 6.0,
    3.0, 7.0,11.0,
    4.0,10.0,16.0,
  };
  double m_cat_ans[] = {
    3.0, 7.0,11.0,
    5.0,13.0,21.0,
    7.0,19.0,31.0,
  };

  zRawMatZero( m, 5, 3, 5 );
  zRawVecDyad( v1, 3, v2, 3, m, 5 );
  zAssert( zRawVecDyad, zRawMatEqual( m, 5, m_ans, 3, 3, 3, zTOL ) );
  set_raw_mat_all( m, 5, 3, 3, 1.0 );
  zRawMatAddDyad( m, 5, v1, 3, v2, 3 );
  zAssert( zRawMatAddDyad, zRawMatEqual( m, 5, m_add_ans, 3, 3, 3, zTOL ) );
  zRawMatCopy( m_ans, 3, m, 5, 3, 3 );
  zRawMatSubDyad( m, 5, v1, 3, v2, 3 );
  zAssert( zRawMatSubDyad, zRawMatIsTiny( m, 5, 3, 3 ) );
  set_raw_mat_all( m, 5, 3, 3, 1.0 );
  zRawMatCatDyad( m, 5, 2.0, v1, 3, v2, 3 );
  zAssert( zRawMatCatDyad, zRawMatEqual( m, 5, m_cat_ans, 3, 3, 3, zTOL ) );
}

int main(void)
{
  zRandInit();
  assert_raw_mat_zero();
  assert_raw_mat_touchup();
  assert_raw_mat_rand();
  assert_raw_mat_copy();
  assert_raw_mat_equal();
  assert_raw_mat_is_tol();
  assert_raw_mat_diag();
  assert_raw_mat_get_put_row_col();
  assert_raw_mat_swap_row_col();
  assert_raw_mat_get_put();
  assert_raw_mat_arith();
  assert_raw_mat_col_arith();
  assert_raw_mat_norm();
  assert_raw_mat_transpose();
  assert_raw_mul_mat_vec();
  assert_raw_mul_mat_mat();
  assert_raw_mat_dyad();
  return 0;
}
