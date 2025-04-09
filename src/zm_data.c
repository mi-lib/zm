/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_data - data analysis.
 */

#include <zm/zm_data.h>
#include <zm/zm_ip.h>

/* check if the specified point is a peak. */
static bool _zDataIsPeak(double src[], size_t n, uint i)
{
  if( i <= 0 || i >= n - 1 ) return false;
  return src[i] > src[i-1] && src[i] > src[i+1];
}

/* pick up peaks of a data sequence. */
zIndex zDataPeak(double src[], size_t n, uint w)
{
  zIndex peakidx = NULL;
  zIntList list;
  zIntListCell *cp;
  uint i;

  zListInit( &list );
  for( i=1; i<n; i++ ){
    if( _zDataIsPeak( src, n, i ) ){
      if( zListIsEmpty(&list) || i > zListHead(&list)->data + w ){
        if( !zIntListAdd( &list, i ) ) goto TERMINATE;
      } else
      if( src[i] > src[zListHead(&list)->data] ){
        zListDeleteHead( &list, &cp );
        free( cp );
        if( !zIntListAdd( &list, i ) ) goto TERMINATE;
      }
    }
  }
  peakidx = zIndexCreateFromList( &list );
 TERMINATE:
  zListDestroy( zIntListCell, &list );
  return peakidx;
}

/* an internal comparison function for zDataSort. */
static int _zDataSortCmp(void *v1, void *v2, void *dummy)
{
  if( *(double*)v1 > *(double*)v2 ) return 1;
  if( *(double*)v1 < *(double*)v2 ) return -1;
  return 0;
}

/* sort an array of double-precision floating-point vector in the descent order. */
double *zDataSort(double data[], size_t n)
{
  zQuickSort( data, n, sizeof(double), _zDataSortCmp, NULL );
  return data;
}

/* an internal comparison function for zDataSortIndex. */
static int _zDataSortIndexCmp(void *i1, void *i2, void *priv)
{
  if( ((double *)priv)[*(int*)i1] < ((double *)priv)[*(int*)i2] ) return 1;
  if( ((double *)priv)[*(int*)i1] > ((double *)priv)[*(int*)i2] ) return -1;
  return 0;
}

/* sort an integer vector in the descent order of the corresponding samples. */
zIndex zDataSortIndex(double data[], size_t n, zIndex index)
{
  zQuickSort( zArrayBuf(index), zIndexSizeNC(index), sizeof(int), _zDataSortIndexCmp, data );
  return index;
}

/* an internal comparison function for zDataSortAbsIndex. */
static int _zDataSortAbsIndexCmp(void *i1, void *i2, void *priv)
{
  double d1, d2;

  d1 = fabs( ((double *)priv)[*(int*)i1] );
  d2 = fabs( ((double *)priv)[*(int*)i2] );
  if( d1 < d2 ) return 1;
  if( d1 > d2 ) return -1;
  return 0;
}

/* sort an integer vector in the descent order of absolute values of the corresponding samples. */
zIndex zDataSortAbsIndex(double data[], size_t n, zIndex index)
{
  zQuickSort( zArrayBuf(index), zIndexSizeNC(index), sizeof(int), _zDataSortAbsIndexCmp, data );
  return index;
}

/* internal function for smoothing a data sequence based on Savitzky-Golay's method */
static bool _zDataSmoothSG(double src[], size_t n, size_t w, uint dim, double dest[], double (*pexip_func)(const zPexIP*,double))
{
  zPexIP pc;
  zVec ts; /* timestamp */
  zVecStruct window;
  uint wh, nw;
  uint i;
  bool result = true;

  ts = zVecAlloc( w );
  if( !zPexIPAlloc( &pc, w-1, dim ) || !ts ){
    result = false;
    goto TERMINATE;
  }
  for( i=0; i<w; i++ ) zVecElemNC(ts,i) = i;
  zVecSetSizeNC( &window, w );
  zVecBufNC(&window) = src;
  if( !zPexIPCreateLSM( &pc, w-1, dim, ts, &window ) ){
    result = false;
    goto TERMINATE;
  }
  wh = w / 2;
  for( i=0; i<wh; i++ ) dest[i] = pexip_func(&pc,i);
  nw = n - wh;
  for( ; i<nw; i++ ){
    zVecBufNC(&window)++;
    if( !zPexIPCreateLSM( &pc, w-1, dim, ts, &window ) ){
      result = false;
      goto TERMINATE;
    }
    dest[i] = pexip_func(&pc,wh-1);
  }
  for( ; i<n; i++ ) dest[i] = pexip_func(&pc,i-nw+wh);

 TERMINATE:
  zVecFree( ts );
  zPexIPFree( &pc );
  return result;
}

/* smooth a data sequence based on Savitzky-Golay's method */
bool zDataSmoothSG(double src[], size_t n, size_t w, uint dim, double dest[])
{
  return _zDataSmoothSG( src, n, w, dim, dest, zPexIPVal );
}

/* smooth a data sequence based on Savitzky-Golay's method */
bool zDataSmoothVelSG(double src[], size_t n, size_t w, uint dim, double dest[])
{
  return _zDataSmoothSG( src, n, w, dim, dest, zPexIPVel );
}

/* pick up peaks of a smoothed data sequence based on Savitzky-Golay's method. */
zIndex zDataPeakSG(double src[], size_t n, size_t w, uint dim)
{
  zIndex peakidx = NULL;
  zIntList list;
  double *g;
  uint i;

  if( !( g = zAlloc( double, n ) ) ) return NULL;
  if( !zDataSmoothVelSG( src, n, w, dim, g ) ) goto TERMINATE;
  zListInit( &list );
  n--;
  for( i=1; i<n; i++ )
    if( g[i-1] > 0 && g[i] < 0 )
      if( !zIntListAdd( &list, i ) ) goto TERMINATE;
  peakidx = zIndexCreateFromList( &list );
 TERMINATE:
  free( g );
  zListDestroy( zIntListCell, &list );
  return peakidx;
}
