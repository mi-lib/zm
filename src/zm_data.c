/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_data - data analysis.
 */

#include <zm/zm_data.h>
#include <zm/zm_ip.h>

/* internal function for smoothing a data sequence based on Savitzky-Golay's method */
static bool _zDataSmoothSG(double src[], size_t n, size_t w, int dim, double dest[], double (*pexip_func)(zPexIP*,double))
{
  zPexIP pc;
  zVec ts; /* timestamp */
  zVecStruct window;
  int wh, nw;
  register int i;
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
bool zDataSmoothSG(double src[], size_t n, size_t w, int dim, double dest[])
{
  return _zDataSmoothSG( src, n, w, dim, dest, zPexIPVal );
}

/* smooth a data sequence based on Savitzky-Golay's method */
bool zDataSmoothVelSG(double src[], size_t n, size_t w, int dim, double dest[])
{
  return _zDataSmoothSG( src, n, w, dim, dest, zPexIPVel );
}

zListClass( int_list_t, int_list_cell_t, int );

/* smooth a data sequence based on Savitzky-Golay's method and pick up peaks. */
zIndex zDataPeakSG(double src[], size_t n, int w, int dim)
{
  zIndex peakidx = NULL;
  int_list_t list;
  int_list_cell_t *cp;
  double *g;
  register int i;

  if( !( g = zAlloc( double, n ) ) ) return NULL;
  if( !zDataSmoothVelSG( src, n, w, dim, g ) ) goto TERMINATE;
  zListInit( &list );
  n--;
  for( i=1; i<n; i++ ){
    if( g[i-1] > 0 && g[i] < 0 ){
      if( !( cp = zAlloc( int_list_cell_t, 1 ) ) ) goto TERMINATE;
      cp->data = i;
      zListInsertHead( &list, cp );
    }
  }
  if( zListSize(&list) > 0 ){
    if( !( peakidx = zIndexCreate( zListSize(&list) ) ) ) goto TERMINATE;
    i=0;
    zListForEach( &list, cp )
      zIndexElemNC(peakidx,i++) = cp->data;
  }
 TERMINATE:
  free( g );
  zListDestroy( int_list_cell_t, &list );
  return peakidx;
}
