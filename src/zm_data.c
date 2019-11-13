/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_data - data analysis.
 */

#include <zm/zm_data.h>
#include <zm/zm_ip.h>

/* smooth a data sequence based on Savitzky-Golay's method */
void zDataSmoothSG(double src[], size_t n, size_t w, int dim, double dest[])
{
  zPexIP pc;
  zVec ts; /* timestamp */
  zVecStruct window;
  int wh, nw;
  register int i;

  ts = zVecAlloc( w );
  if( !zPexIPAlloc( &pc, w-1, dim ) || !ts ) goto TERMINATE;
  for( i=0; i<w; i++ ) zVecElemNC(ts,i) = i;
  zVecSetSizeNC( &window, w );
  zVecBufNC(&window) = src;
  if( !zPexIPCreateLSM( &pc, w-1, dim, ts, &window ) ) goto TERMINATE;
  wh = w / 2;
  for( i=0; i<wh; i++ ) dest[i] = zPexIPVal(&pc,i);
  nw = n - wh;
  for( ; i<nw; i++ ){
    zVecBufNC(&window)++;
    if( !zPexIPCreateLSM( &pc, w-1, dim, ts, &window ) ) goto TERMINATE;
    dest[i] = zPexIPVal(&pc,wh-1);
  }
  for( ; i<n; i++ ) dest[i] = zPexIPVal(&pc,i-nw+wh);

 TERMINATE:
  zVecFree( ts );
  zPexIPFree( &pc );
}
