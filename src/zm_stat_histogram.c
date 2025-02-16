/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_stat_histogram.h
 * \brief statistics: histogram.
 * \author Zhidao
 */

#include <zm/zm_stat.h>

/* allocate bin of a histogram. */
zHistogram *zHistogramAllocBin(zHistogram *histogram, double min, double max, size_t num)
{
  if( !( histogram->bin = zIndexAlloc( num ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  zHistogramBinClear( histogram );
  histogram->min = min;
  histogram->max = max;
  return histogram;
}

/* destroy a histogram. */
void zHistogramDestroy(zHistogram *histogram)
{
  zIndexFree( histogram->bin );
}

/* vote a data and increment bin of a histogram. */
void zHistogramVote(zHistogram *histogram, double val)
{
  int i;

  i = ( val - histogram->min ) / ( histogram->max - histogram->min ) * zIndexSizeNC(histogram->bin);
  i = zLimit( i, 0, zIndexSizeNC(histogram->bin)-1 );
  zIndexElemNC(histogram->bin,i)++;
}

/* create a histogram from data and specified range. */
zHistogram *zHistogramCreate(zHistogram *histogram, const double data[], size_t size, double min, double max, size_t num)
{
  int i;

  if( !zHistogramAllocBin( histogram, min, max, num ) )
    return NULL;
  for( i=0; i<size; i++ )
    zHistogramVote( histogram, data[i] );
  return histogram;
}

/* create a histogram from data by automatically adjust range. */
zHistogram *zHistogramCreateAuto(zHistogram *histogram, const double data[], size_t size, int num)
{
  double min, max, margin;

  zDataMinMax( data, size, &min, NULL, &max, NULL );
  margin = 0.5 * ( max - min ) / ( num - 1 );
  return zHistogramCreate( histogram, data, size, min - margin, max + margin, num );
}

/* minimum value of the range of a bin. */
double zHistogramBinMin(zHistogram *histogram, int i)
{
  i = zLimit( i, 0, zIndexSizeNC(histogram->bin)-1 );
  return histogram->min + ( histogram->max - histogram->min ) * (double)i/zIndexSizeNC(histogram->bin);
}

/* maximum value of the range of a bin. */
double zHistogramBinMax(zHistogram *histogram, int i)
{
  return zHistogramBinMin( histogram, i+1 );
}

/* middle value of the range of a bin. */
double zHistogramBinMid(zHistogram *histogram, int i)
{
  return 0.5 * ( zHistogramBinMax(histogram,i) + zHistogramBinMin(histogram,i) );
}

/* print out a histogram to a file. */
void zHistogramFPrint(FILE *fp, zHistogram *histogram)
{
  int i;

  for( i=0; i<zIndexSizeNC(histogram->bin); i++ ){
    fprintf( fp, "%g %d\n", zHistogramBinMin(histogram,i), zHistogramBinCount(histogram,i) );
  }
}
