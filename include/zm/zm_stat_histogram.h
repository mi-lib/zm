/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_stat_histogram.h
 * \brief statistics: histogram.
 * \author Zhidao
 */

#ifndef __ZM_STAT_HISTOGRAM_H__
#define __ZM_STAT_HISTOGRAM_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \struct \brief histogram class. */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zHistogram ){
  double min; /*!< \brief minimum value */
  double max; /*!< \brief maximum value */
  zIndex bin; /*!< \brief bin */
};

#define zHistogramBinSize(histogram)    zIndexSize( (histogram)->bin )
#define zHistogramBinCount(histogram,i) zIndexElem( (histogram)->bin, i )
#define zHistogramBinClear(histogram)   zIndexZero( (histogram)->bin )

/*! \brief allocate bin of a histogram. */
__ZM_EXPORT zHistogram *zHistogramAllocBin(zHistogram *histogram, double min, double max, size_t num);
/*! \brief destroy a histogram. */
__ZM_EXPORT void zHistogramDestroy(zHistogram *histogram);
/*! \brief vote a data and increment bin of a histogram. */
__ZM_EXPORT void zHistogramVote(zHistogram *histogram, double val);

/*! \brief create a histogram from data and specified range. */
__ZM_EXPORT zHistogram *zHistogramCreate(zHistogram *histogram, const double data[], size_t size, double min, double max, size_t num);
/*! \brief create a histogram from data by automatically adjust range. */
__ZM_EXPORT zHistogram *zHistogramCreateAuto(zHistogram *histogram, const double data[], size_t size, int num);

/* minimum value of the range of a bin. */
__ZM_EXPORT double zHistogramBinMin(zHistogram *histogram, int i);
/* maximum value of the range of a bin. */
__ZM_EXPORT double zHistogramBinMax(zHistogram *histogram, int i);
/* middle value of the range of a bin. */
__ZM_EXPORT double zHistogramBinMid(zHistogram *histogram, int i);

/*! \brief print out a histogram to a file. */
__ZM_EXPORT void zHistogramFPrint(FILE *fp, zHistogram *histogram);

__END_DECLS

#endif /* __ZM_STAT_HISTOGRAM_H__ */
