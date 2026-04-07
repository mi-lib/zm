/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_line - optimization tools: line search.
 */

#ifndef __ZM_OPT_LINE_H__
#define __ZM_OPT_LINE_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief golden ratio (sqrt(5)+1)/2 */
#define zGOLDENRATIO  1.618033988749895
/*! \brief minor segment ratio of golden section (sqrt(5)+1)/(sqrt(5)-1) */
#define zGOLDENRATIO2 0.381966011250105

/*! \brief golden section method */
__ZM_EXPORT double zOptLineGoldenSection(double (*eval)(double,void*), double xmin, double xmax, void *util, int iter);

/*! \brief trisection method. */
__ZM_EXPORT double zOptLineTrisection(double (*eval)(double,void*), double xmin, double xmax, void *util, int iter);

/*! \brief Brent's method */
__ZM_EXPORT double zOptLineBrent(double (*eval)(double,void*), double a, double b, void *util, int iter);

__END_DECLS

#endif /* __ZM_OPT_LINE_H__ */
