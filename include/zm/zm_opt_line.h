/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_line - optimization tools: line search.
 */

#ifndef __ZM_OPT_LINE_H__
#define __ZM_OPT_LINE_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zOptLine
 * line search class
 * ********************************************************** */

/*! \brief golden section proportion */
#define Z_GSEC 0.38196601125911

/*! \brief golden section method */
__ZM_EXPORT double zOptLineGSEC(double (*eval)(double,void*), double a, double b, void *util, int iter);

/*! \brief bisection method */
__ZM_EXPORT double zOptLineBisec(double (*eval)(double,void*), double a, double b, void *util, int iter);

/*! \brief Brent's method */
__ZM_EXPORT double zOptLineBrent(double (*eval)(double,void*), double a, double b, void *util, int iter);

__END_DECLS

#endif /* __ZM_OPT_LINE_H__ */
