/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_nm - optimization tools: Nelder-Mead's downhill simplex method.
 */

#ifndef __ZM_OPT_NM_H__
#define __ZM_OPT_NM_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

typedef struct{
  double (*eval)(zVec,void*); /* evaluator function */
  int num;   /* number of vertices of simplex */
  zVec *e;   /* bases */
  zVec f;    /* evaluation values */
  zVec pin;  /* pin (COG of the bottom face of simplex) */
  zVec test; /* test stick */
  zIndex index;
} zOptNM;

/* downhill simplex/polytope method by Nelder and Mead (1965).
 */
__EXPORT zOptNM *zOptNMCreate(zOptNM *opt, int dim, double (*eval)(zVec,void*));
__EXPORT void zOptNMDestroy(zOptNM *opt);
__EXPORT int zOptNMSolve(zOptNM *opt, zVec var, void *util, double tol, int iter, double *eval);

__END_DECLS

#endif /* __ZM_OPT_NM_H__ */
