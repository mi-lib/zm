/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_pso - optimization tools: Particle Swarm Optimization method.
 */

#ifndef __ZM_OPT_PSO_H__
#define __ZM_OPT_PSO_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* Particle Swarm Optimizaiton method by Kennedy and Eberhart (1995).
 * reference:
 * J. Kennedy and R. Eberhart, Particle swarm optimization, in Proceedings of
 * IEEE International Conference on Neural Networks, pp. 1942-1948, 1995.
 */

#define ZOPT_PSO_DEFAULT_NUM      100
#define ZOPT_PSO_DEFAULT_C1         2.0
#define ZOPT_PSO_DEFAULT_C2         2.0
#define ZOPT_PSO_DEFAULT_VEL_RATE   0.1

/*! \brief solve an optimization problem by Particle Swarm Optimization method.
 *
 * zOptSolvePSO() solves an optimization problem by Particle Swarm Optimization
 * method proposed by Kennedy and Eberhart (1995).
 */
__ZM_EXPORT int zOptSolvePSO(double (* f)(const zVec,void*), void *util, const zVec min, const zVec max, int iter, double tol, int num, double c1, double c2, double vel_rate, zVec ans, double *eval);

/*! \brief solve an optimization problem by Particle Swarm Optimization method.
 *
 * zOptSolvePSODefault() solves an optimization problem by Particle Swarm
 * Optimization method with default parameters.
 */
__ZM_EXPORT int zOptSolvePSODefault(double (* f)(const zVec,void*), void *util, const zVec min, const zVec max, int iter, double tol, zVec ans, double *eval);

__END_DECLS

#endif /* __ZM_OPT_PSO_H__ */
