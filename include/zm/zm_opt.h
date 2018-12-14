/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt - optimization tools.
 */

#ifndef __ZM_OPT_H__
#define __ZM_OPT_H__

#include <zm/zm_le.h>

#include <zm/zm_opt_line.h> /* line search */
#include <zm/zm_opt_lp.h>   /* linear programming */
#include <zm/zm_opt_lcp.h>  /* linear complementary problem */
#include <zm/zm_opt_qp.h>   /* quadratic programming */

/*! \macro small quantity for finite-difference derivative */
#define Z_OPT_EPS ( 1.0e-6 )

#include <zm/zm_opt_nm.h>   /* downhill simplex/polytope method */
#include <zm/zm_opt_ga.h>   /* genetic algorithm */
#include <zm/zm_opt_dm.h>   /* descent method */

#endif /* __ZM_OPT_H__ */
