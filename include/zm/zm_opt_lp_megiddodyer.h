/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_lp_megiddodyer - optimization tools: linear programming based on Megiddo-Dyer algorithm.
 */

#ifndef __ZM_OPT_LP_MEGIDDO_DYER_H__
#define __ZM_OPT_LP_MEGIDDO_DYER_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \struct zLPMegiddoDyerConstraint
 * \brief linear inequality constraint: a1 x1 + a2 x2 <= b -> y ><= m x + y0
 */
typedef struct{
  double m;
  double y0;
} zLPMegiddoDyerConstraint;

/*! \struct zLPMegiddoDyerConstraintList
 * \brief list of inequality constraints.
 */
ZEDA_DEF_LIST_CLASS( zLPMegiddoDyerConstraintList, zLPMegiddoDyerConstraintListCell, zLPMegiddoDyerConstraint );

/*! \struct zLPMegiddoDyerConstraintPtrList
 * \brief list pointers to inequality constraints.
 */
ZEDA_DEF_LIST_CLASS( zLPMegiddoDyerConstraintPtrList, zLPMegiddoDyerConstraintPtrListCell, zLPMegiddoDyerConstraintListCell* );

/*! \struct zLPMegiddoDyerConstraintPair
 * \brief a pair of linear inequality constraints.
 */
typedef struct{
  zLPMegiddoDyerConstraintPtrListCell *cp1;
  zLPMegiddoDyerConstraintPtrListCell *cp2;
  double x;
  double y;
} zLPMegiddoDyerConstraintPair;

/*! \struct zLPMegiddoDyerConstraintPairList
 * \brief list of pairs of linear inequality constraints.
 */
ZEDA_DEF_LIST_CLASS( zLPMegiddoDyerConstraintPairList, zLPMegiddoDyerConstraintPairListCell, zLPMegiddoDyerConstraintPair );

/*! \struct zLPMegiddoDyer
 * \brief a linear programming problem solver by Megiddo-Dyer algorithm.
 *
 * It solves a linear programming problem with two variables (x1,x2) in the following form:
 *  min_(x1,x2) { c1 x1 + c2 x2 | ak1 x1 + ak2 x2 <= bk (k=1,...,N) }
 * A typicall flow to define and solve a linear programming problem is as follows:
 *  zLPMegiddoDyerInit( &lp );
 *  zLPMegiddoDyerSetCostCoefficient( &lp, c1, c2 );
 *  zLPMegiddoDyerAddConstraint( &lp, a1, a2, b );
 *  ...
 *  zLPMegiddoDyerSolve( &lp, &x1, &x2, &cost );
 *  zLPMegiddoDyerDestroy( &lp );
 */
typedef struct{
  double cost_c1; /*!< \brief a coefficient of the first variable of cost function */
  double cost_c2; /*!< \brief a coefficient of the second variable of cost function */
  double x_min;   /*!< \brief minimum bound of x */
  double x_max;   /*!< \brief maximum bound of x */
  /*! \cond */
  zLPMegiddoDyerConstraintList upperbound;
  zLPMegiddoDyerConstraintList lowerbound;
  zLPMegiddoDyerConstraintPtrList upperboundpool;
  zLPMegiddoDyerConstraintPtrList lowerboundpool;
  zLPMegiddoDyerConstraintPairList upperboundpair;
  zLPMegiddoDyerConstraintPairList lowerboundpair;
  /* to check boundedness of feasible region */
  double mu_min; /* minimum gradient of upperbound segments */
  double mu_max; /* maximum gradient of upperbound segments */
  double ml_min; /* minimum gradient of lowerbound segments */
  double ml_max; /* maximum gradient of lowerbound segments */
  ubyte feasibility_flag;
  /*! \endcond */
} zLPMegiddoDyer;

/*! \brief initialize properties for Megiddo-Dyer algorithm. */
__ZM_EXPORT void zLPMegiddoDyerInit(zLPMegiddoDyer *lp);

/*! \brief destroy properties for Megiddo-Dyer algorithm. */
__ZM_EXPORT void zLPMegiddoDyerDestroy(zLPMegiddoDyer *lp);

/*! \brief set coefficients of cost function for Megiddo-Dyer algorithm. */
__ZM_EXPORT bool zLPMegiddoDyerSetCostCoefficient(zLPMegiddoDyer *lp, double c1, double c2);

/*! \brief add a linear inequality constraint to a linear programming problem to be solved by Megiddo-Dyer algorithm. */
__ZM_EXPORT bool zLPMegiddoDyerAddConstraint(zLPMegiddoDyer *lp, double a1, double a2, double b);

/*! \brief solve a linear programming problem by Megiddo-Dyer algorithm. */
__ZM_EXPORT bool zLPMegiddoDyerSolve(zLPMegiddoDyer *lp, double *x1, double *x2, double *cost);

/*! \brief solve a linear programming problem directly by checking all intersections of boundaries of linear inequality constraints. */
__ZM_EXPORT bool zLPMegiddoDyerSolveDirect(zLPMegiddoDyer *lp, double *x1, double *x2, double *cost);

#if DEBUG == 1
/*! \brief read a linear programming problem from file for Megiddo-Dyer algorithm. */
__ZM_EXPORT zLPMegiddoDyer *zLPMegiddoDyerReadFile(zLPMegiddoDyer *lp, const char *filename);
/*! \brief write a linear programming problem to a file for Megiddo-Dyer algorithm. */
__ZM_EXPORT zLPMegiddoDyer *zLPMegiddoDyerWriteFile(zLPMegiddoDyer *lp, const char *filename);
/*! \brief write a linear programming problem to a file in an ASCII format for Megiddo-Dyer algorithm. */
__ZM_EXPORT zLPMegiddoDyer *zLPMegiddoDyerWriteFileASCII(zLPMegiddoDyer *lp, const char *filename);
#endif /* DEBUG */

__END_DECLS

#endif /* __ZM_OPT_LP_MEGIDDO_DYER_H__ */
