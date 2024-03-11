/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_gng.h
 * \brief Growing Neural Gas.
 * \author Zhidao
 */

#ifndef __ZM_GNG_H__
#define __ZM_GNG_H__

#include <zm/zm_vec.h>

__BEGIN_DECLS

/* ********************************************************** */
/*! \defgroup GNG
 * \{ *//* ************************************************** */

struct _zGNGUnit;
typedef struct _zGNGUnit zGNGUnit;

/*! \brief edge class */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zGNGEdgeData ){
  int age;
  zGNGUnit *unit1;
  zGNGUnit *unit2;
};

zListClass( zGNGEdgeList, zGNGEdge, zGNGEdgeData );
zListClass( zGNGEdgePtrList, zGNGEdgePtr, zGNGEdgeData* );

/*! \brief unit class */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zGNGUnitData ){
  zVec v;
  double error;
  zGNGEdgePtrList eplist;
};

zListClass( zGNGUnitList, zGNGUnit, zGNGUnitData );

#define Z_GNG_DEFAULT_BATCH_TRIAL_SIZE         100
#define Z_GNG_DEFAULT_LIFETIME                  50
#define Z_GNG_DEFAULT_EXTENTION_COEFF_CLOSEST    0.2
#define Z_GNG_DEFAULT_EXTENTION_COEFF_NEIGHBOR   0.006
#define Z_GNG_DEFAULT_REDUCTION_RATIO1           0.5
#define Z_GNG_DEFAULT_REDUCTION_RATIO2           0.995

/*! \brief Growing Neural Gas class.
 *
 * The original paper about the algorithm is:
 * Bernd Fritzke, A Growing Neural Gas Learns Topology, in Proceedings of the 7th
 * International Conference on Neural Information Processing Systems, pp. 625-632, 1994.
 */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zGNG ){
  zGNGUnitList unitlist;
  zGNGEdgeList edgelist;
  int batch_trial_size;
  int lifetime;
  double extention_coeff_closest;
  double extention_coeff_neighbor;
  double error_reduction_ratio1;
  double error_reduction_ratio2;
  zVec (* sampler)(zVec,void*);
  zVec _v;
  zVec _vm;
};

#define zGNGSetBatchTrialSize(gng,s)         ( (gng)->batch_trial_size = (s) )
#define zGNGSetLifetime(gng,lt)              ( (gng)->lifetime = (lt) )
#define zGNGSetExtentionCoeffClosest(gng,e)  ( (gng)->extention_coeff_closest = (e) )
#define zGNGSetExtentionCoeffNeighbor(gng,e) ( (gng)->extention_coeff_neighbor = (e) )
#define zGNGSetErrorReductionRatio1(gng,r)   ( (gng)->error_reduction_ratio1 = (r) )
#define zGNGSetErrorReductionRatio2(gng,r)   ( (gng)->error_reduction_ratio2 = (r) )
#define zGNGSetSampler(gng,s)                ( (gng)->sampler = (s) )

/*! \brief initialize Growing Neural Gas. */
zGNG *zGNGInit(zGNG *gng, int dim, zVec (* f)(zVec,void*), void *util);

/*! \brief destroy Growing Neural Gas. */
bool zGNGDestroy(zGNG *gng);

/*! \brief update Growing Neural Gas. */
bool zGNGUpdate(zGNG *gng, void *util);

/* for debug */

void zGNGFPrint(FILE *fp, zGNG *gng);
void zGNGFWrite(FILE *fp, zGNG *gng);

/*! \} */

__END_DECLS

#endif /* __ZM_GNG_H__ */
