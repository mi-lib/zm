/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_gng: Growing Neural Gas.
 */

#include <zm/zm_gng.h>

/* edge */

static void _zGNGEdgeDataFPrint(FILE *fp, zGNGEdgeData *edge)
{
  fprintf( fp, "(edge)\n" );
  fprintf( fp, " age = %d\n", edge->age );
  fprintf( fp, " unit1 = %p\n", edge->unit1 );
  fprintf( fp, " unit2 = %p\n", edge->unit2 );
}

/* unit */

static zGNGUnit *_zGNGUnitAlloc(zVec v)
{
  zGNGUnit *unit;

  if( !( unit = zAlloc( zGNGUnit, 1 ) ) ) return NULL;
  if( !( unit->data.v = zVecClone( v ) ) ){
    free( unit );
    return NULL;
  }
  unit->data.error = 0;
  zListInit( &unit->data.eplist );
  return unit;
}

static bool _zGNGUnitDestroy(zGNGUnit *unit)
{
  if( !zListIsEmpty( &unit->data.eplist ) ){
    ZRUNERROR( ZM_ERR_GNG_NONEMPTY_UNIT );
    return false;
  }
  zVecFree( unit->data.v );
  zFree( unit );
  return true;
}

static zGNGUnit *_zGNGUnitOther(zGNGUnit *unit, zGNGEdgeData *edge)
{
  return edge->unit1 == unit ? edge->unit2 : edge->unit1;
}

static bool _zGNGUnitPurgeEdge(zGNGUnit *unit, zGNGEdgeData *edge)
{
  zGNGEdgePtr *ep;

  zListForEach( &unit->data.eplist, ep ){
    if( ep->data == edge ){
      zListPurge( &unit->data.eplist, ep );
      zFree( ep );
      return true;
    }
  }
  return false; /* edge is not included */
}

static void _zGNGUnitFPrint(FILE *fp, zGNGUnit *unit)
{
  zGNGEdgePtr *ep;

  fprintf( fp, "(unit) %p\n", unit );
  fprintf( fp, " vec: " );
  zVecFPrint( fp, unit->data.v );
  fprintf( fp, " number of neighbors = %d\n", zListSize(&unit->data.eplist) );
  zListForEach( &unit->data.eplist, ep ){
    fprintf( fp, "  %p <- -> %p\n", ep->data->unit1, ep->data->unit2 );
  }
}

/* Growing Neural Gas */

static zVec _zGNGDefaultSampler(zVec v, void *util)
{
  zVecCopy( zVecListSelectRand( (zVecList*)util ), v );
  return v;
}

static void _zGNGPurgeUnit(zGNG *gng, zGNGUnit *unit)
{
  zListPurge( &gng->unitlist, unit );
  _zGNGUnitDestroy( unit );
}

static bool _zGNGPurgeEdge(zGNG *gng, zGNGEdge *ec)
{
  _zGNGUnitPurgeEdge( ec->data.unit1, &ec->data );
  if( zListIsEmpty( &ec->data.unit1->data.eplist ) ) _zGNGPurgeUnit( gng, ec->data.unit1 );
  _zGNGUnitPurgeEdge( ec->data.unit2, &ec->data );
  if( zListIsEmpty( &ec->data.unit2->data.eplist ) ) _zGNGPurgeUnit( gng, ec->data.unit2 );
  zListPurge( &gng->edgelist, ec );
  zFree( ec );
  return true;
}

static bool _zGNGConnect(zGNG *gng, zGNGUnit *unit1, zGNGUnit *unit2)
{
  zGNGEdge *ec;
  zGNGEdgePtr *ep1, *ep2;

  ec  = zAlloc( zGNGEdge, 1 );
  ep1 = zAlloc( zGNGEdgePtr, 1 );
  ep2 = zAlloc( zGNGEdgePtr, 1 );
  if( !ec || !ep1 || !ep2 ){
    ZALLOCERROR();
    zFree( ec );
    zFree( ep1 );
    zFree( ep2 );
    return false;
  }
  ep1->data = ep2->data = &ec->data;
  ec->data.age = 0;
  ec->data.unit1 = unit1;
  ec->data.unit2 = unit2;
  zListInsertHead( &gng->edgelist, ec );
  zListInsertHead( &unit1->data.eplist, ep1 );
  zListInsertHead( &unit2->data.eplist, ep2 );
  return true;
}

static bool _zGNGDisconnect(zGNG *gng, zGNGUnit *unit1, zGNGUnit *unit2)
{
  zGNGEdge *ec;

  zListForEach( &gng->edgelist, ec ){
    if( ( ec->data.unit1 == unit1 && ec->data.unit2 == unit2 ) ||
        ( ec->data.unit1 == unit2 && ec->data.unit2 == unit1 ) ){
      _zGNGPurgeEdge( gng, ec );
      return true;
    }
  }
  return false; /* not connected */
}

/* initialize Growing Neural Gas. */
zGNG *zGNGInit(zGNG *gng, int dim, zVec (* f)(zVec,void*), void *util)
{
  zGNGUnit *unit1, *unit2;

  zListInit( &gng->unitlist );
  zListInit( &gng->edgelist );
  /* default parameters */
  zGNGSetBatchTrialSize( gng, Z_GNG_DEFAULT_BATCH_TRIAL_SIZE );
  zGNGSetLifetime( gng, Z_GNG_DEFAULT_LIFETIME );
  zGNGSetExtentionCoeffClosest( gng, Z_GNG_DEFAULT_EXTENTION_COEFF_CLOSEST );
  zGNGSetExtentionCoeffNeighbor( gng, Z_GNG_DEFAULT_EXTENTION_COEFF_NEIGHBOR );
  zGNGSetErrorReductionRatio1( gng, Z_GNG_DEFAULT_REDUCTION_RATIO1 );
  zGNGSetErrorReductionRatio2( gng, Z_GNG_DEFAULT_REDUCTION_RATIO2 );
  zGNGSetSampler( gng, f ? f : _zGNGDefaultSampler );
  gng->_v  = zVecAlloc( dim );
  gng->_vm = zVecAlloc( dim );
  if( !gng->_v || !gng->_vm ){
    ZALLOCERROR();
    zVecFree( gng->_v );
    zVecFree( gng->_vm );
    return NULL;
  }
  /* first two units */
  gng->sampler( gng->_v, util );
  if( !( unit1 = _zGNGUnitAlloc( gng->_v ) ) ) goto ZGNG_INIT1_ERROR;
  do{
    gng->sampler( gng->_v, util );
  } while( zVecIsEqual( gng->_v, unit1->data.v, zTOL ) );
  if( !( unit2 = _zGNGUnitAlloc( gng->_v ) ) ) goto ZGNG_INIT2_ERROR;
  zListInsertHead( &gng->unitlist, unit1 );
  zListInsertHead( &gng->unitlist, unit2 );
  if( !_zGNGConnect( gng, unit1, unit2 ) ) goto ZGNG_INIT2_ERROR;
  return gng;

 ZGNG_INIT2_ERROR:
  _zGNGUnitDestroy( unit2 );
 ZGNG_INIT1_ERROR:
  _zGNGUnitDestroy( unit1 );
  zVecFree( gng->_v );
  zVecFree( gng->_vm );
  return NULL;
}

/* destroy Growing Neural Gas. */
bool zGNGDestroy(zGNG *gng)
{
  zGNGEdge *ec;

  while( !zListIsEmpty( &gng->edgelist ) ){
    ec = zListHead(&gng->edgelist);
    _zGNGPurgeEdge( gng, ec );
  }
  return true;
}

static void _zGNGFindProximity2(zGNG *gng, zVec v, zGNGUnit **unit1, zGNGUnit **unit2)
{
  zGNGUnit *uc;
  double d, dmin1 = HUGE_VAL, dmin2 = HUGE_VAL;

  *unit1 = *unit2 = zListTail(&gng->unitlist);
  zListForEach( &gng->unitlist, uc ){
    d = zVecSqrDist( v, uc->data.v );
    if( d < dmin1 ){
      *unit2 = *unit1;
      dmin2 = dmin1;
      *unit1 = uc;
      dmin1 = d;
    } else
    if( d < dmin2 ){
      *unit2 = uc;
      dmin2 = d;
    }
  }
}

static zGNGUnit *_zGNGFindPivot(zGNG *gng)
{
  zGNGUnit *uc, *pivot;

  pivot = zListTail(&gng->unitlist);
  zListForEach( &gng->unitlist, uc ){
    if( uc->data.error >= pivot->data.error ) pivot = uc;
  }
  return pivot;
}

static zGNGUnit *_zGNGFindSeed(zGNG *gng, zGNGUnit *pivot)
{
  zGNGUnit *uc, *seed;
  zGNGEdgePtr *ec;

  seed = _zGNGUnitOther( pivot, zListTail(&pivot->data.eplist)->data );
  zListForEach( &pivot->data.eplist, ec ){
    uc = _zGNGUnitOther( pivot, ec->data );
    if( uc->data.error >= seed->data.error ) seed = uc;
  }
  return seed;
}

static void _zGNGAdapt(zGNG *gng, void *util)
{
  zGNGUnit *u1, *u2, *un;
  zGNGEdgePtr *ec;
  bool u2_is_neighbor = false;

  _zGNGFindProximity2( gng, gng->sampler( gng->_v, util ), &u1, &u2 );
  u1->data.error += zVecSqrDist( u1->data.v, gng->_v );
  zVecInterDivDRC( u1->data.v, gng->_v, gng->extention_coeff_closest );
  zListForEach( &u1->data.eplist, ec ){
    un = _zGNGUnitOther( u1, ec->data );
    zVecInterDivDRC( un->data.v, gng->_v, gng->extention_coeff_neighbor );
    if( un == u2 ){
      ec->data->age = 0;
      u2_is_neighbor = true;
    } else
      ec->data->age++;
  }
  if( !u2_is_neighbor )
    _zGNGConnect( gng, u1, u2 );
}

static void _zGNGRemoveAgedEdge(zGNG *gng)
{
  zGNGEdge *ec, *ecp;

  zListForEach( &gng->edgelist, ec ){
    if( ec->data.age > gng->lifetime ){
      ecp = zListCellPrev(ec);
      _zGNGPurgeEdge( gng, ec );
      ec = ecp;
    }
  }
}

static bool _zGNGGrow(zGNG *gng)
{
  zGNGUnit *ur, *u1, *u2;

  u1 = _zGNGFindPivot( gng );
  u2 = _zGNGFindSeed( gng, u1 );
  zVecMid( u1->data.v, u2->data.v, gng->_v );
  if( !( ur = _zGNGUnitAlloc( gng->_v ) ) ) return false;
  u1->data.error *= gng->error_reduction_ratio1;
  u2->data.error *= gng->error_reduction_ratio1;
  ur->data.error = u1->data.error;
  zListInsertHead( &gng->unitlist, ur );
  _zGNGConnect( gng, u1, ur );
  _zGNGConnect( gng, u2, ur );
  _zGNGDisconnect( gng, u1, u2 );
  zListForEach( &gng->unitlist, ur ){
    ur->data.error *= gng->error_reduction_ratio2;
  }
  return true;
}

/* update Growing Neural Gas. */
bool zGNGUpdate(zGNG *gng, void *util)
{
  int i;

  for( i=0; i<gng->batch_trial_size; i++ ){
    _zGNGAdapt( gng, util );
    _zGNGRemoveAgedEdge( gng );
  }
  return _zGNGGrow( gng );
}

/* for debug */

/* print out data of Growing Neural Gas. */
void zGNGFPrint(FILE *fp, zGNG *gng)
{
  zGNGUnit *uc;
  zGNGEdge *ec;

  fprintf( fp, "Number of network cell = %d\n", zListSize(&gng->unitlist) );
  zListForEach( &gng->unitlist, uc ){
    _zGNGUnitFPrint( fp, uc );
  }
  fprintf( fp, "Number of edges = %d\n", zListSize(&gng->edgelist) );
  zListForEach( &gng->edgelist, ec ){
    _zGNGEdgeDataFPrint( fp, &ec->data );
  }
}

/* output edges of Growing Neural Gas. */
void zGNGFWrite(FILE *fp, zGNG *gng)
{
  zGNGEdge *ec;

  zListForEach( &gng->edgelist, ec ){
    zVecDataFPrint( fp, ec->data.unit1->data.v );
    zVecDataFPrint( fp, ec->data.unit2->data.v );
    fprintf( fp, "\n" );
  }
}
