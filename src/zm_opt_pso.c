/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_pso - optimization tools: Particle Swarm Optimization method.
 */
#include <zm/zm_opt.h>

/* particle of a swarm */
typedef struct{
  zVec x;
  zVec v;
  zVec mybest;
  zVec dmybest;
  zVec dbest;
  double inertia;
  double eval;
} zOptPSOParticle;

#define ZOPT_PSO_INERTIA_MIN 0.5
#define ZOPT_PSO_INERTIA_MAX 0.9

/* free memory for a particle */
void _zOptPSOParticleFree(zOptPSOParticle *p)
{
  zVecFreeAtOnce( 5, p->x, p->v, p->mybest, p->dmybest, p->dbest );
}

/* create a particle */
zOptPSOParticle *_zOptPSOParticleCreate(zOptPSOParticle *p, zVec min, zVec max, zVec vmin, zVec vmax, double (* f)(zVec,void*), void *util)
{
  p->x       = zVecAlloc( zVecSizeNC(min) );
  p->v       = zVecAlloc( zVecSizeNC(min) );
  p->mybest  = zVecAlloc( zVecSizeNC(min) );
  p->dmybest = zVecAlloc( zVecSizeNC(min) );
  p->dbest   = zVecAlloc( zVecSizeNC(min) );
  if( !p->x || !p->v || !p->mybest || !p->dmybest || !p->dbest ){
    _zOptPSOParticleFree( p );
    return NULL;
  }
  zVecInterDiv( min, max, zRandF(0,1), p->x );
  zVecInterDiv( vmin, vmax, zRandF(0,1), p->v );
  zVecCopyNC( p->x, p->mybest );
  p->inertia = zRandF( ZOPT_PSO_INERTIA_MIN, ZOPT_PSO_INERTIA_MAX );
  p->eval = f( p->x, util );
  return p;
}

/* update a particle */
bool _zOptPSOParticleUpdate(zOptPSOParticle *p, zVec min, zVec max, zVec vmin, zVec vmax, zVec best, double c1, double c2, double (* f)(zVec,void*), void *util)
{
  int i;
  double eval;

  zVecSubNC( p->mybest, p->x, p->dmybest );
  zVecSubNC( best, p->x, p->dbest );
  for( i=0; i<zVecSizeNC(p->x); i++ ){
    zVecElemNC(p->v,i) = zLimit( p->inertia * zVecElemNC(p->v,i)
                           + c1*zRandF(0,1) * zVecElemNC(p->dmybest,i)
                           + c2*zRandF(0,1) * zVecElemNC(p->dbest,i),
      zVecElemNC(vmin,i), zVecElemNC(vmax,i) );
    zVecElemNC(p->x,i) += zVecElemNC(p->v,i);
    /* cases to bounce */
    if( zVecElemNC(p->x,i) > zVecElemNC(max,i) ){
      zVecElemNC(p->x,i) = zVecElemNC(max,i);
      zVecElemNC(p->v,i) =-zVecElemNC(p->v,i);
    } else
    if( zVecElemNC(p->x,i) < zVecElemNC(min,i) ){
      zVecElemNC(p->x,i) = zVecElemNC(min,i);
      zVecElemNC(p->v,i) =-zVecElemNC(p->v,i);
    }
  }
  if( ( eval = f( p->x, util ) ) < p->eval ){ /* update my best */
    zVecCopyNC( p->x, p->mybest );
    p->eval = eval;
  }
  return zVecIsTiny( p->v );
}

/* swarm of particles */
typedef struct{
  int num;
  zOptPSOParticle *particle;
  double c1;
  double c2;
  zVec best;
  double eval;
} zOptPSOSwarm;

#define _zOptPSOSwarmSetC1(swarm,val) ( (swarm)->c1 = (val) )
#define _zOptPSOSwarmSetC2(swarm,val) ( (swarm)->c2 = (val) )

/* find the best evaluation value of a swarm of particles */
double _zOptPSOSwarmFindBest(zOptPSOSwarm *swarm)
{
  int i;

  swarm->best = swarm->particle[0].mybest;
  swarm->eval = swarm->particle[0].eval;
  for( i=1; i<swarm->num; i++ ){
    if( swarm->particle[i].eval < swarm->eval ){
      swarm->best = swarm->particle[i].mybest;
      swarm->eval = swarm->particle[i].eval;
    }
  }
  return swarm->eval;
}

/* destroy a swarm of particles */
void _zOptPSOSwarmDestroy(zOptPSOSwarm *swarm)
{
  int i;

  for( i=0; i<swarm->num; i++ )
    _zOptPSOParticleFree( &swarm->particle[i] );
  zFree( swarm->particle );
}

/* create a swarm of particles */
bool _zOptPSOSwarmCreate(zOptPSOSwarm *swarm, int num, zVec min, zVec max, zVec vmin, zVec vmax, double (* f)(zVec,void*), void *util, double c1, double c2)
{
  int i;

  if( !( swarm->particle = zAlloc( zOptPSOParticle, num ) ) ) return false;
  swarm->num = num;
  for( i=0; i<num; i++ ){
    if( !_zOptPSOParticleCreate( &swarm->particle[i], min, max, vmin, vmax, f, util ) ){
      _zOptPSOSwarmDestroy( swarm );
      return false;
    }
  }
  _zOptPSOSwarmSetC1( swarm, c1 );
  _zOptPSOSwarmSetC2( swarm, c2 );
  _zOptPSOSwarmFindBest( swarm );
  return true;
}

/* update a swarm of particles */
bool _zOptPSOSwarmUpdate(zOptPSOSwarm *swarm, zVec min, zVec max, zVec vmin, zVec vmax, double (* f)(zVec,void*), void *util)
{
  int i;
  bool is_converged = true;

  _zOptPSOSwarmFindBest( swarm );
  for( i=0; i<swarm->num; i++ )
    if( !_zOptPSOParticleUpdate( &swarm->particle[i], min, max, vmin, vmax, swarm->best, swarm->c1, swarm->c2, f, util ) )
      is_converged = false;
  return is_converged;
}

/* print out particles in a swarm */
#ifdef DEBUG
void zOptPSOSwarmFPrint(FILE *fp, zOptPSOSwarm *swarm)
{
  int i;

  for( i=0; i<swarm->num; i++ ){
    fprintf( fp, "%.10g ", swarm->particle[i].eval );
    zVecValueFPrint( fp, swarm->particle[i].x );
  }
}
#endif

/* solve an optimization problem by Particle Swarm Optimization method. */
int zOptSolvePSO(double (* f)(zVec,void*), void *util, zVec min, zVec max, int iter, double tol, int num, double c1, double c2, double vel_rate, zVec ans, double *eval)
{
  zOptPSOSwarm swarm;
  zVec vmin, vmax;
  int i = -1;

  if( !zVecSizeEqual( min, ans ) || !zVecSizeEqual( max, ans ) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return -1;
  }
  vmin = zVecAlloc( zVecSizeNC(ans) );
  vmax = zVecAlloc( zVecSizeNC(ans) );
  if( !vmin || !vmax ) goto TERMINATE;
  zVecSubNC( max, min, vmax );
  zVecMulNCDRC( vmax, vel_rate );
  zVecRevNC( vmax, vmin );
  if( !_zOptPSOSwarmCreate( &swarm, num, min, max, vmin, vmax, f, util, c1, c2 ) ) return -1;
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    if( _zOptPSOSwarmUpdate( &swarm, min, max, vmin, vmax, f, util ) ) break;
  }
  zVecCopy( swarm.best, ans );
  *eval = swarm.eval;
  _zOptPSOSwarmDestroy( &swarm );
 TERMINATE:
  zVecFree( vmin );
  zVecFree( vmax );
  return i;
}

/* solve an optimization problem by Particle Swarm Optimization method. */
int zOptSolvePSODefault(double (* f)(zVec,void*), void *util, zVec min, zVec max, int iter, double tol, zVec ans, double *eval)
{
  return zOptSolvePSO( f, util, min, max, iter, tol, ZOPT_PSO_DEFAULT_NUM, ZOPT_PSO_DEFAULT_C1, ZOPT_PSO_DEFAULT_C2, ZOPT_PSO_DEFAULT_VEL_RATE, ans, eval );
}
