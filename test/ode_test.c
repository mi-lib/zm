#include <zm/zm_ode.h>

/* ODE quadrature methods to be tested */
typedef struct{
  const char *name;
  zODE *(* create)(zODE*,int,int,zVec (*)(double,zVec,void*,zVec));
  void (* destroy)(zODE*);
  zVec (* update)(zODE*,double,zVec,double,void*);
  int step;  /* step number of history (only for multistep methods) */
  int order; /* theoretical order of accuracy */
} method_t;

#define METHOD(type,step,order) { #type, zODECreate##type, zODEDestroy##type, zODEUpdate##type, step, order }

static const method_t method_onestep[] = {
  METHOD( Euler,  0, 1 ),
  METHOD( Heun,   0, 3 ), /* named after the third-order Runge-Kutta method by Heun */
  METHOD( RK4,    0, 4 ),
  METHOD( RKG,    0, 4 ),
  METHOD( RKF45,  0, 5 ), /* the fifth-order solution is adopted */
  METHOD( CK45,   0, 5 ),
  METHOD( DP45,   0, 5 ),
  METHOD( BEuler, 0, 1 ),
  METHOD( Gauss,  0, 4 ),
  METHOD( Radau,  0, 3 ),
  METHOD( TR,     0, 2 ),
};
#define METHOD_ONESTEP_NUM ( sizeof(method_onestep)/sizeof(method_t) )

static void method_assign(zODE *ode, const method_t *method)
{
  ode->create = method->create;
  ode->destroy = method->destroy;
  ode->update = method->update;
  ode->order = method->order;
  zODEAssignFunc( ode, NULL, NULL );
}

static zODE *method_create(zODE *ode, const method_t *method, int dim, int step, zVec (* f)(double,zVec,void*,zVec))
{
  method_assign( ode, method );
  return zODECreate( ode, dim, step, f );
}

/* x' = -y, y' = x (unit circle) */
zVec circle(double t, zVec x, void *dummy, zVec dx)
{
  zVecSetElemNC( dx, 0,-zVecElemNC(x,1) );
  zVecSetElemNC( dx, 1, zVecElemNC(x,0) );
  return dx;
}

/* error of the solution of the unit circle from the exact one at the time t */
double circle_error(zVec x, double t, double r)
{
  return sqrt( zSqr( zVecElemNC(x,0) - r*cos(t) ) + zSqr( zVecElemNC(x,1) - r*sin(t) ) );
}

/* integrate the unit circle from t=0 to t=T with a one-step method, and return the error at T. */
double circle_integrate(const method_t *method, double dt, double T)
{
  zODE ode;
  zVec x;
  int i, n;
  double err;

  n = (int)( T / dt + 0.5 );
  method_create( &ode, method, 2, method->step, circle );
  x = zVecCreateList( 2, 1.0, 0.0 );
  for( i=0; i<n; i++ )
    zODEUpdate( &ode, i*dt, x, dt, NULL );
  err = circle_error( x, n*dt, 1.0 );
  zVecFree( x );
  zODEDestroy( &ode );
  return err;
}

/* observed order of accuracy from errors with time steps dt and dt/2 */
double observed_order(double err1, double err2)
{
  return log( err1 / err2 ) / log( 2 );
}

#define ORDER_TOL 0.2

/* 1. order of accuracy of one-step methods */
void assert_ode_order_onestep(void)
{
  const double dt = 0.05, T = 2.0;
  int i;
  double order;
  bool result = true;

  for( i=0; i<METHOD_ONESTEP_NUM; i++ ){
    order = observed_order( circle_integrate( &method_onestep[i], dt, T ), circle_integrate( &method_onestep[i], 0.5*dt, T ) );
    if( fabs( order - method_onestep[i].order ) > ORDER_TOL ){
      eprintf( "%s: observed order = %g (expected %d)\n", method_onestep[i].name, order, method_onestep[i].order );
      result = false;
    }
  }
  zAssert( zODEUpdate (order of accuracy of one-step methods), result );
}

/* x' = t^6 for t>=0 and x' = 0 for t<0, whose solution is x = t^7/7 for t>=0. Since x is constant and
 * x' is zero for t<0, the initial histories of Adams's and Gear's methods are consistent with it. */
zVec poly6(double t, zVec x, void *dummy, zVec dx)
{
  zVecSetElemNC( dx, 0, t > 0 ? pow( t, 6 ) : 0 );
  return dx;
}

/* integrate poly6 from t=0 to t=1 with a multistep method, and return the error at t=1. */
double poly6_integrate(const method_t *method, int step, double dt)
{
  zODE ode;
  zVec x;
  int i, n;
  double err;

  n = (int)( 1.0 / dt + 0.5 );
  method_create( &ode, method, 1, step, poly6 );
  x = zVecCreateList( 1, 0.0 );
  if( method->create == zODECreateGear ) zODEInitHistoryGear( &ode, x );
  for( i=0; i<n; i++ )
    zODEUpdate( &ode, i*dt, x, dt, NULL );
  err = fabs( zVecElemNC(x,0) - 1.0/7 );
  zVecFree( x );
  zODEDestroy( &ode );
  return err;
}

/* 2. order of accuracy of multistep methods */
void assert_ode_order_multistep(void)
{
  const method_t adams = METHOD( Adams, 0, 0 ), gear = METHOD( Gear, 0, 0 );
  const double dt = 0.00625;
  int step;
  double order;
  bool result = true;

  for( step=1; step<=6; step++ ){ /* the order of Adams's method with k steps is k. */
    order = observed_order( poly6_integrate( &adams, step, dt ), poly6_integrate( &adams, step, 0.5*dt ) );
    if( fabs( order - step ) > ORDER_TOL ){
      eprintf( "Adams (%d steps): observed order = %g (expected %d)\n", step, order, step );
      result = false;
    }
  }
  /* the order of Gear's method with k steps is k. Since the nonlinear equation is solved with the
   * tolerance zTOL, the error saturates for more steps. */
  for( step=1; step<=4; step++ ){
    order = observed_order( poly6_integrate( &gear, step, dt ), poly6_integrate( &gear, step, 0.5*dt ) );
    if( fabs( order - step ) > ORDER_TOL ){
      eprintf( "Gear (%d steps): observed order = %g (expected %d)\n", step, order, step );
      result = false;
    }
  }
  zAssert( zODEUpdate (order of accuracy of multistep methods), result );
}

/* integrate the unit circle with a large time step, and check the transition of the norm of the state.
 * sgn = 1: the norm monotonically increases, sgn = -1: monotonically decreases, sgn = 0: is preserved. */
bool circle_norm_check(const method_t *method, int sgn)
{
  const double dt = 0.5, T = 100.0, tol = 1.0e-9;
  zODE ode;
  zVec x;
  int i, n;
  double norm, norm_prev = 1.0;
  bool result = true;

  n = (int)( T / dt + 0.5 );
  method_create( &ode, method, 2, method->step, circle );
  x = zVecCreateList( 2, 1.0, 0.0 );
  for( i=0; i<n; i++ ){
    zODEUpdate( &ode, i*dt, x, dt, NULL );
    norm = zVecNorm( x );
    if( ( sgn == 0 && fabs( norm - 1 ) > tol ) ||
        ( sgn > 0 && norm <= norm_prev ) || ( sgn < 0 && norm >= norm_prev ) ){
      eprintf( "%s: norm = %.15g -> %.15g at step %d\n", method->name, norm_prev, norm, i );
      result = false;
      break;
    }
    norm_prev = norm;
  }
  zVecFree( x );
  zODEDestroy( &ode );
  return result;
}

/* 3. qualitative behaviors for a harmonic oscillation */
void assert_ode_circle_norm(void)
{
  const method_t euler = METHOD( Euler, 0, 1 ), beuler = METHOD( BEuler, 0, 1 ), radau = METHOD( Radau, 0, 3 );
  const method_t tr = METHOD( TR, 0, 2 ), gauss = METHOD( Gauss, 0, 4 );
  bool result = true;

  if( !circle_norm_check( &euler, 1 ) ) result = false;  /* explicit Euler method amplifies */
  if( !circle_norm_check( &beuler, -1 ) ) result = false; /* backward Euler method damps */
  if( !circle_norm_check( &radau, -1 ) ) result = false;  /* Radau method damps */
  if( !circle_norm_check( &tr, 0 ) ) result = false;     /* trapezoidal method preserves norm */
  if( !circle_norm_check( &gauss, 0 ) ) result = false;  /* Gauss method preserves norm */
  zAssert( zODEUpdate (norm of harmonic oscillation), result );
}

#define STIFF_A 100.0

/* x' = -A x */
zVec stiff(double t, zVec x, void *dummy, zVec dx)
{
  zVecSetElemNC( dx, 0, -STIFF_A * zVecElemNC(x,0) );
  return dx;
}

/* stability functions R(z) of A-stable one-step methods, with which x <- R(z) x for x' = -A x and z = -A dt */
double stiff_r_beuler(double z){ return 1 / ( 1 - z ); }
double stiff_r_tr(double z){ return ( 1 + z/2 ) / ( 1 - z/2 ); }
double stiff_r_gauss(double z){ return ( 1 + z/2 + z*z/12 ) / ( 1 - z/2 + z*z/12 ); }
double stiff_r_radau(double z){ return ( 1 + z/3 ) / ( 1 - 2*z/3 + z*z/6 ); }

/* 4. stability for a stiff equation */
void assert_ode_stiff(void)
{
  const method_t method[] = { /* A-stable one-step methods */
    METHOD( BEuler, 0, 1 ),
    METHOD( TR,     0, 2 ),
    METHOD( Gauss,  0, 4 ),
    METHOD( Radau,  0, 3 ),
  };
  double (* r[])(double) = { stiff_r_beuler, stiff_r_tr, stiff_r_gauss, stiff_r_radau };
  const method_t gear = METHOD( Gear, 3, 3 );
  const double dt = 0.1; /* A dt = 10, for which explicit methods diverge */
  const double tol = 1.0e-8;
  const int n = 20;
  zODE ode;
  zVec x;
  int i, k;
  double xabs, xabs_prev, x_exact;
  bool result = true;

  /* one-step methods: |x| monotonically decreases, and x coincides with R(z)^i, even after it becomes
   * much smaller than zTOL, since the nonlinear equation is solved with a tolerance relative to x. */
  for( k=0; k<sizeof(method)/sizeof(method_t); k++ ){
    method_create( &ode, &method[k], 1, method[k].step, stiff );
    x = zVecCreateList( 1, 1.0 );
    for( x_exact=xabs_prev=1.0, i=0; i<5*n; i++ ){
      zODEUpdate( &ode, i*dt, x, dt, NULL );
      x_exact *= r[k]( -STIFF_A * dt );
      if( ( xabs = fabs( zVecElemNC(x,0) ) ) >= xabs_prev || fabs( zVecElemNC(x,0) - x_exact ) > tol * fabs( x_exact ) ){
        eprintf( "%s: |x| = %g -> %g at step %d (exact: %g)\n", method[k].name, xabs_prev, xabs, i, fabs( x_exact ) );
        result = false;
        break;
      }
      xabs_prev = xabs;
    }
    zVecFree( x );
    zODEDestroy( &ode );
  }
  /* multistep method: |x| does not monotonically decrease, but is bounded and decays. */
  method_create( &ode, &gear, 1, gear.step, stiff );
  x = zVecCreateList( 1, 1.0 );
  zODEInitHistoryGear( &ode, x );
  for( i=0; i<n; i++ ){
    zODEUpdate( &ode, i*dt, x, dt, NULL );
    if( ( xabs = fabs( zVecElemNC(x,0) ) ) > 1.0 ){
      eprintf( "%s: |x| = %g at step %d\n", gear.name, xabs, i );
      result = false;
      break;
    }
  }
  if( fabs( zVecElemNC(x,0) ) > 1.0e-3 ){
    eprintf( "%s: |x| = %g after %d steps\n", gear.name, fabs( zVecElemNC(x,0) ), n );
    result = false;
  }
  zVecFree( x );
  zODEDestroy( &ode );
  zAssert( zODEUpdate (stability for a stiff equation), result );
}

#define FORCED_C 1.0e-13

/* x' = c - x */
zVec forced(double t, zVec x, void *c, zVec dx)
{
  zVecSetElemNC( dx, 0, *(double *)c - zVecElemNC(x,0) );
  return dx;
}

/* 5. implicit methods starting from the zero state */
void assert_ode_implicit_zero(void)
{
  const method_t method[] = {
    METHOD( BEuler, 0, 1 ),
    METHOD( TR,     0, 2 ),
    METHOD( Gauss,  0, 4 ),
    METHOD( Radau,  0, 3 ),
    METHOD( Gear,   1, 1 ), /* Gear's method with one step is identical with backward Euler method */
  };
  double (* r[])(double) = { stiff_r_beuler, stiff_r_tr, stiff_r_gauss, stiff_r_radau, stiff_r_beuler };
  const double cs[] = { 0, 1.0, FORCED_C }, dt = 0.1, tol = 1.0e-8;
  const int n = 10;
  zODE ode;
  zVec x;
  int i, j, k;
  double c, x_exact;
  bool result = true;

  /* since y = x - c satisfies y' = -y, the solution is x = c ( 1 - R(-dt)^i ). The tolerance of the residual
   * of the nonlinear equation is zero at the zero state, which has to be updated properly even if c dt is
   * much smaller than zTOL, and has to remain zero if c = 0. */
  for( j=0; j<sizeof(cs)/sizeof(double); j++ ){
    c = cs[j];
    for( k=0; k<sizeof(method)/sizeof(method_t); k++ ){
      method_create( &ode, &method[k], 1, method[k].step, forced );
      x = zVecCreateList( 1, 0.0 );
      if( method[k].create == zODECreateGear ) zODEInitHistoryGear( &ode, x );
      for( x_exact=1, i=0; i<n; i++ ){
        zODEUpdate( &ode, i*dt, x, dt, &c );
        x_exact *= r[k]( -dt );
        if( fabs( zVecElemNC(x,0) - c * ( 1 - x_exact ) ) > tol * fabs( c * ( 1 - x_exact ) ) ){
          eprintf( "%s: x = %.17g at step %d (exact: %.17g)\n", method[k].name, zVecElemNC(x,0), i, c * ( 1 - x_exact ) );
          result = false;
          break;
        }
      }
      zVecFree( x );
      zODEDestroy( &ode );
    }
  }
  zAssert( zODEUpdate (implicit methods from the zero state), result );
}

/* 6. rejection of the zero time step */
void assert_ode_zero_dt(void)
{
  zODE ode;
  zVec x;
  int k;
  bool result = true;

  for( k=0; k<METHOD_ONESTEP_NUM; k++ ){
    x = zVecCreateList( 2, 1.0, 0.0 );
    method_create( &ode, &method_onestep[k], 2, 0, circle );
    if( zODEUpdate( &ode, 0, x, 0, NULL ) ) result = false;
    zODEDestroy( &ode );
    method_assign( &ode, &method_onestep[k] );
    zODECreateDC( &ode, 2, 0, circle );
    if( zODEUpdateDC( &ode, 0, x, 0, NULL ) ) result = false;
    zODEDestroyDC( &ode );
    if( zVecElemNC(x,0) != 1.0 || zVecElemNC(x,1) != 0.0 ){
      eprintf( "%s: x = (%g, %g) after a zero time step\n", method_onestep[k].name, zVecElemNC(x,0), zVecElemNC(x,1) );
      result = false;
    }
    zVecFree( x );
  }
  zAssert( zODEUpdate (rejection of the zero time step), result );
}

/* x'' = -x */
zVec oscillator(double t, zVec x, zVec dx, void *dummy, zVec ddx)
{
  zVecSetElemNC( ddx, 0, -zVecElemNC(x,0) );
  return ddx;
}

/* energy error of a harmonic oscillation x'' = -x with x(0) = 1 and x'(0) = 0 */
double oscillator_energy_error(zVec x, zVec dx)
{
  return 0.5 * ( zSqr( zVecElemNC(x,0) ) + zSqr( zVecElemNC(dx,0) ) ) - 0.5;
}

/* integrate a harmonic oscillation for a long period, and return the maximum absolute energy errors in
 * the first and the second halves of the period. If drift is not null, it is set for true if the energy
 * monotonically decreases. */
void oscillator_integrate(zODE2 *ode, bool is_leapfrog, double *err1, double *err2, bool *drift)
{
  const double dt = 0.1;
  const int n = 10000;
  zVec x, dx;
  int i;
  double e, e_prev = 0;

  zODE2Create( ode, 1, 0, oscillator );
  x = zVecCreateList( 1, 1.0 );
  dx = zVecCreateList( 1, 0.0 );
  if( is_leapfrog ) zODE2InitHistoryLeapfrog( ode, x, dx, dt );
  *err1 = *err2 = 0;
  if( drift ) *drift = true;
  for( i=0; i<n; i++ ){
    zODE2Update( ode, i*dt, x, dx, dt, NULL );
    e = oscillator_energy_error( x, dx );
    if( i < n/2 ){
      if( fabs( e ) > *err1 ) *err1 = fabs( e );
    } else{
      if( fabs( e ) > *err2 ) *err2 = fabs( e );
    }
    if( drift && e >= e_prev ) *drift = false;
    e_prev = e;
  }
  zVecFree( x );
  zVecFree( dx );
  zODE2Destroy( ode );
}

/* 7. symplecticity of second-order ODE solvers */
void assert_ode2_symplectic(void)
{
  zODE2 ode;
  double err1, err2;
  bool drift, result = true;

  /* symplectic methods: the energy error is bounded. */
  zODE2Assign( &ode, Leapfrog, NULL, NULL, NULL, NULL );
  oscillator_integrate( &ode, true, &err1, &err2, NULL );
  if( err2 > 1.1 * err1 ){
    eprintf( "Leapfrog: energy error %g in the first half, %g in the second half\n", err1, err2 );
    result = false;
  }
  zODE2Assign( &ode, Symplectic, NULL, NULL, NULL, NULL );
  oscillator_integrate( &ode, false, &err1, &err2, NULL );
  if( err2 > 1.1 * err1 ){
    eprintf( "Symplectic: energy error %g in the first half, %g in the second half\n", err1, err2 );
    result = false;
  }
  /* a non-symplectic method: the energy monotonically drifts. */
  zODE2Assign( &ode, Regular, NULL, NULL, NULL, NULL );
  zODE2AssignRegular( &ode, RK4 );
  oscillator_integrate( &ode, false, &err1, &err2, &drift );
  if( !drift || err2 < 1.5 * err1 ){
    eprintf( "RK4: energy error %g in the first half, %g in the second half (%s)\n", err1, err2, drift ? "monotonic" : "not monotonic" );
    result = false;
  }
  zAssert( zODE2Update (symplecticity), result );
}

/* rejection of the zero time step by a second-order ODE solver */
bool ode2_zero_dt_check(zODE2 *ode, const char *name)
{
  zVec x, v;
  bool result = true;

  zODE2Create( ode, 1, 0, oscillator );
  x = zVecCreateList( 1, 1.0 );
  v = zVecCreateList( 1, 0.0 );
  if( zODE2Update( ode, 0, x, v, 0, NULL ) || zVecElemNC(x,0) != 1.0 || zVecElemNC(v,0) != 0.0 ){
    eprintf( "%s: x = %g, v = %g after a zero time step\n", name, zVecElemNC(x,0), zVecElemNC(v,0) );
    result = false;
  }
  if( zODE2Update( ode, 0, x, v, 0.1, NULL ) != x ){
    eprintf( "%s: zODE2Update does not return x\n", name );
    result = false;
  }
  zVecFree( x );
  zVecFree( v );
  zODE2Destroy( ode );
  return result;
}

/* 8. rejection of the zero time step by second-order ODE solvers */
void assert_ode2_zero_dt(void)
{
  zODE2 ode;
  bool result = true;

  zODE2Assign( &ode, Regular, NULL, NULL, NULL, NULL );
  zODE2AssignRegular( &ode, RK4 );
  if( !ode2_zero_dt_check( &ode, "Regular" ) ) result = false;
  zODE2Assign( &ode, Symplectic, NULL, NULL, NULL, NULL );
  if( !ode2_zero_dt_check( &ode, "Symplectic" ) ) result = false;
  zODE2Assign( &ode, Leapfrog, NULL, NULL, NULL, NULL );
  if( !ode2_zero_dt_check( &ode, "Leapfrog" ) ) result = false;
  zAssert( zODE2Update (rejection of the zero time step), result );
}

/* 9. accuracy of deferred correction */
void assert_ode_dc_accuracy(void)
{
  const double dt = 0.5, T = 10.0;
  zODE ode;
  zVec x;
  int i, k, n;
  double err, err_plain, tol;
  bool result = true;

  n = (int)( T / dt + 0.5 );
  /* the error in each time step is controlled to be in the order of ZM_ODE_TOL. */
  tol = 10 * n * ZM_ODE_TOL;
  /* multistep methods are not available for deferred correction, since their histories are
   * disturbed by trial updates. */
  for( k=0; k<METHOD_ONESTEP_NUM; k++ ){
    err_plain = circle_integrate( &method_onestep[k], dt, T );
    method_assign( &ode, &method_onestep[k] );
    zODECreateDC( &ode, 2, 0, circle );
    x = zVecCreateList( 2, 1.0, 0.0 );
    for( i=0; i<n; i++ )
      zODEUpdateDC( &ode, i*dt, x, dt, NULL );
    if( ( err = circle_error( x, n*dt, 1.0 ) ) > tol || err > err_plain ){
      eprintf( "%s: error = %g with deferred correction (tolerance %g, error without it %g)\n", method_onestep[k].name, err, tol, err_plain );
      result = false;
    }
    zVecFree( x );
    zODEDestroyDC( &ode );
  }
  zAssert( zODEUpdateDC (accuracy), result );
}

/* integrate the unit circle from t=0 to t=T with deferred correction, and return the error at T. */
double circle_integrate_dc(const method_t *method, double dt, double T)
{
  zODE ode;
  zVec x;
  int i, n;
  double err;

  n = (int)( T / dt + 0.5 );
  method_assign( &ode, method );
  zODECreateDC( &ode, 2, 0, circle );
  x = zVecCreateList( 2, 1.0, 0.0 );
  for( i=0; i<n; i++ )
    zODEUpdateDC( &ode, i*dt, x, dt, NULL );
  err = circle_error( x, n*dt, 1.0 );
  zVecFree( x );
  zODEDestroyDC( &ode );
  return err;
}

/* 10. order of accuracy of deferred correction */
void assert_ode_dc_order(void)
{
  /* Richardson's extrapolation raises the order of accuracy p by one, or by two for symmetric methods
   * (trapezoidal method). The time steps are chosen so that they are not divided by the step size control
   * and the errors are not dominated by the round-off error or the tolerance of the nonlinear equation
   * solver. */
  const struct{
    method_t method;
    double dt;
    int order;
  } test[] = {
    { METHOD( Euler,  0, 1 ), 0.001, 2 },
    { METHOD( BEuler, 0, 1 ), 0.001, 2 },
    { METHOD( Heun,   0, 3 ), 0.04,  4 },
    { METHOD( RK4,    0, 4 ), 0.04,  5 },
    { METHOD( RKG,    0, 4 ), 0.04,  5 },
    { METHOD( TR,     0, 2 ), 0.02,  4 },
    { METHOD( Radau,  0, 3 ), 0.04,  4 },
  };
  const double T = 2.0;
  int k;
  double order;
  bool result = true;

  for( k=0; k<sizeof(test)/sizeof(test[0]); k++ ){
    order = observed_order( circle_integrate_dc( &test[k].method, test[k].dt, T ), circle_integrate_dc( &test[k].method, 0.5*test[k].dt, T ) );
    if( fabs( order - test[k].order ) > ORDER_TOL ){
      eprintf( "%s: observed order = %g with deferred correction (expected %d)\n", test[k].method.name, order, test[k].order );
      result = false;
    }
  }
  zAssert( zODEUpdateDC (order of accuracy), result );
}

/* x' = -y, y' = x (circle) with counting the number of calls */
zVec circle_count(double t, zVec x, void *count, zVec dx)
{
  (*(int *)count)++;
  return circle( t, x, NULL, dx );
}

/* integrate a circle with the radius r from t=0 to t=1 with an adaptive step size control, and return the
 * relative error at t=1. The number of calls of the differential function is stored where count points. */
double circle_count_integrate(zODE *ode, bool is_dc, double r, int *count)
{
  const double dt = 0.1;
  const int n = 10;
  zVec x;
  int i;
  double err;

  if( is_dc ) zODECreateDC( ode, 2, 0, circle_count );
  else        zODECreate( ode, 2, 0, circle_count );
  x = zVecCreateList( 2, r, 0.0 );
  for( *count=0, i=0; i<n; i++ )
    is_dc ? zODEUpdateDC( ode, i*dt, x, dt, count ) : zODEUpdate( ode, i*dt, x, dt, count );
  err = circle_error( x, n*dt, r ) / r;
  zVecFree( x );
  if( is_dc ) zODEDestroyDC( ode ); else zODEDestroy( ode );
  return err;
}

/* 11. scale invariance of adaptive step size control */
void assert_ode_adaptive_scale(void)
{
  const double r[] = { 1.0e6, 1.0e12 };
  zODE ode;
  int i, count1, count;
  double err1, err;
  bool is_dc, result = true;

  /* the tolerance is relative for large states, so that the time step is divided in the same way as
   * for the unit circle. */
  for( is_dc=false; ; is_dc=true ){
    if( is_dc ) zODEAssign( &ode, RK4, NULL, NULL ); else zODEAssign( &ode, DP45, NULL, NULL );
    err1 = circle_count_integrate( &ode, is_dc, 1.0, &count1 );
    for( i=0; i<sizeof(r)/sizeof(double); i++ ){
      err = circle_count_integrate( &ode, is_dc, r[i], &count );
      if( count > count1 || err > 10 * err1 ){
        eprintf( "%s: radius %g, relative error = %g, %d calls (radius 1: %g, %d calls)\n", is_dc ? "RK4 with DC" : "DP45", r[i], err, count, err1, count1 );
        result = false;
      }
    }
    if( is_dc ) break;
  }
  zAssert( zODEUpdate (scale invariance of adaptive step size control), result );
}

/* integrate the unit circle from t=0 to t=n dt with an adaptive step size control, and return the error at
 * the end. The number of calls of the differential function is stored where count points. */
double circle_integrate_adaptive(const method_t *method, bool is_dc, double dt, int n, int *count)
{
  zODE ode;
  zVec x;
  int i;
  double err;

  method_assign( &ode, method );
  if( is_dc ) zODECreateDC( &ode, 2, 0, circle_count );
  else        zODECreate( &ode, 2, 0, circle_count );
  x = zVecCreateList( 2, 1.0, 0.0 );
  for( *count=0, i=0; i<n; i++ )
    is_dc ? zODEUpdateDC( &ode, i*dt, x, dt, count ) : zODEUpdate( &ode, i*dt, x, dt, count );
  err = circle_error( x, n*dt, 1.0 );
  zVecFree( x );
  if( is_dc ) zODEDestroyDC( &ode ); else zODEDestroy( &ode );
  return err;
}

/* 12. adaptive step size control for backward integration */
void assert_ode_adaptive_backward(void)
{
  const method_t erk[] = { METHOD( RKF45, 0, 5 ), METHOD( CK45, 0, 5 ), METHOD( DP45, 0, 5 ) };
  const double dt = 0.5;
  const int n = 20;
  int k, count_f, count_b;
  double err_f, err_b;
  bool result = true;

  /* since the circle is symmetric with respect to the time reversal, the time step has to be divided in
   * the same way for the backward integration as for the forward integration. */
  for( k=0; k<METHOD_ONESTEP_NUM+sizeof(erk)/sizeof(method_t); k++ ){
    const method_t *method = k < METHOD_ONESTEP_NUM ? &method_onestep[k] : &erk[k-METHOD_ONESTEP_NUM];
    bool is_dc = k < METHOD_ONESTEP_NUM;
    err_f = circle_integrate_adaptive( method, is_dc, dt, n, &count_f );
    err_b = circle_integrate_adaptive( method, is_dc, -dt, n, &count_b );
    if( count_b != count_f || fabs( err_b - err_f ) > 1.0e-6 * err_f ){
      eprintf( "%s%s: backward error = %g, %d calls (forward error = %g, %d calls)\n", method->name, is_dc ? " with DC" : "", err_b, count_b, err_f, count_f );
      result = false;
    }
  }
  zAssert( zODEUpdate (adaptive step size control for backward integration), result );
}

int main(void)
{
  assert_ode_order_onestep();
  assert_ode_order_multistep();
  assert_ode_circle_norm();
  assert_ode_stiff();
  assert_ode_implicit_zero();
  assert_ode_zero_dt();
  assert_ode2_symplectic();
  assert_ode2_zero_dt();
  assert_ode_dc_accuracy();
  assert_ode_dc_order();
  assert_ode_adaptive_scale();
  assert_ode_adaptive_backward();
  return EXIT_SUCCESS;
}
