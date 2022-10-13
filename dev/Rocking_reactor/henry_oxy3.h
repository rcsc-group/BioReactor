
//# Advection/diffusion of a soluble tracer

attribute {
  double D1, D2, alpha;  //D1 for f = 1, D2 for f = 0
  scalar phi1, phi2; // private
}

extern scalar * stracers;

/**
## Defaults

On trees we need to ensure conservation of the tracer when
refining/coarsening. */

#if TREE
event defaults (i = 0)
{
  for (scalar s in stracers) {
#if EMBED
    s.refine = s.prolongation = refine_embed_linear;
#else
    s.refine  = refine_linear;
#endif
    s.restriction = restriction_volume_average;
    s.dirty = true;
  }
}
#endif // TREE


// Advection

static scalar * phi_tracers = NULL;

event vof (i++)
{

  phi_tracers = f.tracers;

  ///*
  for (scalar c in stracers) {

    scalar phi1 = new scalar, phi2 = new scalar;
    //scalar aa = new scalar;
    
    c.phi1 = phi1, c.phi2 = phi2;

    // scalar_clone: src/grid/cartesian-common.h
    scalar_clone (phi1, c);
    scalar_clone (phi2, c);
    phi2.inverse = true;
    // phi1.inverse = true;
    
    f.tracers = list_append (f.tracers, phi1);
    f.tracers = list_append (f.tracers, phi2);

    foreach() {
      double a = c[]/(f[]*c.alpha + (1. - f[]));
      //aa[]   = c[]/(f[]*c.alpha + (1. - f[]));
      phi1[] = a*f[]*c.alpha;
      phi2[] = a*(1. - f[]);
      //fprintf(stdout,"%.3g, %.3g, %.3g %.3g %.3g\n",phi1[],phi2[],f[],c[],aa[]);
    }
    
    // added
    //vof_advection({phi1},i);
    //vof_advection({phi2},i);
  }
  //*/
}

// Diffusion
///*
struct HDiffusion {
  face vector D;     // alpha
  face vector beta;  // newly added for the second term of diffusion
#if EMBED
  double (* embed_flux) (Point, scalar, vector, double *);
#endif
};

static void h_relax (scalar * al, scalar * bl, int l, void * data)
{
  // similar to relax in poisson.h
  // alpha = D; no lambda
  scalar a = al[0], b = bl[0];
  struct HDiffusion * p = (struct HDiffusion *) data;
  face vector D = p->D, beta = p->beta;

  scalar c = a;
  foreach_level_or_leaf (l) {
    // -lambda[] = cm[]/dt
    double n = - sq(Delta)*b[], d = cm[]/dt*sq(Delta);
    foreach_dimension() {  
      n += D.x[1]*a[1] + D.x[]*a[-1] +
	Delta*(beta.x[1]*a[1] - beta.x[]*a[-1])/2.; // added terms in henry.h
      d += D.x[1] + D.x[] -
	Delta*(beta.x[1] - beta.x[])/2.; // added terms in henry.h
    }

    ///*      
#if EMBED
    if (p->embed_flux){
      double c;
      double e = embed_flux (point, a, D, &c);
      n -= c*sq(Delta);
      d += e*sq(Delta);
    }
    if (!d)
      c[] = b[] = 0.;
    else
#endif // EMBED
    //*/
    
    c[] = n/d;
  }
}
//*/

///*
static double h_residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct HDiffusion * p = (struct HDiffusion *) data;
  // similar to residual in poisson.h
  // alpha = D; lambda = beta;
  face vector D = p->D, beta = p->beta;
  double maxres = 0.;
#if TREE
  /// conservative coarse/fine discretisation (2nd order) //
  face vector g[];
  foreach_face()
    // second term is added: beta.x[]
    g.x[] = D.x[]*face_gradient_x (a, 0) + beta.x[]*face_value (a, 0);
  
  foreach (reduction(max:maxres)) {
    // -lambda[] = cm[]/dt;
    res[] = b[] + cm[]/dt*a[];
    foreach_dimension()
      res[] -= (g.x[1] - g.x[])/Delta;

    ///*
    //EMBED
#if EMBED
    if (p->embed_flux){
      double c, e = embed_flux (point, a, D, &c);
      res[] += c - e*a[];
    }
#endif // EMBED
    //*/
    
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#else // !TREE
  // "naive" discretisation (only 1st order on trees) //
  foreach (reduction(max:maxres)) {
    res[] = b[] + cm[]/dt*a[];
    foreach_dimension()
      res[] -= (D.x[1]*face_gradient_x (a, 1) -
		D.x[0]*face_gradient_x (a, 0) +
		beta.x[1]*face_value (a, 1) -
		beta.x[0]*face_value (a, 0))/Delta;
    ///*
    //EMBED
#if EMBED
    if (p->embed_flux) {
      double c, e = embed_flux (point, a, D, &c);
      res[] += c - e*a[];
    }
#endif // EMBED
    //*/
    
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#endif // !TREE    
  return maxres;
}
//*/

event tracer_diffusion (i++)
{
  free (f.tracers);
  f.tracers = phi_tracers;

  for (scalar c in stracers) {

    /**
    The advected concentration is computed from $\phi_1$ and $\phi_2$ as
    $$
    c = \phi_1 + \phi_2
    $$
    and these fields are then discarded. */
    
    scalar phi1 = c.phi1, phi2 = c.phi2, r[];

    foreach() {
      c[] = phi1[] + phi2[];
      r[] = - cm[]*c[]/dt;
    }
    delete ({phi1, phi2});
        
    /**
    The diffusion equation for $c$ is then solved using the multigrid
    solver and the residual and relaxation functions defined above. */
    ///*
    face vector D[], beta[];
    foreach_face() {
      double ff = (f[] + f[-1])/2.;
      D.x[] = fm.x[]*c.D1*c.D2/(c.D1*(1. - ff) + ff*c.D2);
      beta.x[] = - D.x[]*(c.alpha - 1.)/
	(ff*c.alpha + (1. - ff))*(f[] - f[-1])/Delta;
    }
  
    restriction ({D, beta, cm});
    struct HDiffusion q;
    q.D = D;
    q.beta = beta;

    // before multigrid solver
    //foreach_face()
      // face vector q: no errors; no NaN
      //fprintf(stdout,"%.3g, %.3g, %.3g, %.3g \n",D.x[],D.y[],beta.x[],beta.y[]);

      // cell-centered: tracer and residual
      //fprintf(stdout,"before %.3g %.3g\n",c[],r[]);      

    // from mgstats poisson
    ///*
    scalar aa = c;
#if EMBED
    if (!q.embed_flux && aa.boundary[embed] != symmetry)
      q.embed_flux = embed_flux;
#endif //EMBED
    //*/
    
    mg_solve ({c}, {r}, h_residual, h_relax, &q);
    //poisson ({c}, {r}, h_residual, h_relax, &q);

    // Poisson solver
    // return poisson (f,r,p.D,lambda);

    // after multigrid solver
    //foreach()
      // cell-centered: tracer and residual      
      //fprintf(stdout,"after %.3g %.3g\n",c[],r[]);
    
    //*/
  }
}

/**
## References

~~~bib
@article{haroun2010,
  title = {Volume of fluid method for interfacial reactive mass transfer: 
           application to stable liquid film},
  author = {Haroun, Y and Legendre, D and Raynal, L},
  journal = {Chemical Engineering Science},
  volume = {65},
  number = {10},
  pages = {2896--2909},
  year = {2010},
  doi = {10.1016/j.ces.2010.01.012}
}

@hal{farsoiya2021, hal-03227997}
~~~
*/
