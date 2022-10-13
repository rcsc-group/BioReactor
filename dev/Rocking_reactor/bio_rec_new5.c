#include "embed.h"
#include "navier-stokes/centered.h"
#define mu(f) ( 1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2)  )
//#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/conserving.h"
//#include "tracer.h"
#include "utils.h"
//#include "henry.h"
//#include "/gpfs/scratch/mkim79/basilisk/henry_oxy2.h"
#include "henry_oxy3.h"
//#include "view.h"

const double NN      = 128;      // resolution
const double t_change= 11.2278;  // change in rocking motion
//const double th_cont = 90;     // contact angle
const double t_mix   = 100.;    // time for releasing tracers
const double t_dump  = 100.;    // save the dump file
const double t_end_file = 100.; // t<t_end_file; saving freq=dt_file
const double t_end   = 500.;    // final time
const double dt_file = 0.1123;   // saving frq before t_end_file
const double dt_file_refine = 0.0281; // saving frq after t_end_file
const double dt_video= 0.0281;

///*
//&&&&&& Bioreactor parameters &&&&&//
const double L_bio  = 0.4;     // m - reference length scale
const double Th_max = 0.122;   // rad - maximum angle of rotation - 7 degree
const double T_per  = 1.5;     // sec - 40rpm-1.5; 20rpm-3;
const double D_in   = 10.0e-3;  // m
const double U_in   = 0.0053;  // 0.3316; m/s
//*/

/*
// &&&&& Mesh refinement &&&&&//
int    MINLEVEL = 8;
int    MAXLEVEL = 10;
double F_MAX    = 1e-6;
double U_MAX    = 1e-2;
*/

double LL = 1.0;    // width
double Ly = 0.5;   // height
double L_piv = 0.25; // distance to pivot point
int    i_fig = 500;

scalar ft[], c[], oxy[];   // for tracer and oxygen transfer
scalar * tracers = {ft}, * stracers = {c,oxy};
char buf1[50], buf2[50], buf3[20], buf4[20], buf5[20], buf6[20];
double (* gradient) (double, double, double) = minmod2;
double U0, Ub, Re_w, Re_a, We_w, Fr, rhor, mur, Pe_tracer_1, Pe_tracer_2, Pe_oxy_1, Pe_oxy_2, Th, Th_d, Th_2d, U_bio, w_bio, w_bio_st, T_bio, Th_max2, D_in_non, U_in_non;

//&&&&& Material properties &&&&&//
const double rho_w = 1.0e3;   // kg/m^3
const double rho_a = 1.225; // 0.01e3;   // kg/m^3
const double mu_w  = 1.0e-3;  // Pa-s @ 20C
const double mu_a  = 1.81e-5; // Pa-s - why multiplied 10?
const double grav  = 9.8;     // m/s^2
const double sigma = 0.0728;  // N/m
const double D_tracer_1 = 0.44e-9; // concentration = 1; diffusion in water
const double D_tracer_2 = 1.0e-30; // concentration = 0; diffusion in air
const double D_oxy_1 = 1.98e-5; // m^2/s; concentration = 1; diffusion in air
const double D_oxy_2 = 1.90e-9; // m^2/s; concentration = 0; diffusion in water
const double c_tracer_alpha = 1.0e30; // concentration 1 to concentration 0
const double c_oxy_alpha    = 30;     // c_l=alpha*c_g; c1=alpha*c2; solubility; concentration 1 to concentration 0; air to water = 1/30;

//&&&&& Contact angle &&&&&//
/*
vector h[];
h.t[left]  = contact_angle(th_cont*pi/180);
h.t[right] = contact_angle(th_cont*pi/180);
*/

int main(){
  
  //***** Initial conditions *****// 
  L0 = LL;                 // domain size
  //mask (y > 0.9 ? top:none);
  
  origin (-L0/2., -L0/2.); // origin
  init_grid(NN);           // resolution
  
  double H_bio,V_bio;
  H_bio  = L_bio*Ly;
  V_bio  = L_bio/4*(H_bio + 0.5*L_bio*tan(Th_max));
  
  //U_bio  = L_bio*Th_max/T_per;  // reference velocity m/s
  U_bio  = V_bio/(H_bio*0.5)/T_per;
  T_bio  = L_bio/U_bio;
  w_bio  = 2*pi/T_per;          // angular velocity rad/s
  w_bio_st = w_bio*T_bio;       // dimensionless w_bio = 2*pi
  U0     = w_bio_st*Th_max;     // initial velocity (clockwise)
  Ub     = 0.;                  // Boundary velocity

  D_in_non = D_in/L_bio;        // oxygen support by two inlets
  U_in_non = U_in/U_bio;

  /*
  dt_file  = (T_per/T_bio)/5;   // saving frq before t_end_file
  dt_file_refine = (T_per/T_bio)/10; // saving frq after t_end_file
  dt_video = (T_per/T_bio)/20;
  */

  //***** Dimensionless number *****//
  Re_w = rho_w*U_bio*L_bio/mu_w;
  Re_a = rho_a*U_bio*L_bio/mu_a;
  We_w = rho_w*U_bio*U_bio*L_bio/sigma;
  Fr   = U_bio/sqrt(grav*L_bio);
  rhor = rho_a/rho_w;
  mur  = mu_a/mu_w;
  Pe_tracer_1 = U_bio*L_bio/D_tracer_1;
  Pe_tracer_2 = U_bio*L_bio/D_tracer_2;
  Pe_oxy_1    = U_bio*L_bio/D_oxy_1;
  Pe_oxy_2    = U_bio*L_bio/D_oxy_2;
    
  rho1 = 1.0;
  rho2 = rho1*rhor;
  mu1  = 1./Re_w;
  mu2  = mur*mu1;
  f.sigma = 1./We_w;

  // Contact angle
  //f.height = h;
  
  // D1 for c=1; D2 for c=0; c1=alpha*c2
  c.D1 = 1./Pe_tracer_1;
  c.D2 = 1./Pe_tracer_2;
  
  //c.alpha = 1.0e30;  // c1 = alpha*c2
  c.alpha = c_oxy_alpha;

  oxy.D1 = 1./Pe_oxy_1;
  oxy.D2 = 1./Pe_oxy_2;
  oxy.alpha = c_oxy_alpha;

  //&&&&& Boundary conditions &&&&&//
  u.n[left]  = dirichlet(0.);
  u.t[left]  = dirichlet(0.);
  u.n[right] = dirichlet(0.);
  u.t[right] = dirichlet(0.);
  u.n[top] = dirichlet(0.);
  u.t[top] = dirichlet(0.);
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);

  /*
  oxy[left]  = (y>=0.1) && (y<=0.15) && (f==1)		\
                ? dirichlet(1.) : dirichlet(0.);
  //oxy[right] = (y>=0.1) && (y<=0.15)			\
      //              ? dirichlet(1.) : dirichlet(0.);
  u.n[left]  = (y>=0.1) && (y<=0.15)	\
                ? dirichlet(1.) : dirichlet(0.);
  //u.n[right] = (y>=0.1) && (y<=0.15)			\
      //              ? dirichlet(-0.01) : dirichlet(0.);
  */

  u.n[embed] = dirichlet(0.);
  u.t[embed] = dirichlet(0.);

  //***** Tracers: slope limiter *****//
  //ft.gradient = minmod2;
  c.gradient = minmod2;
  oxy.gradient = minmod2;

  // Iterations
  NITERMAX = 500;
  TOLERANCE = 1.0e-5;
  
  CFL = 0.02;

  run();
}

//***** Initial conditions  *****//
event init (t = 0)
{
  // Initial phases: y>0 -> volume fraction=0 (air); y<0 -> volume fraction=1
  fraction (f, -y);

  // solid
  solid (cs,fs, intersection( -(y-0.5*Ly), -(-y-0.5*Ly) ) );

  // Tracer: positive -> ft=1; negative -> ft=0
  //fraction (ft, -(sq(x-0) + sq(y+0.25*LL) - sq(LL*0.1)) );

  // Oxygen: positive -> c=1; negative -> c=0
  // c.phi1 -> f(volume fraction); c.phi2 -> 1-f(1-volume fraction)
  // c.D1 -> volume fraction=1 -> water in bioreactor
  // c.D2 -> volume fraction=0 -> air in bioreactor
  // fraction (c, -(sq(x-0) + sq(y+0.25*LL) - sq(0.1*LL)) );  
  
  // initial velocity field
  // solid vol fraction = 1 (no solid) -> Ux0, solid vol fraction = 0 (w/solid) -> u.x = 0
  //foreach(){
  //u.x[] =  U0*(y+L_piv*cos(Th_max));
  // u.y[] = -U0*(x+L_piv*sin(Th_max));
    //u.x[] =  U0*y;
    //u.y[] = -U0*x;
  //}
}

event tracer(t = t_mix){
  //c.alpha = c_tracer_alpha;

  fraction(c, -(sq(x-0) + sq(y+Ly*0.5*0.5) - sq(0.084*Ly)) );
  /*
  foreach(){
    if (f[] < 1.)
      oxy[] = 1-f[];
  }
  */
}
///*
event oxygen(i++){
  if (t>=t_mix){
    foreach(){
      //if (f[] < 1.)
      if (f[] == 0.)
	oxy[] = 1-f[];
    }
  }
  /*
  if (t>=t_mix){
    u.n[left]  = (y>=0.20) && (y<=0.225) ? dirichlet(U_in_non) : dirichlet(0.);
    u.n[right] = (y>=0.20) && (y<=0.225) ? dirichlet(-U_in_non): dirichlet(0.);
    foreach(){
      //if ((y>=0.20) && (y<=0.225) && (x<=(-0.5*LL+3./NN)) ){
      if ((y>=0.20) && (y<=0.225) && ((x<=(-0.5*LL+3./NN)) || (x>=(0.5*LL-3./NN)) ) ){
	//u.x[] = 0.01;
	oxy[] = 1.;
      }
    }
  }
  */
  /*
  oxy[left]  = (y>=(Ly*0.5*0.5-0.5*D_in_non)) &&     \
               (y<=(Ly*0.5*0.5+0.5*D_in_non)) \
                ? dirichlet(1.) : dirichlet(0.);
  oxy[right] = (y>=(Ly*0.5*0.5-0.5*D_in_non)) &&  \
               (y<=(Ly*0.5*0.5+0.5*D_in_non)) \
                ? dirichlet(1.) : dirichlet(0.);
  u.n[left]  = (y>=(Ly*0.5*0.5-0.5*D_in_non)) &&        \
               (y<=(Ly*0.5*0.5+0.5*D_in_non))           \
                ? dirichlet(U_in_non) : dirichlet(0.);
  u.n[right] = (y>=(Ly*0.5*0.5-0.5*D_in_non)) && \
               (y<=(Ly*0.5*0.5+0.5*D_in_non)) \
                ? dirichlet(-U_in_non) : dirichlet(0.);
  */
}
//*/
//&&&&& acceleration &&&&&&//
event acceleration(i++)
{
  ///*
  // Rotation
  if (t >= t_change){
    Th    = Th_max*sin(w_bio_st*t);
    Th_d  = w_bio_st*Th_max*cos(w_bio_st*t);
    Th_2d = -w_bio_st*w_bio_st*Th_max*sin(w_bio_st*t);
  }
  else if (t < t_change){
    Th_max2 = (Th_max-0)/(t_change-0)*t;
    Th    = Th_max2*sin(w_bio_st*t);
    Th_d  = w_bio_st*Th_max2*cos(w_bio_st*t);
    Th_2d = -w_bio_st*w_bio_st*Th_max2*sin(w_bio_st*t);
  }
  //*/
  ///* 
  face vector av = a;
  // 1st: gravitational force, 2nd: Coriolis force
  // 3rd: centrifugal force,   4th: azimuthal force, 5th: no traslational force
  foreach_face(x)
    av.x[] = -sin(Th)/(Fr*Fr) + 2*Th_d*face_value(u.y,0)	\
    + Th_d*Th_d*(x+L_piv*sin(Th)) + Th_2d*(y+L_piv*cos(Th));
  foreach_face(y)
    av.y[] = -cos(Th)/(Fr*Fr) - 2*Th_d*face_value(u.x,0)	\
    + Th_d*Th_d*(y+L_piv*cos(Th)) - Th_2d*(x+L_piv*sin(Th));
  a = av;
  //*/
}

/*
//&&&&& Mesh refinement  &&&&&//
event adapt(i++){
  adapt_wavelet ({f,u}, (double []){F_MAX,U_MAX,U_MAX,U_MAX},maxlevel=MAXLEVEL,minlevel=MINLEVEL);
}
*/

event dump(t=t_dump){
  dump(file="dump");
  snprintf(buf1, sizeof(buf1), "Dump_%d_%g_%d.txt",N,t,pid());
  FILE * out_all = fopen(buf1,"w");
  foreach(){
    fprintf(out_all,"%g %g %g %g %g %g \n",x,y,u.x[],u.y[],f[],c[]);
  }
  fclose(out_all);
}

//***** Log files  *****//
event pro_logfile(i++){
  fprintf(stderr, "%d %d %g %d %d \n", i, N, t, mgp.i, mgu.i);    
}

//***** Animations
event movies(t += dt_video; t<=t_end)
{
  scalar omega[], m[];
  vorticity (u, omega);
  
  output_ppm(omega, file = "vorticity.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
	     min = -50, max = 50, linear = true, map = cool_warm);
  output_ppm(u.x, file = "x-vel.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
	     min = -0.1, max = 0.1, linear = true, map = cool_warm);
  output_ppm(u.y, file = "y-vel.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  	     min = -0.1, max = 0.1, linear = true, map = cool_warm);  
  output_ppm(f, file = "vol_fraction.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
	     min = -0.1, max = 1.1, linear = true, map = cool_warm);
  //output_ppm(ft, file = "tracer.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  //	     linear = true, map = cool_warm);
  output_ppm(c, file = "tracer.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  	     min = 0, max = 1, linear = true, map = cool_warm);  
  output_ppm(oxy,file = "oxygen.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
	     min = 0, max = 1, linear = true, map = cool_warm);
}

/*
event animationU (t += dt_video)
{
  view (tx = -0.4, ty = -0.4,  width = 500, height = 300, quat = {0, 0, 0.0, 0.0});
  clear();
  squares ("u.y", spread = -1, linear = true, map = cool_warm );   //draws colormap of field
  draw_vof ("f", lc = {0.0,0.0,0.0}, lw=1);                        //draw interface (line colour and linewidth?)
  box();
  cells();
  save ("HorizontalVelocity.mp4");
}
*/

event figures(i += i_fig; t <= t_end)
{
  scalar omega[];
  vorticity (u,omega);
  
  output_ppm(c, file = "tracer.png");
  output_ppm(f, file = "vol_frac.png");
  //output_ppm(ft, file = "tracer.png");
  output_ppm(omega, file="vorticity.png");
  output_ppm(u.x, file="x-vel.png");
  output_ppm(u.y, file="y-vel.png");
  //output_ppm( cs, file="solid.png");
  output_ppm(oxy, file="oxygen.png");
}
/*
event out_files(t+=dt_file; t<=t_end_file)
{
  snprintf(buf1, sizeof(buf1), "Data_all/Data_all_%d_%.9g_%d.txt",N,t,pid());
  //snprintf(buf2, sizeof(buf2), "alpha_%g.txt",t);
  FILE * out_all = fopen(buf1,"wb");
  //FILE * interf = fopen(buf2,"w");  
    
  foreach(){
    fprintf(out_all,"%g %g %g %g %g %g %g %g\n",x,y,u.x[],u.y[],f[],c[],cs[],oxy[]);
  }
  //output_facets(f,interf);
  
  fclose(out_all);  //fclose(interf);
}
*/
/*
event out_files_refine(t+=dt_file_refine; t<=t_end)
{
  if (t >= t_end_file){
    snprintf(buf1, sizeof(buf1), "Data_all/Data_all_%d_%.9g_%d.txt",N,t,pid());
    //snprintf(buf2, sizeof(buf2), "alpha_%g.txt",t);
    FILE * out_all = fopen(buf1,"wb");
    //FILE * interf = fopen(buf2,"w");  
    
    foreach(){
      fprintf(out_all,"%g %g %g %g %g %g %g %g\n",x,y,u.x[],u.y[],f[],c[],cs[],oxy[]);
    }
    //output_facets(f,interf);
  
    fclose(out_all);  //fclose(interf);
  }
}
*/
