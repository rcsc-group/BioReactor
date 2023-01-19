#include "embed.h"
#include "navier-stokes/centered.h"
#define mu(f) ( 1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2)  )
//#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/conserving.h"
#include "henry_oxy2.h"
#include "view2.h"
#include "tag.h"

#define EMBED            1
#define CONTACT          0
#define OXYGEN           1
#define OXYGEN_CIRCLE    0
#define OXYGEN_AIR       1
#define TRACER           1
#define ACCELERATION     1
#define AMR              1
#define REMOVE_DROP      1
#define CFL_COND         0 
#define DUMP             0
#define FIGURES          0
#define VIDEOS_original  0
#define VIDEOS_new       1
#define OUT_FILES        1
#define OUT_FILES_REFINE 0

//&&& Simulation setup &&&//
const double NN      = 512;     // resolution
const double t_change= 12.1487; // change in rocking motion
const double th_cont = 90;      // contact angle
const double t_mix   = 25;      // time for releasing tracers
const double t_dump  = 25;      // save the dump file
const double t_end_file = 100; // t<t_end_file; saving freq=dt_file
const double t_end   = 100;    // final time
const double dt_file = 0.6074/2;   // saving frq before t_end_file
const double dt_file_refine = 0.0281; // saving frq after t_end_file
const double dt_video= 0.6074/10;  // 0.0281
const double dt_oxy  = 0.001;
const int    i_fig   = 500;
const double CFL_num = 0.1;    // CFL number
const double N_output= 128;    // output resolution
const double remove_minsize   = 20;
const double remove_threshold = 1.0e-4;


//&&& Reactor geometry &&&//
double LL = 1.0;        // width
double Ly = 0.286;      // 0.50; height
double L_piv = 0.143;   // 0.25; distance to pivot point

//&&&&& Material properties &&&&&//
const double rho_w = 1.0e3;   // kg/m^3
const double rho_a = 1.225; // 0.01e3;   // kg/m^3
const double mu_w  = 1.0e-3;  // Pa-s @ 20C
const double mu_a  = 1.81e-5; // Pa-s - why multiplied 10?
const double grav  = 9.8;     // m/s^2
const double sigma = 0.0728;  // N/m
const double D_tracer_1 = 0.44e-9; // 0.44e-9; concentration = 1; diffusion in water
const double D_tracer_2 = 1.0e-30; // concentration = 0; diffusion in air
const double D_oxy_1 = 1.90e-9;    // 1.98e-5; m^2/s; phase = 1; diffusion in water
const double D_oxy_2 = 1.98e-5;    // 1.90e-9; m^2/s; phase = 0; diffusion in air
const double c_tracer_alpha = 1.0e30;    // concentration 1 to concentration 0
const double c_oxy_alpha    = 1./30;     // c_l=alpha*c_g; c1=alpha*c2; solubility; concentration 1 to concentration 0; air to water = 1/30;

//&&&&& Contact angle &&&&&//
#if CONTACT
vector h[];
h.t[left]  = contact_angle(th_cont*pi/180);
h.t[right] = contact_angle(th_cont*pi/180);
#endif

double Th_max, T_per, R_tr, x_tr, y_tr;
scalar ft[], c[], oxy[];   // for tracer and oxygen transfer
scalar * tracers = {ft}, * stracers = {c,oxy};
char buf1[50], buf2[50], buf3[20], buf4[20], buf5[20], buf6[20];
double (* gradient) (double, double, double) = minmod2;
double U0, Ub, Re_w, Re_a, We_w, Fr, rhor, mur, Pe_tracer_1, Pe_tracer_2, Pe_oxy_1, Pe_oxy_2, Th, Th_d, Th_2d, U_bio, w_bio, w_bio_st, T_per_st, T_bio, Th_max2, D_in_non, U_in_non;
int MINLEVEL, MAXLEVEL;
FILE * fp_stats;

int main(int argc, char * argv[]){

  // Operating condition
  double L_bio = atof(argv[1]);   // m - reference length scale
  double ANGLE = atof(argv[2]);   // degree
  double RPM   = atof(argv[3]);   // RPM

  L0 = LL;                 // domian size
  origin (-L0/2., -L0/2.); // origin
#if AMR == 0
  init_grid(NN);           // uniform grid
#endif

  // Mesh refinement
#if AMR
  MINLEVEL = atof(argv[4]);
  MAXLEVEL = atof(argv[5]);
  double F_MAX = 1e-6;
  double U_MAX = 0;
#endif

  //&&& Tracer conditions &&&//
  R_tr  = 0.0084/L_bio;    // tracer radius
  x_tr  = 0.;              // tracer position
  y_tr  = -Ly*0.5*0.5;

  Th_max = ANGLE*pi/180;   // radian
  T_per  = 60./RPM;        // sec
  
  double H_bio,V_bio;
  H_bio  = L_bio*Ly;
  V_bio  = L_bio/4*(H_bio + 0.5*L_bio*tan(Th_max));
  U_bio  = V_bio/(H_bio*0.5)/T_per; // Characteristic velocity
  T_bio  = L_bio/U_bio;             // Characteristic time
  w_bio  = 2*pi/T_per;          // angular velocity rad/s
  w_bio_st = w_bio*T_bio;       // dimensionless w_bio = 2*pi
  T_per_st = T_per/T_bio;       // dimensionless rocking period
  U0     = w_bio_st*Th_max;     // initial velocity (clockwise)
  Ub     = 0.;                  // Boundary velocity

  //dt_file = T_per_st/2;        // saving frq before t_end_file
  //dt_file_refine = T_per_st/5; // saving frq after t_end_file
  //dt_video= T_per_st/10;       // 0.0281

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
  
  //  Contact angle
#if CONTACT
  f.height = h;
#endif

#if TRACER  
  // D1 for c=1; D2 for c=0; c1=alpha*c2
  c.D1 = 1./Pe_tracer_1;
  c.D2 = 1./Pe_tracer_2;
  
  //c.alpha = 1.0e30;  // c1 = alpha*c2
  c.alpha = c_tracer_alpha;
#endif

#if OXYGEN
  oxy.D1 = 1./Pe_oxy_1;
  oxy.D2 = 1./Pe_oxy_2;
  oxy.alpha = c_oxy_alpha;
#endif

  //&&&&& Boundary conditions &&&&&//
  u.n[left]  = dirichlet(0.);
  u.t[left]  = dirichlet(0.);
  u.n[right] = dirichlet(0.);
  u.t[right] = dirichlet(0.);
  u.n[top] = dirichlet(0.);
  u.t[top] = dirichlet(0.);
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);

#if EMBED
  u.n[embed] = dirichlet(0.);
  u.t[embed] = dirichlet(0.);
#endif

  //***** Tracers: slope limiter *****//
  c.gradient = minmod2;
  oxy.gradient = minmod2;
  
  char name[200];
  sprintf(name, "logstats.dat");
  fp_stats = fopen(name, "w");
  
  // Iterations
  NITERMAX = 300; 
  TOLERANCE = 1.0e-4;
  
  run();
  
  fclose(fp_stats);
}

//***** Initial conditions  *****//
event init (t = 0)
{
  // Initial phases: y>0 -> volume fraction=0 (air); y<0 -> volume fraction=1
  fraction (f, -y);

  // solid
#if EMBED
  solid (cs,fs, intersection( -(y-0.5*Ly), -(-y-0.5*Ly) ) );
#endif
}

#if REMOVE_DROP
event remove_drop(i++){
  remove_droplets(f,remove_minsize,remove_threshold,false);  // remove droplets
  remove_droplets(f,remove_minsize,remove_threshold,true);   // remove bubbles
}
#endif

#if CFL_COND
event CFL_cond(i++){
  CFL = CFL_num;
}
#endif

#if TRACER
event tracer(t = t_mix){
  fraction(c, -(sq(x-x_tr) + sq(y-y_tr) - sq(R_tr)) );
}
#endif

#if OXYGEN
event oxygen (t=t_mix; i++){

#if OXYGEN_CIRCLE
  fraction(oxy, -(sq(x-0) + sq(y-Ly*0.5*0.5) - sq(0.084*Ly)) );
#endif

#if OXYGEN_AIR
  foreach(){
    //if (abs(f[]) < 1e-4)
    if ((f[] == 0) && (cs[]==1))
      oxy[] = 1.;
  }
#endif
}
#endif

#if ACCELERATION
//&&&&& acceleration &&&&&&//
event acceleration(i++)
{
  // initial stage
  if (t >= t_change){
    Th    = Th_max*sin(w_bio_st*t);
    Th_d  = w_bio_st*Th_max*cos(w_bio_st*t);
    Th_2d = -w_bio_st*w_bio_st*Th_max*sin(w_bio_st*t);
  }
  // regular state
  else if (t < t_change){
    Th_max2 = (Th_max-0)/(t_change-0)*t;
    Th    = Th_max2*sin(w_bio_st*t);
    Th_d  = w_bio_st*Th_max2*cos(w_bio_st*t);
    Th_2d = -w_bio_st*w_bio_st*Th_max2*sin(w_bio_st*t);
  }
   
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
}
#endif

#if AMR
event adapt( t=0 ){  
  // Refine the domain inside bioreactor
  refine(level<MAXLEVEL && (  (y > -0.7*Ly) && (y < 0.7*Ly) ));
}
#endif

#if DUMP
event dump(t=t_dump){
  dump(file="dump");
  snprintf(buf1, sizeof(buf1), "Dump_%d_%g_%d.txt",N,t,pid());
  FILE * out_all = fopen(buf1,"w");
  foreach(){
    fprintf(out_all,"%g %g %g %g %g %g \n",x,y,u.x[],u.y[],f[],c[]);
  }
  fclose(out_all);
}
#endif

event logstats (t+=0.1; t <= t_end) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    if (pid() == 0){
      fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
      fflush(fp_stats);
    }
}

/*
//event gfsview (t += 50.0) {
event gfsview (t += 0.1) {
    char name_gfs[200];
    sprintf(name_gfs,"Slices/Bioreactor-%0.1f.gfs",t);

    FILE* fp_gfs = fopen (name_gfs, "w");
    output_gfs(fp_gfs);
    fclose(fp_gfs);
}
*/
/*
event saveInterfaces (t += 10.0) {

    char nameInterfaces1[200];

    sprintf(nameInterfaces1,"Interfaces/interfacesLiquid-%0.1f.dat",t);

    FILE * fp1 = fopen(nameInterfaces1, "w");
    output_facets (f, fp1);	
    fclose(fp1);
}
*/

#if VIDEOS_original 
//***** Animations
event movies(t += dt_video; t<=t_end)
{
  scalar omega[], cc[], oxyy[];
  vorticity (u, omega);
  
  output_ppm(omega, file = "vorticity.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
	     min = -50, max = 50, linear = true, map = cool_warm);
  //output_ppm(u.x, file = "x-vel.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  //	     min = -0.1, max = 0.1, linear = true, map = cool_warm);
  //output_ppm(u.y, file = "y-vel.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  //  	     min = -0.1, max = 0.1, linear = true, map = cool_warm);  
  //output_ppm(f, file = "vol_fraction.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
//	     min = 0.0, max = 1.0, linear = true, map = cool_warm);
  //output_ppm(ft, file = "tracer.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  //	     linear = true, map = cool_warm);
  //output_ppm(cc, file = "tracer.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
  //	     min = 0, max = 1, linear = true, map = cool_warm);  
  //output_ppm(oxyy,file = "oxygen.mp4", box = { {-L0/2,-L0/2},{L0,L0} },
//	     min = 0, max = 0.04, linear = true, map = cool_warm);
}
#endif

#if VIDEOS_new
event movies_upgrade(t += dt_video; t<=t_end)
{
  scalar omega[], ff[], cc[], oxyy[];
  char timestring[100];

  vorticity (u,omega);

  // vorticity
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("omega",map=cool_warm,min=-50.0,max=50.0);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save("vorticity3.mp4");
  save("vorticity.png");

  // volume fraction
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("f",map=cool_warm,min=0.0,max=1.0);
  draw_vof("cs","fs");
  cells();
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save("volume_fraction3.mp4");
  save("volume_fraction.png");

  // tracer
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("c",map=cool_warm,min=0.0,max=0.1);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save("tracer3.mp4");
  save("tracer.png");

  // oxygen
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("oxy",map=cool_warm,min=0.0,max=0.04);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save("oxygen3.mp4");
  save("oxygen.png");
}
#endif

#if FIGURES
event figures(i += i_fig; t <= t_end)
{
  scalar omega[];
  vorticity (u,omega);
  
  output_ppm(c, file = "tracer.png");
  output_ppm(f, file = "vol_frac.png");
  //output_ppm(ft, file = "tracer.png");
  output_ppm(omega, file="vorticity.png");
  //output_ppm(u.x, file="x-vel.png");
  //output_ppm(u.y, file="y-vel.png");
  //output_ppm( cs, file="solid.png");
  output_ppm(oxy, file="oxygen.png");
}
#endif

#if OUT_FILES
event out_files(t+=dt_file; t<=t_end_file)
{
  ///*
  snprintf(buf1, sizeof(buf1), "Data_all/Data_all_%d_%.9g_%d.txt",N,t,pid());
  FILE * out_all = fopen(buf1,"wb");  
    
  foreach(){
    fprintf(out_all,"%g %g %g %g %g %g %g %g\n",x,y,u.x[],u.y[],f[],c[],cs[],oxy[]);
  }
  fclose(out_all);
  //*/
/*
  char name[80];
  sprintf(name,"Data_all/Output_%d_%.4g.dat",N,t);
  FILE * out_all2 = fopen(name,"w");
  output_field({u.x,u.y,f,c,cs,oxy},out_all2,n=N_output,linear=true);
  fclose(out_all2);
*/
}
#endif

#if OUT_FILES_REFINE
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
#endif
