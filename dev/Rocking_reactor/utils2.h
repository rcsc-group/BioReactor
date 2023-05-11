/**
# Various utility functions

## Default parameters and variables

The default maximum timestep and CFL number. */

/**
The *statsf()* function returns the minimum, maximum, volume sum,
standard deviation and volume for field *f*. */

typedef struct {
  double min, max, sum, sum2, stddev, volume;
} stats2;

stats2 statsf2 (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
  foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
	  reduction(max:max) reduction(min:min)) 
    if (dv() > 0. && f[] != nodata) {
      volume += dv();
      sum    += dv()*f[];
      sum2   += dv()*sq(f[]);
      if (f[] > max) max = f[];
      if (f[] < min) min = f[];
    }
  stats2 s;
  s.min = min, s.max = max, s.sum = sum, s.sum2 = sum2, s.volume = volume;
  if (volume > 0.)
    sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}
