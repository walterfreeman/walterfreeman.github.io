#include <stdio.h>
#include <math.h>

int main(void)
{
  double x,y,z;
  double GM=4*M_PI*M_PI; // this is the value for GM_sun when time is measured in years and distance in AU
  double vx,vy,vz;
  double dt=1e-5;
  double t=0;
  double r;
  int step=0, drawinterval=1000; // only draw a frame every 1000 timesteps
  x=1;
  y=0;
  z=0;
  vx=0;
  vy=2*M_PI * 0.8;
  vz=0;

  while (1)
  {
    // leapfrog update
    step++;
    t+=dt;
    x+=vx*dt/2;
    y+=vy*dt/2;
    z+=vz*dt/2;
    
    r=sqrt(x*x + y*y + z*z);
    vx-=GM*x/(r*r*r)*dt;
    vy-=GM*y/(r*r*r)*dt;
    vz-=GM*z/(r*r*r)*dt;
    
    x+=vx*dt/2;
    y+=vy*dt/2;
    z+=vz*dt/2;
    if (step % drawinterval == 0) // only animate one every drawinterval steps
    {
      printf("C 1 1 0\n"); // set color to yellow
      printf("c3 0 0 0 0.1\n"); // draw a big ball at the origin for the Sun
      printf("C 0.5 0.5 1\n"); // set color to blue
      printf("ct3 0 %e %e %e\n",x,y,z); // draw a smaller ball with a trail for the planet
      printf("F\n");
    }
  }
}
