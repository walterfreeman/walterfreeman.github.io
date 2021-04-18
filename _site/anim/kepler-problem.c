#include <stdio.h>
#include <math.h>



void zeromomentum(double& v1x, double& v1y, double& v2x, double& v2y, double m1, double m2)
{
  double v_center_x;
  double v_center_y;
  
  v_center_x = (v1x*m1 + v2x*m2) / (m1+m2);
  v_center_y = (v1y*m1 + v2y*m2) / (m1+m2);

  v1x -= v_center_x;
  v2x -= v_center_x;
  v1y -= v_center_y;
  v2y -= v_center_y;
}












int main(void)
{
  double x2,y2;
  double v2x,v2y;
  double x1,y1;
  double v1x,v1y;
  double G=4*M_PI*M_PI; // this is the value for GM_sun when time is measured in years and distance in AU
  double dt=1e-5;
  double t=0;
  double r;
  double m1=0.1, m2=1;
  int step=0, drawinterval=1000; // only draw a frame every 1000 timesteps

  x1=5;
  y1=0;
  v1x=0;
  v1y=2*M_PI * 0.3;
  x2=0;
  y2=0;
  v2x=v2y=0;
 
  zeromomentum(v1x, v1y, v2x, v2y, m1, m2);










































  while (1)
  {
    // leapfrog update
    step++;
    t+=dt;
    x1+=v1x*dt/2;
    y1+=v1y*dt/2;
    x2+=v2x*dt/2;
    y2+=v2y*dt/2;
    
    r=hypot(x1-x2, y1-y2);
    v1x-=G*m2*(x1-x2)/(r*r*r)*dt;
    v1y-=G*m2*(y1-y2)/(r*r*r)*dt;
    v2x+=G*m1*(x1-x2)/(r*r*r)*dt;
    v2y+=G*m1*(y1-y2)/(r*r*r)*dt;
    
    x1+=v1x*dt/2;
    y1+=v1y*dt/2;
    x2+=v2x*dt/2;
    y2+=v2y*dt/2;
  

     if (step % drawinterval == 0) // only animate one every drawinterval steps
    {
      printf("C 1 1 0\n"); // set color to yellow
      printf("ct3 0 %e %e 0 %e\n",x1,y1,sqrt(m1)*.1); // draw a smaller ball with a trail for the planet
      printf("C 0.5 0.5 1\n"); // set color to blue
      printf("ct3 1 %e %e 0 %e\n",x2,y2,sqrt(m2)*.1); // draw a smaller ball with a trail for the planet
      
      printf("T -0.9 -0.9\nTotal momentum: (% .3e, % .3e)\n",m1*v1x+m2*v2x,m2*v2y+m1*v1y);

      printf("F\n");
    }
  }
}
