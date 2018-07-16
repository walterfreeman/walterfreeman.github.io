import numpy as np

dt=0.02
g=9.8
L=3
N=20
t=0

thetas=2*np.logspace(0,-2,N)
omegas=np.zeros(N)
periods=np.zeros(N)

while 1:
   # Aspel-Euler-Cromer update
   t=t+dt
   thetas=thetas+omegas*dt
   omegas=omegas-g/L*np.sin(thetas)*dt

   for i in range(N):
       # fancy colors, using some trig to pick vectors along the color wheel
       print ("C",0.4*np.sin(i)+0.6,0.4*np.sin(i+2)+0.6,0.4*np.sin(i+4)+0.6)
       x=L*np.sin(thetas[i])
       y=-L*np.cos(thetas[i])
       z=-i*0.2
       print ("l3",0,0,z,x,y,z)
       print ("c3",x,y,z,0.04)
   print ("F\n")
