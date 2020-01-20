import numpy as np

dt=0.03
g=9.8
L=3
theta=1
omega=0
frame=0
t=0

while True:
    # Euler-Cromer-Aspel update
    theta=theta+omega*dt
    omega=omega-g/L*np.sin(theta)*dt
    frame+=1;
    t += dt 
    # Animation
    x=L*np.sin(theta)
    y=-L*np.cos(theta)
    print ("l",0,0,x,y)
    print ("c",x,y,0.05)
    print ("F")
