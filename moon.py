import numpy as np
import matplotlib.pyplot as plt

r=[np.array([0,2])]
v=[np.array([1.2,0.0])]
time_step=0.01
earth_radius=1
earth_r=np.array([0,0])
gamma_earth=1.
moon_radius=0.1
moon_r=np.array([10,0])
gamma_moon=5



def get_acceleration(r1,r2,gamma):
    dr=r2-r1
    return gamma/(np.linalg.norm(dr)**2)*dr

def add_speed( r1,r2,gamma,time_step):
    return get_acceleration(r1,r2,gamma)*time_step

def next_point(r,v,time_step):
    return r+v*time_step

N=1000
theta=0.0001
for i in range(N):
    moon_r = np.dot(np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]]),moon_r)
    v_new=v[-1]+add_speed(r[-1],earth_r,gamma_earth,time_step)+add_speed(r[-1],moon_r,gamma_moon,time_step)
    v.append(v_new)
    r.append(next_point(r[-1],v[-1],time_step=time_step))

r=np.array(r)
v=np.array(v)
fig, ax = plt.subplots()
plt.plot(r[:,0],r[:,1],'r')
plt.plot(r[-1,0],r[-1,1],'b.')
earth = plt.Circle(earth_r, earth_radius, color='grey')
moon = plt.Circle(moon_r, moon_radius, color='grey')
ax.add_artist(earth)
ax.add_artist(moon)
limit=np.max([11,np.max(np.linalg.norm(r,axis=1))*1.1])
plt.xlim([-limit,limit])
plt.ylim([-limit,limit])
ax.set_aspect('equal')
plt.show()
t=np.linspace(0,1,N+1)
plt.plot(t,np.linalg.norm(v,axis=1))
plt.show()
