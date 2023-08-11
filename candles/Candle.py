import numpy as np
import matplotlib.pyplot as plt
import math

u_1=[1.]
u_2=[2.]
v_1=[1.]
v_2=[1.]

time=[0.]
time_step=0.001
end_time=10.

eps=1.
au=1.
av=1.
c=1.
sigma0=1.
u0=1.

while time[-1] < end_time:
    du1= 1/eps*(-u_1[-1]+au*v_1[-1]*np.exp(u_1[-1]/(1+u_1[-1]/c)))\
         -sigma0*(1+u_1[-1]/c)**4+sigma0*u0*(1+u_2[-1]/c)**4
    du2 = 1/eps*(-u_2[-1]+au*v_2[-1]*np.exp(u_2[-1]/(1+u_2[-1]/c)))\
         -sigma0*(1+u_2[-1]/c)**4+sigma0*u0*(1+u_1[-1]/c)**4
    dv1 = 1-v_1[-1]-av*v_1[-1]*np.exp(u_1[-1]/(1+u_1[-1]/c))
    dv2 = 1-v_2[-1]-av*v_2[-1]*np.exp(u_2[-1]/(1+u_2[-1]/c))
    u_1.append(u_1[-1]+du1*time_step)
    u_2.append(u_2[-1]+du2*time_step)
    v_1.append(v_1[-1]+dv1*time_step)
    v_2.append(v_2[-1]+dv2*time_step)

    time.append(time[-1]+time_step)

time=np.array(time)
u_1=np.array(u_1)
u_2=np.array(u_2)
v_1=np.array(v_1)
v_2=np.array(v_2)

time_mask= np.logical_and(time>1,time<100)


fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('time')
ax1.set_ylabel('u', color=color)
ax1.plot(time, u_1, color=color)
ax1.plot(time, u_2, color=color, linestyle='dotted')
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('v', color=color)  # we already handled the x-label with ax1
ax1.plot(time, v_1, color=color)
ax1.plot(time, v_2, color=color, linestyle='dotted')
ax2.tick_params(axis='y', labelcolor=color)

#plt.xlim([2,10])
ax1.set_ylim([0,3])
ax2.set_ylim([0,2])

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
