# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 17:33:11 2018
@author: ddsch (Python 3.6)
David Schmidt

Description: Simple Ray Traceing Project

This code will propagate rays for field angles of 0 and ang. More rays then the other code which allows for more information on focal spot locations
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as animation

#Houcekeeping to set plot parameters and close any open figure
plt.close()
font = {'weight' : 'bold',
        'size'   : 22}

mpl.rc('font', **font)

#%% Lens Data (input as surfaces)
"""
Area to input lens system paramters. Two examples are given to show how to make the system. 
Assumes even number of surfaces.
"""
surf_pow=[0.015071,0,-0.011662,-0.016327,0.005736,0.01424]
surf_dis=[40,8.74,11.05,2.78,7.63,9.54,100]
surf_ind=[1,1.617,1,1.649,1,6,1]
surf_tot=6

lenstot=int(surf_tot/2)
lens_pos=[]
for ii in range(lenstot):
    lens_pos.append([sum(surf_dis[0:(ii*2+1)]),30,surf_dis[ii*2+1],surf_ind[ii*2+1]])

"""
surf_pow=[0.00125,0.00125]
surf_dis=[10,4,450]
surf_ind=[1,1.516824,1]
surf_tot=2

lenstot=1
lens_pos=[[surf_dis[0],40,surf_dis[1]]]
"""
#%% Ray data

#Ray data is populated below by giving amount of rays and an angle

nray_m=80
nray_a=80
angle=14
h_m=np.linspace(-20,20,nray_m)
h_a=np.linspace(-angle*np.pi/180*surf_dis[0]-20,-angle*np.pi/180*surf_dis[0]+20,nray_a)
ang_m=np.zeros(nray_m)
ang_a=np.ones(nray_a)*angle*np.pi/180

tot=len(h_m)+len(h_a)
h=np.reshape([h_m,h_a],(1,tot))[0]
ang=np.reshape([ang_m,ang_a],(1,tot))[0]


#%% paraxial equations

def refract(nin,nout,ang,y,lpow):
    return (nin*ang-y*lpow)/nout

def prop(y,ang,dis):
    return (ang*dis+y)

#%% generate data for lines
dx=1
tot_length=sum(surf_dis)
x=np.arange(0,tot_length,dx) 

rays=len(h)
ray_data=np.zeros((len(x),rays*2))
for ii in range(0,rays):
    ray_data[0,ii*2]=h[ii]
    ray_data[0,(ii*2)+1]=ang[ii]

sii=0
surfnext=surf_dis[sii]
for ii in range(len(x)-1):
    if x[ii]>=surfnext:
        sii=sii+1
        surfnext=surfnext+surf_dis[sii]
        for kk in range(0,rays):
            ray_data[ii,kk*2+1]=refract(surf_ind[sii-1],surf_ind[sii],ray_data[ii,kk*2+1],ray_data[ii,kk*2],surf_pow[sii-1])
    for kk in range(0,rays):
        ray_data[ii+1,kk*2]=prop(ray_data[ii,kk*2],ray_data[ii,kk*2+1],dx)
        ray_data[ii+1,kk*2+1]=ray_data[ii,kk*2+1]        

y_set=ray_data[:,0::2]
#%% Generate figure
xlimlow=0
xlimhigh=sum(surf_dis)+5
ylimlow=-(max(h)+25)
ylimhigh=-ylimlow

fig = plt.figure(figsize=(20, 10), dpi=80, facecolor='w', edgecolor='k')
ax = fig.add_subplot(111, aspect='equal', xlim=(xlimlow, xlimhigh), ylim=(ylimlow, ylimhigh),ylabel="Ray height (mm)",xlabel="Distance (mm)",title="Ray Tracing System")
ax.grid()
#%% Place paraxial lenses
for ii in range(0,lenstot):
    lens_use=lens_pos[ii]
    xl=lens_use[0]
    yl=lens_use[1]
    wl=lens_use[2]
    ax.add_patch(mpl.patches.Rectangle((xl,-yl), width=wl, height=2*yl, angle=0.0, color='b',alpha=lens_use[3]/max(surf_ind)))
    

#%% generate lines for animation from ray plot
plotcols=["blue"]*nray_m+["red"]*nray_a

lines = []
for index in range(rays):
    lobj = ax.plot([],[],lw=2,color=plotcols[index])[0]
    lines.append(lobj)

def init():
    for line in lines:
        line.set_data([],[])
    return lines

def animate(i):
    xu=[]
    yu=[]
    for ll in range(rays):
        xu.append(x[0:i])
        yu.append(y_set[0:i,ll])
    for lnum,line in enumerate(lines):
        line.set_data(xu[lnum], yu[lnum]) # set data for each line separately. 
    return lines

ani = animation.FuncAnimation(fig, animate, interval=1, blit=True, init_func=init)

plt.show()