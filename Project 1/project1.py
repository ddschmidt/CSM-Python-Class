# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 17:33:11 2018
@author: ddsch (Python 3.6)
David Schmidt

Description: Simple Ray Traceing Project.

This code is meant for the optical engineer that wishes to trace a marginal ray and chief ray through the system
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as animation

#Housekeeping to set plot parameters and close any open figure
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
surf_ind=[1,1.617,1,1.649,1,1.617,1]
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
#%% Ray data [marginal ray, chief ray]
h=[14.93368,-14*np.pi/180*surf_dis[0]-6.9574]
ang=[0,14*np.pi/180]

#h=[20,-5*np.pi/180*surf_dis[0]]
#ang=[0,5*np.pi/180]


#%% paraxial equations

def refract(nin,nout,ang,y,lpow):
    return (nin*ang-y*lpow)/nout

def prop(y,ang,dis):
    return (ang*dis+y)

#%% generate data for lines
dx=.1
tot_length=sum(surf_dis)
x=np.arange(0,tot_length,dx) 
marg_info=np.zeros((len(x),2))
chief_info=np.zeros((len(x),2))

marg_info[0,0]=h[0]
marg_info[0,1]=ang[0]
chief_info[0,0]=h[1]
chief_info[0,1]=ang[1]

sii=0
surfnext=surf_dis[sii]
for ii in range(len(x)-1):
    if x[ii]>=surfnext:
        sii=sii+1
        surfnext=surfnext+surf_dis[sii]
        marg_info[ii,1]=refract(surf_ind[sii-1],surf_ind[sii],marg_info[ii,1],marg_info[ii,0],surf_pow[sii-1])
        chief_info[ii,1]=refract(surf_ind[sii-1],surf_ind[sii],chief_info[ii,1],chief_info[ii,0],surf_pow[sii-1])
    marg_info[ii+1,0]=prop(marg_info[ii,0],marg_info[ii,1],dx)
    chief_info[ii+1,0]=prop(chief_info[ii,0],chief_info[ii,1],dx)
    marg_info[ii+1,1]=marg_info[ii,1]
    chief_info[ii+1,1]=chief_info[ii,1]
        


y_marg=marg_info[:,0]
y_chief=chief_info[:,0]
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
plotlays, plotcols = [2], ["black","red"]
lines = []
for index in range(2):
    lobj = ax.plot([],[],lw=2,color=plotcols[index])[0]
    lines.append(lobj)

def init():
    for line in lines:
        line.set_data([],[])
    return lines

def animate(i):
    xu = [x[0:i],x[0:i]]
    yu = [y_marg[0:i],y_chief[0:i]]
    for lnum,line in enumerate(lines):
        line.set_data(xu[lnum], yu[lnum]) # set data for each line separately. 
    return lines

ani = animation.FuncAnimation(fig, animate, interval=1, blit=True, init_func=init,frames=len(x),repeat=False)

plt.show()