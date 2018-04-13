import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from scipy.spatial.distance import pdist, squareform
import scipy.integrate as integrate
import matplotlib.animation as animation

#constants
radiusH=0.1#120*10e-12
massH=1#1.6737*10e-27
box_length=100*radiusH
#this one is actually important, how many particles we put in the box
numparticles=100
###########################################################################################
#### area and physics ###
#### x,y,z,vx,vy,vz,m ###
class physicsvolume:
    def __init__(self,
                cornors=[-box_length/2,box_length/2,-box_length/2,box_length/2,-box_length/2,box_length/2],
                init_state=[[0,0,0,0,0,0,1],[0,0,0,0,0,0,0,1]],
                radius=radiusH):

        #things internal to the box
        self.init_state=np.asarray(init_state,dtype=float)
        self.state=self.init_state.copy()
        self.radius=radiusH
        self.time_elapsed=0
        self.cornors=cornors
        #self.ForceAtTime=0
        self.kineticEnergy=0
        #self.Temperature=0
        #self.Volume=box_length**3
        #self.Area=(box_length**2)*6
        #self.measuredpressure=0
        #self.pressure=0
        #self.magVelocityinit
        #self.magVelocityfinal

    def step(self, dt):
        #move time
        #self.ForceAtTime=0
        self.time_elapsed +=dt

            #move particles
        self.state[:,:3]+=dt*self.state[:,3:6]

    #Simple Collision Detection thing
            #finds distance between particles
                #produces an NxN array of distances between particles
        D=squareform(pdist(self.state[:,:3]))
            #checks if distance is less than the radius
        #for i in range(50):
        #    for j in range(50):
        #        if D[i,j]<2*self.radius+5 and i!=j:
        #            print("hey",str(i),str(j),str(D[i,j]))
        #            print(2*self.radius)

        #this is weird but it creates an array of bool values that just say if this particle is hitting another particle or not
        ind1,ind2 = np.where(D<2*self.radius)
        unique = (ind1<ind2)
        ind1 = ind1[unique]
        ind2 = ind2[unique]

            #calculate collisions
        for i1,i2 in zip(ind1,ind2):
            m1=massH
            m2=massH
            # positions
            r1=self.state[i1,:3]
            r2=self.state[i2,:3]
            # velocities
            v1=self.state[i1,3:6]
            v2=self.state[i2,3:6]
            #relative location and velocities
            r_rel=r1-r2
            v_rel=v1-v2
            #momentum of com
            v_cm=(m1*v1+m2*v2)/(m1+m2)
            #collisions of spheres
            rr_rel = np.dot(r_rel,r_rel)
            vr_rel = np.dot(v_rel,r_rel)
            v_rel=2*r_rel*vr_rel / rr_rel - v_rel

            #assign new velocities
            self.state[i1,3:6]=(v_cm+v_rel*m2/(m1+m2))
            self.state[i2,3:6]=(v_cm-v_rel*m1/(m1+m2))
            #print((v_cm+v_rel*m2/(m1+m2)))
        #boundary crossing
        #this just makes sure the balls don't clip through the wall
        crossed_x1 = (self.state[:,0]<self.cornors[0] + self.radius)
        crossed_x2 = (self.state[:,0]>self.cornors[1] - self.radius)
        crossed_y1 = (self.state[:,1]<self.cornors[2] + self.radius)
        crossed_y2 = (self.state[:,1]>self.cornors[3] - self.radius)
        crossed_z1 = (self.state[:,2]<self.cornors[4] + self.radius)
        crossed_z2 = (self.state[:,2]>self.cornors[5] - self.radius)

#        #print(crossed_x1)
#        #### THIS ADDS A FORCE ON THE WALL ####
#        # force from a collition F=m*dv/dt            dv=2*v dt=step length
#        if len(self.state[crossed_x1 | crossed_x2])==1:
#            #print(abs(self.state[crossed_x1|crossed_x2,3]))
#            self.ForceAtTime=2*float(abs(self.state[crossed_x1|crossed_x2,3]))*abs(self.state[crossed_x1|crossed_x2,6])/dt+self.ForceAtTime
#
#        if len(self.state[crossed_y1 | crossed_y2])==1:
#            #print(abs(self.state[crossed_x1|crossed_x2,3]))
#            self.ForceAtTime=2*float(abs(self.state[crossed_y1|crossed_y2,4]))*abs(self.state[crossed_y1|crossed_y2,6])/dt+self.ForceAtTime
#        if len(self.state[crossed_z1 | crossed_z2])==1:
#            #print(abs(self.state[crossed_x1|crossed_x2,3]))
#            self.ForceAtTime=2*float(abs(self.state[crossed_z1|crossed_z2,5]))*abs(self.state[crossed_z1|crossed_z2,6])/dt+self.ForceAtTime

        #this means something is really broken
        #if self.ForceAtTime<0:
            #print("ahh")
        #calculates pressure from the FORCE
        #self.measuredpressure=self.ForceAtTime*numparticles/self.Area

        self.state[crossed_x1,0]=self.cornors[0] + self.radius
        self.state[crossed_x2,0]=self.cornors[1] - self.radius

        self.state[crossed_y1,1]=self.cornors[2]+self.radius
        self.state[crossed_y2,1]=self.cornors[3]-self.radius

        self.state[crossed_z1,2]=self.cornors[4]+self.radius
        self.state[crossed_z2,2]=self.cornors[5]-self.radius

        #makes the balls "bounce" back the way they came from
        self.state[crossed_x1 | crossed_x2,3]*=-1
        self.state[crossed_y1 | crossed_y2,4]*=-1
        self.state[crossed_z1 | crossed_z2,5]*=-1

################################################################################################

#this creates an initial state of numparticle number of balls with an initial position velocity and mass
np.random.seed(12)
red=-0.5+np.random.random((numparticles,7))
#makes sure they are in the box
red[:,:3]*=box_length-0.1
#sets all the masses equal
red[:,6]=massH
#red[:,2]=0
#red[:,5]=0
#print("hey")
#print(0.5*red[:,6]*(red[:,3]**2+red[:,4]**2+red[:,5]**2))
#red=[[0,0,0,0,0,0,0],[10,10,0,-1,-1,0,0]]
#says let there be a box
box = physicsvolume(init_state=red,radius=1)
dt=1./30

#%% Animation in 3D
################################################################################################
#plots for animation
fig = plt.figure()
ax= p3.Axes3D(fig)

#plots for the particles
particles, = ax.plot([],[],[],'bo',ms=1)

# Setting the axes properties
ax.set_xlim3d([-box_length/2, box_length/2])
ax.set_xlabel('X')

ax.set_ylim3d([-box_length/2, box_length/2])
ax.set_ylabel('Y')

ax.set_zlim3d([-box_length/2, box_length/2])
ax.set_zlabel('Z')

########################## animation stuff

def init():
    global box, rect
    particles.set_data([],[])
    particles.set_3d_properties([])
    
    return particles,

def animate(i):
    global box, rect, dt, ax, fig, avgpress
    box.step(dt)
    ms = int(fig.dpi*box.radius*fig.get_figwidth()/np.diff(ax.get_xbound())[0])

    particles.set_data(box.state[:,0],box.state[:,1])
    particles.set_3d_properties(box.state[:,2])
    particles.set_markersize(ms)
    
    return particles, 

ani=animation.FuncAnimation(fig,animate,frames=200,interval=0.1,blit=True,init_func=init)

#ani.save('particle_box.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
plt.show()
