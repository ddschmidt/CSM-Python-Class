import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from scipy.spatial.distance import pdist, squareform
import matplotlib.animation as animation


#constants
radiusH=0.1#120*10e-12
massH=1#1.6737*10e-27
box_length=200*radiusH
#this one is actually important, how many particles we put in the box
numparticles=20
#groupsize=6+1
groupdis=20
sepdist=2

weight_align=10
weight_cohesion=10
weight_seperation=1
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
        self.centerofmass=[]
        self.groupvel=[]
        self.neighbormass=0
        self.seperation=[]
        self.cohesion=0
        self.alignment=[]

    def step(self, dt):
        #move time
        self.time_elapsed +=dt

            #move particles
        self.state[:,:3]+=dt*self.state[:,3:6]

    #Simple Collision Detection thing
            #finds distance between particles
                #produces an NxN array of distances between particles
        D=squareform(pdist(self.state[:,:3]))
        V=self.state[:,3:6]
       
        ###################################################
        #Flocking behavior (figure out what birds are within area of influence for each other)
        bird1,bird2 = np.where(D<groupdis*self.radius)
        test = [[bird1[i], bird2[i]] for i in range(len(bird1))]
        
        for i1 in bird1:
            m1=self.state[i1,6]
            r1=np.array(self.state[i1,:3])
            v1=np.array(V[i1])
            self.centerofmass=np.array((0,0,0))
            self.seperation=np.array((0,0,0))
            self.alignment=np.array((0,0,0))
            self.groupvel
            self.neighbormass=1
            for i2 in bird2:
                if i1==i2:
                    pass
                elif [i1,i2] in test:
                #neighbor properties
                    #mass
                    m2=self.state[i2,6]
                    # position
                    r2=np.array(self.state[i2,:3])
                    # velocity
                    v2=np.array(V[i2])
                    
                #neighbor calculations
                    self.neighbormass=self.neighbormass+self.state[i2,6]
                    #center of mass
                    self.centerofmass=self.centerofmass+m2*r2
                    
                    #velocity pointing
                    self.alignment=self.alignment+v2

                    dis=abs(np.linalg.norm(r1-r2))

                    if dis<=sepdist:
                        self.seperation=self.seperation-(r2-r1)
                    
            #cohesion
            com=self.centerofmass/self.neighbormass
            comnorm=np.linalg.norm(com)
            if comnorm==0.0:
                self.cohesion=np.array((0,0,0))
            else:
                self.cohesion=self.cohesion/comnorm
            #print(com)
            #time.sleep(2)
            #self.cohesion=com/np.linalg.norm(com)
            #seperation
            sepnorm=np.linalg.norm(self.seperation)
            if sepnorm==0.0:
                self.seperation=np.array((0,0,0))
            else:
                self.seperation=self.seperation/sepnorm
                
            #self.alignment=self.alignment/np.linalg.norm(self.alignment)
            alignnorm=np.linalg.norm(self.alignment)
            if alignnorm==0.0:
                self.alignment=np.array((0,0,0))
            else:
                self.alignment=self.alignment/alignnorm
  
            #assign new velocities
            velfull=v1+(self.alignment*weight_align+self.cohesion*weight_cohesion+self.seperation*weight_seperation)*dt
            veln=np.linalg.norm(velfull)
            self.state[i1,3:6]=np.ndarray.tolist(velfull/veln)

        ########################################################
            #calculate collisions
        ind1,ind2 = np.where(D<2*self.radius)
        unique = (ind1<ind2)
        ind1 = ind1[unique]
        ind2 = ind2[unique] 
        for i1,i2 in zip(ind1,ind2):
            m1=self.state[i1,6]
            m2=self.state[i2,6]
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
            self.state[i1,3:6]=(v_cm+v_rel*m2/(m1+m2))/np.linalg.norm((v_cm+v_rel*m2/(m1+m2)))
            self.state[i2,3:6]=(v_cm-v_rel*m1/(m1+m2))/np.linalg.norm((v_cm-v_rel*m2/(m1+m2)))
            #print((v_cm+v_rel*m2/(m1+m2)))
        ###############################################################
        #boundary crossing
        #this just makes sure the balls don't clip through the wall
        crossed_x1 = (self.state[:,0]<self.cornors[0] + self.radius)
        crossed_x2 = (self.state[:,0]>self.cornors[1] - self.radius)
        crossed_y1 = (self.state[:,1]<self.cornors[2] + self.radius)
        crossed_y2 = (self.state[:,1]>self.cornors[3] - self.radius)
        crossed_z1 = (self.state[:,2]<self.cornors[4] + self.radius)
        crossed_z2 = (self.state[:,2]>self.cornors[5] - self.radius)

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
        
        for i in range(numparticles):
            veln=np.linalg.norm(self.state[i1,3:6])
            self.state[i1,3:6]=np.ndarray.tolist(self.state[i1,3:6]/veln)

        ###############################################################

################################################################################################

#this creates an initial state of numparticle number of balls with an initial position velocity and mass
np.random.seed()
red=np.random.rand(numparticles,7)-0.5
#makes sure they are in the box
red[:,:3]*=box_length*.1
#red[:,3:6]*=100
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

D=squareform(pdist(box.state[:,:3]))
bird1,bird2 = np.where(D<groupdis*box.radius)
test = [[bird1[i], bird2[i]] for i in range(len(bird1))]
#%% Animation in 3D
################################################################################################
#plots for animation
fig = plt.figure()
ax= p3.Axes3D(fig)

#plots for the particles
particles, = ax.plot([],[],[],'bo',ms=1,animated=True)

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

ani=animation.FuncAnimation(fig,animate,frames=100,interval=0.1,blit=True,init_func=init)

#ani.save('MovWave.mpeg', writer="ffmpeg")

#ani.save('particle_box.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
plt.show()
