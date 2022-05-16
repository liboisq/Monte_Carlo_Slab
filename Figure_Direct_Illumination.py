#!/usr/bin/python
# -*-coding:utf-8 -*

from pylab import *
from scipy import integrate
#from refractive_index import refice
import sys
sys.path.append("/home/liboisq/RTMC")
from MonteCarlo3D import Photon
import sys
import scipy.special as sc
from numpy import frompyfunc
import pickle


#-------------Latex Style--------------
rc('font',**{'family':'sans-serif','sans-serif':'Computer Modern Sans Serif','size':'26'})
rc('text', usetex=True)
rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}',r'\sansmath']
rc('xtick',labelsize=26) 
rc('ytick',labelsize=26)

# Slab characteristics
H = 2. # m
tau = 2. # optical thickness
le = H/tau # extinction length

# Scattering properties
g = 0      # asymmetry parameter
omega = 1  # single scattering albedo

# Simulation characteristics
mu0 = 0.5 # incidenct angle
N = 10000000 # number of photons launched

# Several photons launched, starting within the medium to avoid useless unscattered photons
# Total distance saved
if 0:  
    all_l = []        # saving all paths
    all_l_ref = []
    all_l_trans = []
    all_l_dir = []
    all_weight = []   # saving all paths weights
    all_weight_ref = []
    all_weight_trans = []
    all_weight_dir = []
    
    for n in range(N):
        events = 0  # number of scattering events to date
        
        # Initial direction 
        muz = mu0
        mux = sqrt(1-mu0**2)
        muy = 0                # phi taken as 0 always
        direction_incident = [mux, muy, muz] 
        
        # starting point within the medium       
        d0 = H*random()  
        l0 = d0/muz      # distance travelled to this point
        #weight = H/muz*1/le*exp(-l0/le) # to account for probability of being scattered at depth z0, in average equals 1-exp(-tau/muz)
        weight = H/(muz*le)*exp(-l0/le) # to account for probability of being scattered at depth z0, cf Monte Carlo estimation of a mean function (1/H is the uniform density of Z)
                                         # in average equals 1-exp(-tau/muz)
        
        # For each photon launched, some complementary part is directly transmitted
        all_l_dir+= [H/muz]
        all_weight_dir+= [exp(-tau/muz)]
        
        # Start at depth
        position_start= [0,0,d0]
    
        # Initialisation
        p0 = Photon(position_start,direction_incident,omega,le,g)
        p0.move_forward()
        direction_start = p0.direction # just to get a direction
        p = Photon(position_start,direction_start,omega,le,g)
              
        # Continue while still in the slab (always true first) and keep only single-scattering
        while p.position[2]>=0 and p.position[2]<=H and events<1: 
            events+=1
            
            # Check only photons that escaped slab in maximum 1 scattering event                 
            p.move_forward()  # updates position and direction, continues if still in the slab
      
        # Terminated -> sorting photons   
        # Reflected photons
        if p.position[2]<0:
            # remove part of the last path travelled outside the slab
            # last_l = lin + lout
            # lin = -last_z/muz
            p.all_distance-= p.last_move
            muz = p.last_direction[2]
            p.all_distance+= -p.last_position[2]/muz
            
            all_l_ref+= [p.all_distance + l0] # L0 to account for distance travelled before first scattering
            all_weight_ref+= [weight]


        # Transmitted photons    
        elif p.position[2]>H:     
            # remove part of the last path travelled outside the slab
            # last_l = lin + lout
            # lin = (H-last_z)/muz
            p.all_distance-= p.last_move
            muz = p.last_direction[2]
            p.all_distance+= (H-p.last_position[2])/muz

            all_l_trans+= [p.all_distance + l0]
            all_weight_trans+= [weight]


    all_l_ref = array(all_l_ref)
    all_l_trans = array(all_l_trans)
    all_l_dir = array(all_l_dir)
    all_weight_dir = array(all_weight_dir)
    all_weight_ref = array(all_weight_ref)
    all_weight_trans = array(all_weight_trans)
    
    pickle.dump(all_l_ref, open("direct_all_l_ref.p", "wb" ))
    pickle.dump(all_l_trans, open("direct_all_l_trans.p", "wb" ))
    pickle.dump(all_l_dir, open("direct_all_l_dir.p", "wb" ))
    pickle.dump(all_weight_ref, open("direct_all_weight_ref.p", "wb" ))
    pickle.dump(all_weight_trans, open("direct_all_weight_trans.p", "wb" ))
    pickle.dump(all_weight_dir, open("direct_all_weight_dir.p", "wb" ))


else:
    all_l_ref = pickle.load(open("direct_all_l_ref.p", "rb" ))
    all_l_trans = pickle.load(open("direct_all_l_trans.p", "rb" ))
    all_l_dir = pickle.load(open("direct_all_l_dir.p", "rb" ))
    all_weight_ref = pickle.load(open("direct_all_weight_ref.p", "rb" ))
    all_weight_trans = pickle.load(open("direct_all_weight_trans.p", "rb" ))
    all_weight_dir = pickle.load(open("direct_all_weight_dir.p", "rb" ))


N_dir = sum(all_weight_dir)
N_ref = sum(all_weight_ref)
N_trans = sum(all_weight_trans)

print("N_trans",N_trans)
print("N_ref",N_ref)
print("N_dir",N_dir)
print("N sum",N_trans+N_dir+N_ref)

pdir = exp(-tau/mu0)
print("-----------------")
print("Proportion of directly transmitted")
print("Theory",pdir)
print("Simulation",N_dir/N)
print("-----------------")

x0 = np.linspace(H/mu0+1e-4,10*(H/mu0),10000)
x00 = np.linspace(H+1e-4,H/mu0-1e-4,1000)
p1t = mu0/(2*le)*np.exp(-x0/le)*np.heaviside(x0-H/mu0,1)*(-np.log(1-H/(mu0*x0))-H/(mu0*x0))
p2t = mu0/(2*le)*np.exp(-x00/le)*np.heaviside(x00-H,1)*np.heaviside(H/mu0-x00,1)*(-np.log((H-x00*mu0)/((1-mu0)*x00))+(x00-H)/(mu0*x00))
ptrans = integrate.trapz(p1t,x0) + integrate.trapz(p2t,x00)

print("Proportion of transmitted, 1 scattering")
print("Theory",ptrans)
print("Simulation",N_trans/N)
print("-----------------")

x1 = np.linspace(0+1e-15,H*(1+1/mu0)-1e-15,10000)
x2 = np.linspace(H*(1+1/mu0)+1e-15,10*H*(1+1/mu0),10000)

p1r1 = mu0/(2*le)*np.exp(-x1/le)*np.heaviside(x1,1)*np.heaviside(H*(1+1/mu0)-x1,1)*(1/mu0+np.log(mu0/(1+mu0)))
pref1 = integrate.trapz(p1r1,x1)                   
p1r2 = mu0/(2*le)*np.exp(-x2/le)*np.heaviside(x2-H*(1+1/mu0),1)*(H/(mu0*x2-H)+np.log(1-H/(mu0*x2)))
pref2 = integrate.trapz(p1r2,x2)

print("Proportion of reflected, 1 scattering")
print("Theory",pref1+pref2)
print("Simulation",N_ref/N)
print("-----------------")


fig1,(ax1,ax2) = subplots(2,1,sharex=True,figsize=(14,12))
fig1.subplots_adjust(top=0.93)
fig1.suptitle(r"$H=%s~\mathrm{m} \quad \tau=%s \quad \mu_0=%s$"%(H,tau,mu0))
xlim(0,5*H)
dx = 5*H/500
ax2.hist(all_l_ref, weights= all_weight_ref/N/dx, bins=500,range=(0,5*H),density=False,color="0.8")
#ax1.plot(x0,N/N_ref_short*pr0,"r")
ax2.plot(x1,p1r1,"k")
ax2.plot(x2,p1r2,"k")   
ax2.set_yscale('log')
#ax2.set_xlim(0,3)
#ax1.plot(x1,N/N_ref_long*(pr3),"k")

#pt1 = 1./(2*x**4*l)*exp(-x/l)*(L**4*log(x/L-1)-x**4*log(1-L/x)+(L**3*x - L*x**3))
ax1.hist(all_l_trans, weights=all_weight_trans/N/dx,bins=500,range=(0,5*H),density=False,color="0.8")

ax1.plot(x0,p1t,"k")
ax1.plot(x00,p2t,"k")
#ax2.set_xlim(0,3)

#ax3.hist(all_l_dir, weights=all_weight_dir,bins=100,range=(0,5*L),density=True)
#ax3.plot(x0,pdir*N/N_dir)
ax2.set_xlabel(r"$l\mathrm{(m)}$",size=26)
ax1.set_ylabel(r"$p_{1,t}(l)$",size=26)
ax2.set_ylabel(r"$p_{1,r}(l)$",size=26)
#fig1.savefig("Figures/Direct_Illumination_PDF.pdf",format = "pdf", dpi=300)

show()


    
  