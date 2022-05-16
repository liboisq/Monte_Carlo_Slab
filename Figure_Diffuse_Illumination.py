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
H = 1. # m
tau = 1. # optical thickness
le = H/tau # extinction length

# Scattering properties
g = 0      # asymmetry parameter
omega = 1  # single scattering albedo

# Simulation characteristics
N = 10000000 # number of photons launched

# Multiple photons launched, starting within the medium to avoid useless unscattered photons
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
        events = 0
        # Isotropic illumination
        # p(theta) = 2 cos(theta)*sin(theta) = sin(2 theta)
        # cumul = 1/2 (1- cos(2 theta)) = sin2(theta) = 1- cos2(theta)
        # donc theta = acos(sqrt(1-random())) = acos(sqrt(random()))   
        mu0 = sqrt(random()) # diffuse illumination
        
        #  # Initial direction 
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
    
    pickle.dump(all_l_ref, open("diffuse_all_l_ref.p", "wb" ))
    pickle.dump(all_l_trans, open("diffuse_all_l_trans.p", "wb" ))
    pickle.dump(all_l_dir, open("diffuse_all_l_dir.p", "wb" ))
    pickle.dump(all_weight_ref, open("diffuse_all_weight_ref.p", "wb" ))
    pickle.dump(all_weight_trans, open("diffuse_all_weight_trans.p", "wb" ))
    pickle.dump(all_weight_dir, open("diffuse_all_weight_dir.p", "wb" ))


else:
    all_l_ref = pickle.load(open("diffuse_all_l_ref.p", "rb" ))
    all_l_trans = pickle.load(open("diffuse_all_l_trans.p", "rb" ))
    all_l_dir = pickle.load(open("diffuse_all_l_dir.p", "rb" ))
    all_weight_ref = pickle.load(open("diffuse_all_weight_ref.p", "rb" ))
    all_weight_trans = pickle.load(open("diffuse_all_weight_trans.p", "rb" ))
    all_weight_dir = pickle.load(open("diffuse_all_weight_dir.p", "rb" ))


N_dir = sum(all_weight_dir)
N_ref = sum(all_weight_ref)
N_trans = sum(all_weight_trans)

print("N_trans",N_trans)
print("N_ref",N_ref)
print("N_dir",N_dir)
print("N sum",N_trans+N_dir+N_ref)

nmax = 5
xall = linspace(0,100,10000)

pdir = 2*H**2*exp(-xall/le)/xall**3
pdir = integrate.trapz(pdir,xall)

print("-----------------")
print("Directly transmitted")
print(pdir)
print(N_dir/N)
print("-----------------")

x0 = np.linspace(H+1e-5,10*H,10000)
p1t = 2./(3*x0**3*le)*exp(-x0/le)*(H**3*log(x0/H-1)-x0**3*log(1-H/x0) + H**2*x0 - H*x0**2)

ptrans = integrate.trapz(p1t,x0)

print("Transmitted, 1 scattering")
print(ptrans)
print(N_trans/N)
print("-----------------")

x1 = np.linspace(0+1e-15,2*H-1e-15,10000)
x2 = np.linspace(2*H+1e-15,10*H,10000)

p1r1 = 1./le*exp(-x1/le)*2./3*(1-log(2))
pref1 = integrate.trapz(p1r1,x1)
                        
#p1r2 = 1./l*exp(-x2/l)*1./3*(2*log(1-L/x2)+4*(L/x2)**3*log(x2/L-1) + 2*L/x2**3*(1./2*((x2-L)**2-(L**2/(x2-L))**2) + 2*L*(x2**2/(x2-L))) -1./2*((x2/(x2-L))**2-1) + 2*L/(x2-L) + 3./2*(L/(x2-L))**2)   
p1r2 = 1./le*exp(-x2/le)*2./3*(log(1-H/x2)+2*(H/x2)**3*log(x2/H-1) + H/x2**3*(1./2*((x2-H)**2-(H**2/(x2-H))**2) + \
                             2*H*(x2-H*x2/(x2-H))) -1./4*((x2/(x2-H))**2-1) + H/(x2-H) + 3./4*(H/(x2-H))**2)

pref2 = integrate.trapz(p1r2,x2)

print("Reflected, 1 scattering")
print(pref1+pref2)
print(N_ref/N)
print("-----------------")

fig1,(ax1,ax2) = subplots(2,1,sharex=True,figsize=(13,12))
fig1.subplots_adjust(top=0.93)
fig1.suptitle(r"$H=%s~\mathrm{m} \quad \tau=%s$"%(H,tau))
xlim(0,5*H)
dx = 5*H/500
ax2.hist(all_l_ref, weights= all_weight_ref/N/dx, bins=500,range=(0,5*H),density=False,color="0.8")
#ax1.plot(x0,N/N_ref_short*pr0,"r")
ax2.plot(x1,p1r1,"k")
ax2.plot(x2,p1r2,"k")   
#ax1.plot(x1,N/N_ref_long*(pr3),"k")

#pt1 = 1./(2*x**4*l)*exp(-x/l)*(L**4*log(x/L-1)-x**4*log(1-L/x)+(L**3*x - L*x**3))
ax1.hist(all_l_trans, weights=all_weight_trans/N/dx,bins=500,range=(0,5*H),density=False,color="0.8")
ax1.plot(x0,p1t,"k")
ax2.set_xlabel(r"$l\mathrm{(m)}$",size=26)
ax1.set_ylabel(r"$p_{1,t}(l)$",size=26)
ax2.set_ylabel(r"$p_{1,r}(l)$",size=26)

#ax3.hist(all_l_dir, weights=all_weight_dir,bins=100,range=(0,5*L),density=True)
#ax3.plot(x0,pdir*N/N_dir)

fig1.savefig("Figures/Diffuse_Illumination_PDF.pdf",format = "pdf", dpi=300)
show()



    
  