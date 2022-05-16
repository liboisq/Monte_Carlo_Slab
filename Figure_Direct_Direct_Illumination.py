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
tau = 1. # optical thickness
le = H/tau # extinction length

# Scattering properties
g = 0      # asymmetry parameter
omega = 1  # single scattering albedo

# Simulation characteristics
theta0 = 10./180*pi
mu0 = 0.9 # incidenct angle
s0 = 1./mu0
theta = 10./180*pi # 10 for 10 M points
mur = 0.8
mut = 0.8
theta=arccos(mur)
sr = 1./mur
st = 1./mut
dmu = 0.02
N = 100000000 # number of photons launched

# Several photons launched, starting within the medium to avoid useless unscattered photons
# Total distance saved
if 0:  
    all_l = []        # saving all paths
    all_l_ref1 = []
    all_l_ref2 = []
    all_l_trans1 = []
    all_l_trans2 = []
    all_l_dir = []
    all_weight = []   # saving all paths weights
    all_weight_ref1 = []
    all_weight_ref2 = []
    all_weight_trans1 = []
    all_weight_trans2 = []
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
#        print(p.all_distance)
#        input()
        # Continue while still in the slab (always true first) and keep only single-scattering
        while p.position[2]>=0 and p.position[2]<=H and events<2: 
            events+=1
            
            # Check only photons that escaped slab in maximum 1 scattering event                 
            p.move_forward()  # updates position and direction, continues if still in the slab
      
        # Terminated -> sorting photons   
        # Reflected photons      
        if p.position[2]<0:
            if p.last_direction[2]>-mur-dmu and p.last_direction[2]<-mur+dmu:
                p.all_distance-= p.last_move
                muz = p.last_direction[2]
                p.all_distance+= -p.last_position[2]/muz
                if events ==1 :
                    all_l_ref1+= [p.all_distance + l0] # L0 to account for distance travelled before first scattering
                    all_weight_ref1+= [weight]
                elif events == 2:    
                    all_l_ref2+= [p.all_distance + l0] # L0 to account for distance travelled before first scattering
                    all_weight_ref2+= [weight]

        # Transmitted photons    
        elif p.position[2]>H:    
            if p.last_direction[2]>mut-dmu and p.last_direction[2]<mut+dmu:
                p.all_distance-= p.last_move
                muz = p.last_direction[2]
                p.all_distance+= (H-p.last_position[2])/muz
                if events == 1:
                    all_l_trans1+= [p.all_distance + l0]
                    all_weight_trans1+= [weight]
                elif events == 2:
                    all_l_trans2+= [p.all_distance + l0]
                    all_weight_trans2+= [weight]



    all_l_ref1 = array(all_l_ref1)
    all_l_trans1 = array(all_l_trans1)
    all_l_ref2 = array(all_l_ref2)
    all_l_trans2 = array(all_l_trans2)    
    all_l_dir = array(all_l_dir)
    all_weight_dir = array(all_weight_dir)
    all_weight_ref1 = array(all_weight_ref1)
    all_weight_trans1 = array(all_weight_trans1)
    all_weight_ref2 = array(all_weight_ref2)
    all_weight_trans2 = array(all_weight_trans2)    
    
    pickle.dump(all_l_ref1, open("direct_direct_all_l_ref1_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "wb" ))
    pickle.dump(all_l_trans1, open("direct_direct_all_l_trans1_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "wb" ))
    pickle.dump(all_l_ref2, open("direct_direct_all_l_ref2_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "wb" ))
    pickle.dump(all_l_trans2, open("direct_direct_all_l_trans2_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "wb" ))
    pickle.dump(all_l_dir, open("direct_direct_all_l_dir_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "wb" ))
    pickle.dump(all_weight_ref1, open("direct_direct_all_weight_ref1_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "wb" ))
    pickle.dump(all_weight_trans1, open("direct_direct_all_weight_trans1_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "wb" ))
    pickle.dump(all_weight_ref2, open("direct_direct_all_weight_ref2_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "wb" ))
    pickle.dump(all_weight_trans2, open("direct_direct_all_weight_trans2_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "wb" ))
    pickle.dump(all_weight_dir, open("direct_direct_all_weight_dir_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "wb" ))


else:
    all_l_ref1 = pickle.load(open("direct_direct_all_l_ref1_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "rb" ))
    all_l_trans1 = pickle.load(open("direct_direct_all_l_trans1_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "rb" ))
    all_l_ref2 = pickle.load(open("direct_direct_all_l_ref2_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "rb" ))
    all_l_trans2 = pickle.load(open("direct_direct_all_l_trans2_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "rb" ))
    all_l_dir = pickle.load(open("direct_direct_all_l_dir_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "rb" ))
    all_weight_ref1 = pickle.load(open("direct_direct_all_weight_ref1_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "rb" ))
    all_weight_trans1 = pickle.load(open("direct_direct_all_weight_trans1_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "rb" ))
    all_weight_ref2 = pickle.load(open("direct_direct_all_weight_ref2_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "rb" ))
    all_weight_trans2 = pickle.load(open("direct_direct_all_weight_trans2_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "rb" ))
    all_weight_dir = pickle.load(open("direct_direct_all_weight_dir_mu0_%s_mur_%s_N_%s.p"%(mu0,mur,int(N)), "rb" ))


N_dir = sum(all_weight_dir)
N_ref1 = sum(all_weight_ref1)
N_trans1 = sum(all_weight_trans1)
N_ref2 = sum(all_weight_ref2)
N_trans2 = sum(all_weight_trans2)

print("N_trans1",N_trans1)
print("N_ref1",N_ref1)
print("N_trans2",N_trans2)
print("N_ref2",N_ref2)
print("N_dir",N_dir)


x1r = np.linspace(1e-4,H*(s0+sr)-1e-4,10000)
x1t = np.linspace(H*min(s0,st)+1e-4,H*max(s0,st)-1e-4,10000)
x2r = np.linspace(H*(s0+sr)+1e-4,10*(H/mu0),10000)
x2t = np.linspace(H*(1+s0+st)+1e-4,10*(H/mu0),10000)
x2t1 = np.linspace(H*(s0+st)+1e-4,H*(1+s0+st)-1e-4,10000)
x2r1 = np.linspace(1e-4,2*H-1e-4,10000)
x2r1 = np.linspace(1e-4,min(H*(s0+1),H*(sr+1))-1e-4,10000)
x2r2 = np.linspace(2*H+1e-4,min(H*(s0+1),H*(sr+1))-1e-4,10000)

#x2t = x2r
print(H*(1/mu0 + 1/mur),H*max(s0,st))
Ar = 1./(4*le**2)*exp(-x2r/le)
At = 1./(4*le**2)*exp(-x2t/le)
Ar1 = 1./(4*le**2)*exp(-x2r1/le)*s0

print(sin(theta))
p1t = 1./(2*mu0*le)*1./abs(st-s0)*exp(-x1t/le)*sin(theta)
#p1t = 1./(mu0*le)*1./abs(st-s0)*exp(-x1t/le)
p1r = 1./(2*mu0*le)*1./(sr+s0)*exp(-x1r/le)*sin(theta)#sin(theta)*

p2t = At*sin(theta)*mut*(-(x2t-H*s0)*log(x2t-H*s0)-(x2t-H*st)*log(x2t-H*st)+(x2t-H*(st+s0))*log(x2t-H*(st+s0))+x2t*log(x2t))
p2r = Ar*sin(theta)*mur*(-(x2r-H*s0)*log(x2r-H*s0)-(x2r-H*sr)*log(x2r-H*sr)+(x2r-H*(sr+s0))*log(x2r-H*(sr+s0))+x2r*log(x2r))
d1 = x2r1/(s0+sr)
dm = x2r1/(1+s0)
d12 = x2r2/(s0+sr)
d22 = (x2r2-H*(sr+1))/(s0-1)
dm2 = x2r2/(1+s0)

#p2r1 = Ar1*sin(theta)*mur*(d1 * log(1+mur)+(dm-d1)*log(1-mur) \
#       + mu0/(mu0+mur)*(-mur*x2r1*(log(mur*x2r1)-1)) + mu0*x2r1*(log(x2r1)-1) \
#       - mu0/(mu0+mur)*(1-mur)/(1+mu0)*mu0*x2r1*(log((1-mur)/(1+mu0)*mu0*x2r1)-1) \
#       - mu0*mu0*x2r1/(1+mu0)*(log(mu0*x2r1/(1+mu0))-1))
##
#p2r2 = Ar1*sin(theta)*mur*((d12-d22)*log(1+mur)+(dm2-d22)*log(1-mur) \
#        + mu0*(mu0/mur*(H*(1-mur/mu0)-mur*x2r2)/(1-mu0)*(log(mu0/mur*(H*(1-mur/mu0)-mur*x2r2)/(1-mu0))-1) - (x2r2-H*sr)*(log(x2r2-H*sr)-1)) )
##        + mu0*x2r2*(log(x2r2)-1) \
#        - mu0/(mu0+mur)*(H*(1+mu0/mur)-mu0*x2r2*(1+mur))/(1-mu0)*(log((H*(1+mu0/mur)-mu0*x2r2*(1+mur))/(1-mu0))-1) \
#        - mu0/(mu0+mur)*(1-mur)/(1+mu0)*mu0*x2r2*(log((1-mur)/(1+mu0)*mu0*x2r2)-1) \
#        - mu0*mu0*x2r2/(1+mu0)*(log(mu0*x2r2/(1+mu0))-1) )

p2t1 = s0*At*sin(theta)*mut*(H*log(1+mut)-1./s0*((x2t1-H*s0)*log(x2t1-H*s0)-x2t1*log(x2t1))-1./(1-mut*s0)*(mut*(x2t1-H*s0)*log(mut*(x2t1-H*s0))-(mut*x2t1-H)*log(mut*x2t1-H)))
       
print(p2t1)
print(amax(log(x2r-H*s0)))
print(amax(log(x2r-H*sr)))
print(amax(log(x2t-H*(st+s0))))


fig1,(ax1,ax2) = subplots(2,1,sharex=True,figsize=(13,12))
fig1.subplots_adjust(top=0.93)
fig1.suptitle(r"$H=%s~\mathrm{m} \quad \tau=%s \quad \mu_0=%s \quad \mu_r=\mu_t=%s$"%(H,tau,mu0,mur))

nbins = 200
lmax = 5*H
dx = lmax/nbins
xlim(0,lmax)

#ax2.hist(all_l_ref1, weights= all_weight_ref1/N/dx, bins=nbins,range=(0,lmax),density=False,color="0.5")
ax2.hist(all_l_ref2, weights= all_weight_ref2/N/dx, bins=nbins,range=(0,lmax),density=False,color="0.8")
#ax2.set_yscale('log')
#ax1.plot(x0,N/N_ref_short*pr0,"r")
#ax2.plot(x1r,2*dmu/sin(theta)*p1r,"k") # p(theta)dtheta = p(mu) d(mu) -> 
ax2.plot(x2r,2*dmu/sin(theta)*p2r,"k")   
#ax2.plot(x2r1,2*dmu/sin(theta)*p2r1,"r") 
#ax2.plot(x2r2,2*dmu/sin(theta)*p2r2,"r") 
#ax2.set_yscale('log')
#ax2.set_xlim(0,3)
#ax1.plot(x1,N/N_ref_long*(pr3),"k")

#pt1 = 1./(2*x**4*l)*exp(-x/l)*(L**4*log(x/L-1)-x**4*log(1-L/x)+(L**3*x - L*x**3))

#ax1.hist(all_l_trans1, weights=all_weight_trans1/N/dx,bins=nbins,range=(0,lmax),density=False,color="0.5")
ax1.hist(all_l_trans2, weights=all_weight_trans2/N/dx,bins=nbins,range=(0,lmax),density=False,color="0.8")
#ax1.set_yscale('log')
#ax1.plot(x1t,2*dmu/sin(theta)*p1t,"k")
ax1.plot(x2t,2*dmu/sin(theta)*p2t,"k")
#ax1.plot(x2t1,2*dmu/sin(theta)*p2t1,"r")
#ax2.set_xlim(0,3)

#ax3.hist(all_l_dir, weights=all_weight_dir,bins=100,range=(0,5*L),density=True)
#ax3.plot(x0,pdir*N/N_dir)
ax2.set_xlabel(r"$l\mathrm{(m)}$",size=26)
ax1.set_ylabel(r"$p_{2,t}(l)$",size=26)
ax2.set_ylabel(r"$p_{2,r}(l)$",size=26)
fig1.savefig("Figures/Direct_Direct_Illumination_PDF.pdf",format = "pdf", dpi=300)

show()


    
  
