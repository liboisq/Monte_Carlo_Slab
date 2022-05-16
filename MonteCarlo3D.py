#!/usr/bin/python
# -*-coding:utf-8 -*

from pylab import *
#from refractive_index import refice
import sys

sys.path.append("/home/qlibois/These/Tartes/tartes")
#from tartes import tartes

#-------------Latex Style--------------
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Times'],'size':'25'})
#--------------------------------------

rc('xtick',labelsize=25) 
rc('ytick',labelsize=25) 

# le domaine 3D semi-infini s'étend sur -inf<x<inf, -inf<y<inf et z>=0. Un photon sort s'il 

c = 3e8 # light speed  
ni = 1.32 # ice index in the visible
theta = linspace(-pi,pi,1000)

class Photon: 
  def __init__(self,position,direction,omega,l,g):
      self.position = position
      self.last_position = position
      self.omega = omega	
      self.l = l	
      self.direction = direction # mux, muy, muz
      self.last_direction = direction
      self.all_distance = 0 #sqrt(position[0]**2+position[1]**2+position[2]**2) # valid only for entrance at (0, 0, 0)
      self.g = g
      self.last_move = 0
    
  def move_forward(self):      
      # Starting from a given point with a given direction, update next direction
      
      # 1 - Draw a distance to next scattering event
      # Explanations: p(l) = 1/lint exp(-l/lint) has cumulative distribution 1-exp(-l/lint)
      # inverse of this is lint ln(1-y)
      # uniform y gives x with exponential distribution   
      l_true = -log(random())*self.l #distribution des step sizes   
      self.last_move = l_true
      self.last_direction = self.direction
      
      # New position
      new_pos = [self.position[0]+l_true*self.direction[0],self.position[1]+l_true*self.direction[1],self.position[2]+l_true*self.direction[2]] 
      
      # Select a new scattering direction with Henyey-Greenstein phase function
      a = random()   

      # Numerical approach if g != 0
#      hg = 1./(4*pi)*(1-self.g**2)/(1+self.g**2-2*self.g*cos(theta))**(3./2)
#      cum = trapz(hg,theta)
#      cum/=cum[-1]
#      k = searchsorted(cum,a)
      
      # Analytical approach = inverse of Henyey-Greenstein distribution
      if self.g == 0:
          costheta = 1-2*a
      else:    
          costheta = 1./(2*self.g)*(1+self.g**2-((1-self.g**2)/(1-self.g+2*self.g*a))**2) # -g or +g in last term?
      
      sintheta = sqrt(1-costheta**2)
#      t = theta[k] #choix de l'angle de diffusion
      
      phi = 2*pi*random()      # isotropic equivalent
      
      
      # Find new direction after this event   
      if self.direction[2] == 1:
          mux = sintheta*cos(phi)
          muy = sintheta*sin(phi)
          muz = costheta
#          print("muz",muz)
	
      elif self.direction[2] == -1:
          mux = sintheta*cos(phi)
          muy = sintheta*sin(phi) # wang et al., 1995 different here than Wikipedia (and that's correct)
          muz = -costheta
#          print("muz",muz)
          
      else:
          mux = sintheta*(self.direction[0]*self.direction[2]*cos(phi)-self.direction[1]*sin(phi))/sqrt(1-self.direction[2]**2)+self.direction[0]*costheta
          muy = sintheta*(self.direction[1]*self.direction[2]*cos(phi)+self.direction[0]*sin(phi))/sqrt(1-self.direction[2]**2)+self.direction[1]*costheta
#          print(self.direction[2],sintheta,cos(phi),costheta)
          muz = -sqrt(1-self.direction[2]**2)*sintheta*cos(phi)+self.direction[2]*costheta 
#          mux = sintheta*cos(phi)
#          muy = sintheta*sin(phi)
#          muz = costheta
          
      # checked that mux**2+muy**2+muz**2=1
      
      self.last_position = self.position
      self.position = new_pos
      self.direction = [mux,muy,muz]
      self.all_distance+= l_true


if __name__ == "__main__":
    #1-Albédo spectral
    g1 = 0.86
    g = 0.86
    B = 1.6
    rho_ice=917.
    rho=300.
    SSA=15.
    sigma_e=rho*SSA/2
    l0=1./(sigma_e*(1-g1))
    l0=1./(sigma_e)
    l=1./(sigma_e)
    #l=1./(sigma_e*(1-g1))
    lint = 4*B/(rho_ice*SSA)
    direction_start=[0,0,1]
    position_start=[0,0,l0]
    
    # HG random choice
    
    
    
    all_time=loadtxt("all_time_ref.dat")
    # Patterson 1989
    g1 = 0
    wl = 800
    ki = refice(wl*1e-9)[1]
    gamma = 4*pi*ki/(wl*1e-9) 
    mus = rho*SSA/2
    musg = (1-g1)*mus
    mua = B*gamma*rho/rho_ice
    c = c*(4*B+2*rho_ice/rho)/(4*B*ni+2*rho_ice/rho)
    print(mua,musg)
    
    mua = 17.6
    musg = 850
    c = 3e8
    
    D = 1./(3*(mua+musg))
    z0 = 1./musg
    
    t_patt = linspace(1e-12,550e-12,5000)
    r_patt = (4*pi*D*c)**(-1./2)*z0*t_patt**(-3./2)*exp(-mua*c*t_patt)*exp(-z0**2/(4*D*c*t_patt))
    integ = trapz(r_patt,t_patt)
    
    semilogy(1e12*t_patt,1e-12*r_patt,'r-')
    x,b = histogram(1e12*all_time,bins=arange(0,550,10),density='True')
    semilogy(b[1:],x)
    #gca().set_yscale("log")
    ylim(1e-6,0.1)
    xlim(0,550)
    show()
        
    #
    #fig=figure(15,figsize=(10,7))
    #
    #plot(wls,alb,'o-',color="0",linewidth=2,label=r"Monte Carlo") 
    #plot(wls,exp(-8*9./7*sqrt(B*4*pi*refice(wls*1e-9)[1]/(wls*1e-9)/(3*rho_ice*SSA*2*(1-g)))),'o-',color="0.5",linewidth=2,label=r"Th\'eorie ART")
    #xlabel(r"Longueur d'onde (nm)",size=30)
    #ylabel(r"Alb\'edo",size=30)
    #legend(loc=0,numpoints=1)
    #show()
    
    #fig.savefig("/home/qlibois/These/Ecrit/Manuscrit/Figures/SnowOptics/albedo_MC.pdf",dpi=300,format="pdf")
    
    tt=1
    if tt:
        wl = 500
        alb=[]
        ki=0.1*refice(wl*1e-9)[1]
        print(wl)
        gamma=4*pi*ki/(wl*1e-9) 
        omega=1-2*B*gamma/(SSA*rho_ice)
        out=0.
        N=20000 # number of photons to be launched
        all_time = []
        all_l = []
        all_l2 = []
    #    mua = 17.6
    #    musg = 850
        mua = 0.
        mus = 200.
        g = 0.2
    #    g1 = 0.7
    #    mus = musg/(1-g1)
        mue = mus + mua
    #    mue_norm = mus*(1-g1) + mua
        l=1./(mue)
    
        lint = 0
        z0 = 1./musg
        omega=mus/mue
    #    position_start=[0,0,l]
        direction_start=[0,0,1]
        
        for n in range(N):
          time = 0  
          print(n)
          l0 =-log(random())*l
          position_start=[0,0,l0]
          p=Photon(position_start,direction_start,omega,l0,lint,g)
          #traj=[position_start]
          while random()>(1-omega) and p.position[2]>=0 and p.position[2]<=0.05 :
              p.move_forward()  
              #traj+=[p.position]
          if p.position[2]<0:
              out+=1
              all_time+=[p.tms]
              all_l+=[p.all_distance]
          if p.position[2]<0 or p.position[2]>0.05:
              all_l2+=[p.all_distance]
              
        alb+=[out/N] 
          #traj=array(traj)
          #plot(traj[:,0],traj[:,1],'-')
        
        savetxt("all_time_ref.dat",array(all_time))
        savetxt("all_l_ref.dat",array(all_l))
        print("Mean distance reflected", mean(all_l))
        print("Mean distance all", mean(all_l2))
        hist(all_time)
        print(mean(tms))
        show()
    
    
    tt=0
    if tt==1:
      wls=arange(600,1320,2)
      alb=[]
      for wl in wls:
        ki=refice(wl*1e-9)[1]
        print(wl)
        gamma=4*pi*ki/(wl*1e-9) 
        omega=1-2*B*gamma/(SSA*rho_ice)
        out=0.
        N=1000 # number of photons to be launched
        all_time = []
        for n in range(N): 
          print(wl,n)
          p=Photon(position_start,direction_start,omega,l,lint)
          #traj=[position_start]
          while random()>(1-omega) and p.position[2]>=0:
              p.move_forward()  
              #traj+=[p.position]
          if p.position[2]<0:
              out+=1
    
        alb+=[out/N] 
          #traj=array(traj)
          #plot(traj[:,0],traj[:,1],'-')
        
      savetxt("wls2_MC.dat",wls)
      savetxt("alb2_MC.dat",alb)
      plot(wls,alb,label="Monte Carlo") 
      plot(wls,exp(-8*9./7*sqrt(B*4*pi*refice(wls*1e-9)[1]/(wls*1e-9)/(3*rho_ice*SSA*2*(1-g)))),label="Kokha")
      legend(loc=0)
      show()
    
    #2-Transmittance
    trans_mc=loadtxt("trans_MC.dat")
    eps=range(1,11)
    wl=700
    
    tt=0
    if tt==1:
      ki=refice(wl*1e-9)[1]
      gamma=4*pi*ki/(wl*1e-9) 
      omega=1-2*B*gamma/(SSA*rho_ice)
      trans_mc=[]
      N=10000
      eps=range(1,11)
      for d in eps:
        trans=0
        for n in range(N):
          print(d,n)
          p=Photon(position_start,direction_start,omega,l,all_distance)
          #traj=[position_start]
          while random()>(1-omega) and p.position[2]<=0.01*d and p.position[2]>=0:
            	p.move_forward()  
            	if p.position[2]>0.01*d:
            	  trans+=1
    	  
        trans_mc+=[1./N*trans]
          
      plot(eps,trans_mc)  
      show()
      savetxt("trans_MC.dat",trans_mc) 
      
        
    #transmission tartes
    
    SSA=array([SSA])
    dens=array([rho])
    wavelengths=array([wl])
    l=len(SSA)
    
    fdiff=zeros_like(wavelengths)
    fdir=ones_like(wavelengths)
    mudir=cos(0./180*pi)
    soilalbedo=zeros_like(wavelengths)
    
    #Valeurs de Kokha 2004
    #B0=1.225*ones([l]) #Kokha p.61,69
    B0=B*ones([l])
    g0=g1*ones([l])
    y0=0.728*ones([l])
    W0=0.0611*ones([l])
    Impurities=list()
    for k in range(l):
        impurities=dict()
        impurities["soot"]=[1800.,0]
        Impurities+=[impurities]
          
    trans_tartes=[]
    for d in eps:
      depth=array([0.01*d])	
      albedo,NRJ,soilabs=tartes(SSA,dens,depth,g0,y0,W0,B0,Impurities,soilalbedo,wavelengths*1e-9,fdiff,fdir,mudir)
      trans_tartes+=[soilabs/(fdir*mudir)]
      
    fig=figure(15,figsize=(10,7))
    
    plot(eps,100*array(trans_mc),'o-',color="0",linewidth=2,label=r"Monte Carlo") 
    plot(eps,100*array(trans_tartes),'o-',color="0.5",linewidth=2,label=r"$\delta$-Eddington")
    xlabel(r"Epaisseur (cm)",size=30)
    ylabel(r"Transmittance ($\%$)",size=30)
    legend(loc=0,numpoints=1)
    show()
    
    fig.savefig("/home/qlibois/These/Ecrit/Manuscrit/Figures/SnowOptics/trans_MC.pdf",dpi=300,format="pdf")
    
    wl=700.
    ki=refice(wl*1e-9)[1]
    gamma=4*pi*ki/(wl*1e-9) 
    omega=1-2*B*gamma/(SSA*rho_ice)
    
    c=16./(3*(1-g)*sigma_e)*(9./7)**2
    lm=linspace(0,0.05,100)
    proba=sqrt(c)/(2*sqrt(pi))*lm**(-3./2)*exp(-c/(4*lm))
    
    plot(lm,proba)
    #show()
    
    position_start=[0,0,l0]
    
    N=1000
    diff=[]
    for n in range(N):
      print(n)
      p=Photon(position_start,direction_start,omega,l,all_distance)
      pos_temp=position_start
      l_diff=0  
      while random()>(1-omega) and p.position[2]>=0:
        p.move_forward()   
        if p.position[2]<0 or p.all_distance>0.5:
          diff+=[p.all_distance]
        
    b=linspace(0,0.05,100) 
    val,e=histogram(diff,bins=b,density=True)
    
    plot((b[1:]+b[:-1])/2,val)
    show()
    