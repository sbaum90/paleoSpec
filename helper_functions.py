import numpy as np

from scipy.special import erf

def eta(vmin, vrel, sigv, vesc):
   """ 
   returns mean inverse speed for a Maxwell-Boltzmann
   distribution as a function of vmin for
      vrel - bulk velocity wrt the distribution
      sigv - velocity dispersion
      vesc - escape velocity
   units: enter all speeds in the same units,
      then, eta is returned in the inverse units
   """
   amin = vmin/(np.sqrt(2)*sigv)
   arel = vrel/(np.sqrt(2)*sigv)
   aesc = vesc/(np.sqrt(2)*sigv)
   # prepare format of speeds
   l = 1
   try:
      l = np.max([l,len(amin)])
   except:
      pass
   try:
      l = np.max([l,len(arel)])
   except:
      pass
   try:
      l = np.max([l,len(aesc)])
   except:
      pass
   try:
      if l == len(amin):
         pass
   except:
      amin = amin*np.ones(l)
   try:
      if l == len(arel):
         pass
   except:
      arel = arel*np.ones(l)
   try:
      if l == len(aesc):
         pass
   except:
      aesc = aesc*np.ones(l)
   # calculate the integral
   Nesc = erf(aesc)-2./np.sqrt(np.pi)*aesc*np.exp(-aesc**2)
   N = 1./( 2.*vrel*Nesc )
   vel_int = np.zeros(l)
   vel_int += (aesc<arel)*(amin<np.abs(arel-aesc))*2.*Nesc
   vel_int += (aesc>arel)*(amin<np.abs(arel-aesc))*( erf(amin+arel) - erf(amin-arel) - 4./np.sqrt(np.pi)*arel*np.exp(-aesc**2) )
   vel_int += (np.abs(arel-aesc)<amin)*(amin<(arel+aesc))*( erf(aesc) - erf(amin-arel) - 2./np.sqrt(np.pi)*(arel+aesc-amin)*np.exp(-aesc**2) )
   return N*vel_int

def F2helm(q2, A):
   """
   returns helm form factor.
   input:
      q2 - squared momentum transfer in keV^2
      A - atomic mass number
   """
   hbarc=1.97e5 # keV fm
   s=0.9/hbarc
   a=0.52/hbarc
   c=(1.23*A**(1./3.)-0.6)/hbarc
   rn=np.sqrt(c**2+7./3.*np.pi**2*a**2-5*s**2)
   q=np.sqrt(q2)
   F2 = (3*(np.sin(q*rn)-q*rn*np.cos(q*rn))/(q*rn)**3*np.exp(-q2*s**2/2))**2
   F2 = np.clip(F2, 0, 1e30)
   return F2
