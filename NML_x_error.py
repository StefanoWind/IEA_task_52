# -*- coding: utf-8 -*-
"""
Created on Wed Sep 24 12:07:19 2025

@author: sletizia
"""
import os
cd=os.path.dirname(__file__)
from matplotlib import pyplot as plt
import numpy as np


#%% Inputs
alpha=15 #[deg] lidar half-opening angle
Uhor=10 #[m/s]  #horizontal wind speed at hub height
theta=0 #[deg] #relative wind direction 
ShearBins=np.arange(-0.05,0.41,0.05) #shear exponents
DeltaZBeamsBins=np.arange(-40,41,2) #[m] height differences
ZHH=60 #[m] hub height

#%% Initialization
S,DZ=np.meshgrid(ShearBins,DeltaZBeamsBins)

#%% Main

#samples Uhor a the two beams
Uhor_Los1=(1-DZ/2/ZHH)**S*Uhor
Uhor_Los2=(1+DZ/2/ZHH)**S*Uhor

#LOS velocities
VLos1=Uhor_Los1*np.cos(np.radians(alpha-theta))
VLos2=Uhor_Los2*np.cos(np.radians(alpha+theta))

Vx_WFR=(VLos1+VLos2)/(2*np.cos(np.radians(alpha)))
Vy_WFR=(VLos1-VLos2)/(2*np.sin(np.radians(alpha)))

Vh_WFR=(Vx_WFR**2+Vy_WFR**2)**0.5

DeltaWS=Vh_WFR-Uhor

#error-free LOS velocities
VLos1_0=Uhor*np.cos(np.radians(alpha-theta))
VLos2_0=Uhor*np.cos(np.radians(alpha+theta))

#compact WFR form
c=np.cos(np.radians(alpha))
s=np.sin(np.radians(alpha))
A=1/c**2+1/s**2
B=2*(1/c**2-1/s**2)

Vh_WFR2=0.5*(A*(VLos1**2+VLos2**2)+B*VLos1*VLos2)**0.5

#first order derivatives
dVh_dVLos1=1/8*(2*A*VLos1_0+B*VLos2_0)/Uhor
dVh_dVLos2=1/8*(2*A*VLos2_0+B*VLos1_0)/Uhor

#second order derivatives
d2Vh_dVLos1_2=A/(4*Uhor)-(2*A*VLos1_0+B*VLos2_0)**2/(64*Uhor**3)
d2Vh_dVLos2_2=A/(4*Uhor)-(2*A*VLos2_0+B*VLos1_0)**2/(64*Uhor**3)
d2Vh_dVLos1dVLos2=B/(8*Uhor)-(2*A*VLos1_0+B*VLos2_0)*(2*A*VLos2_0+B*VLos1_0)/(64*Uhor**3)

#Taylor decomposition
DeltaWS_lin=dVh_dVLos1*(VLos1-VLos1_0)+dVh_dVLos2*(VLos2-VLos2_0)
DeltaWS_sq=d2Vh_dVLos1_2*(VLos1-VLos1_0)**2+d2Vh_dVLos2_2*(VLos2-VLos2_0)**2+2*d2Vh_dVLos1dVLos2*(VLos1-VLos1_0)*(VLos2-VLos2_0)

#%% Plots
plt.close('all')
plt.figure(figsize=(16,10))
ax=plt.subplot(2,1,1)
plt.plot(DeltaZBeamsBins,DeltaWS/Uhor*100)
plt.gca().set_prop_cycle(None)
plt.plot(DeltaZBeamsBins,DeltaWS_lin/Uhor*100,'--')
plt.plot(0,0,'-k',label='Laura\'s equations')
plt.plot(0,0,'--k',label='1st order model')
plt.legend()
plt.grid()
plt.ylim([-2,12])

ax=plt.subplot(2,1,2)
plt.plot(DeltaZBeamsBins,DeltaWS/Uhor*100)
plt.gca().set_prop_cycle(None)
plt.plot(DeltaZBeamsBins,(DeltaWS_lin+DeltaWS_sq/2)/Uhor*100,'--')
plt.plot(0,0,'-k',label='Laura\'s equations')
plt.plot(0,0,'--k',label='2nd order model')
plt.legend()
plt.grid()
plt.ylim([-2,12])

