# -*- coding: utf-8 -*-
"""
NML-x error due to terrain slope
"""
import os
cd=os.path.dirname(__file__)
from matplotlib import pyplot as plt
import numpy as np

#%% Inputs
alpha=30 #[deg] lidar half-opening angle
Uhor=10 #[m/s]  #horizontal wind speed at hub height
theta=0 #[deg] #relative wind direction 
ShearBins=np.arange(-0.05,0.41,0.05) #shear exponents
DeltaZBeamsBins=np.arange(-40,41,2) #[m] height differences
DVLosBins=np.arange(-2,2.1,0.1)
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

#wind flow reconstruction
Vx_WFR=(VLos1+VLos2)/(2*np.cos(np.radians(alpha)))
Vy_WFR=(VLos1-VLos2)/(2*np.sin(np.radians(alpha)))

WS_WFR=(Vx_WFR**2+Vy_WFR**2)**0.5
WD_WFR=np.degrees(np.arctan(Vy_WFR/Vx_WFR))

#error
DeltaWS=WS_WFR-Uhor
DeltaWD=WD_WFR-theta

#error-free LOS velocities
VLos1_0=Uhor*np.cos(np.radians(alpha-theta))
VLos2_0=Uhor*np.cos(np.radians(alpha+theta))
VLOS1,VLOS2=np.meshgrid(DVLosBins+VLos1_0,DVLosBins+VLos2_0)

#compact WFR form
c=np.cos(np.radians(alpha))
s=np.sin(np.radians(alpha))
A=1/c**2+1/s**2
B=2*(1/c**2-1/s**2)
WS_WFR2=0.5*(A*(VLos1**2+VLos2**2)+B*VLos1*VLos2)**0.5
WD_WFR2=np.degrees(np.arctan((VLos1-VLos2)/(VLos1+VLos2)*c/s))

#first order derivatives
dWS_dVLos1=1/8*(2*A*VLos1_0+B*VLos2_0)/Uhor
dWS_dVLos2=1/8*(2*A*VLos2_0+B*VLos1_0)/Uhor

dWD_dVLos1=1/(((VLos1_0-VLos2_0)/(VLos1_0+VLos2_0)*c/s)**2+1)*2*VLos2_0/(VLos1_0+VLos2_0)**2*c/s
dWD_dVLos2=1/(((VLos1_0-VLos2_0)/(VLos1_0+VLos2_0)*c/s)**2+1)*(-2*VLos1_0)/(VLos1_0+VLos2_0)**2*c/s

#second order derivatives
d2WS_dVLos1_2=A/(4*Uhor)-(2*A*VLos1_0+B*VLos2_0)**2/(64*Uhor**3)
d2WS_dVLos2_2=A/(4*Uhor)-(2*A*VLos2_0+B*VLos1_0)**2/(64*Uhor**3)
d2WS_dVLos1dVLos2=B/(8*Uhor)-(2*A*VLos1_0+B*VLos2_0)*(2*A*VLos2_0+B*VLos1_0)/(64*Uhor**3)

#Taylor decomposition
DeltaWS_lin=dWS_dVLos1*(VLos1-VLos1_0)+dWS_dVLos2*(VLos2-VLos2_0)
DeltaWS_sq=d2WS_dVLos1_2*(VLos1-VLos1_0)**2+d2WS_dVLos2_2*(VLos2-VLos2_0)**2+2*d2WS_dVLos1dVLos2*(VLos1-VLos1_0)*(VLos2-VLos2_0)

DeltaWD_lin=np.degrees(dWD_dVLos1*(VLos1-VLos1_0)+dWD_dVLos2*(VLos2-VLos2_0))

#surface
DeltaWS_surf=0.5*(A*(VLOS1**2+VLOS2**2)+B*VLOS1*VLOS2)**0.5-Uhor
DeltaWD_surf=np.degrees(np.arctan((VLOS1-VLOS2)/(VLOS1+VLOS2)*c/s))-theta

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
plt.xlabel('Zlos2-Zlos1')
plt.ylabel(r'$\Delta$ Uhor/Uhor [%]')

ax=plt.subplot(2,1,2)
plt.plot(DeltaZBeamsBins,DeltaWS/Uhor*100)
plt.gca().set_prop_cycle(None)
plt.plot(DeltaZBeamsBins,(DeltaWS_lin+DeltaWS_sq/2)/Uhor*100,'--')
plt.plot(0,0,'-k',label='Laura\'s equations')
plt.plot(0,0,'--k',label='2nd order model')
plt.legend()
plt.grid()
plt.ylim([-2,12])
plt.xlabel('Zlos2-Zlos1')
plt.ylabel(r'$\Delta$ Uhor/Uhor [%]')

plt.figure(figsize=(16,5))
plt.plot(DeltaZBeamsBins,DeltaWD)
plt.gca().set_prop_cycle(None)
plt.plot(DeltaZBeamsBins,DeltaWD_lin,'--')
plt.plot(0,0,'-k',label='Laura\'s equations')
plt.plot(0,0,'--k',label='1st order model')
plt.legend()
plt.grid()
plt.ylim([-30,30])
plt.xlabel('Zlos2-Zlos1')
plt.ylabel(r'$\Delta$ theta [degrees]')

fig=plt.figure(figsize=(16,10))
ax=fig.add_subplot(111,projection='3d')
ax.plot_surface(VLOS1-VLos1_0,VLOS2-VLos2_0,DeltaWS_surf/Uhor*100,cmap='Grays',alpha=1)
ax.set_prop_cycle(None)
for i in range(len(ShearBins)):
    plt.plot(VLos1[:,i]-VLos1_0,VLos2[:,i]-VLos2_0,DeltaWS[:,i]/Uhor*100,'.',markersize=3,zorder=10)
ax.set_xlabel(r'$\Delta$ $V_L$ [m s$^{-1}$]')
ax.set_ylabel(r'$\Delta$ $V_R$ [m s$^{-1}$]')
ax.set_zlabel(r'$\Delta$ Uhor/Uhor [%]')
ax.view_init(elev=30, azim=-112)

fig=plt.figure(figsize=(16,10))
ax=fig.add_subplot(111,projection='3d')
ax.plot_surface(VLOS1-VLos1_0,VLOS2-VLos2_0,DeltaWD_surf,cmap='Grays',alpha=1)
ax.set_prop_cycle(None)
for i in range(len(ShearBins)):
    plt.plot(VLos1[:,i]-VLos1_0,VLos2[:,i]-VLos2_0,DeltaWD[:,i],'.',markersize=3,zorder=10)
ax.set_xlabel(r'$\Delta$ $V_L$ [m s$^{-1}$]')
ax.set_ylabel(r'$\Delta$ $V_R$ [m s$^{-1}$]')
ax.set_zlabel(r'$\Delta$ theta [degrees]')
ax.view_init(elev=30, azim=-112)



