# -*- coding: utf-8 -*-
"""
NML-x error due to terrain slope for report
"""
import os
cd=os.path.dirname(__file__)
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import xarray as xr
import numpy as np
import matplotlib as mpl
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.size'] = 12
mpl.rcParams['savefig.dpi']=500

#%% Inputs
gamma=15 #[deg] lidar half-opening angle
V_H=10 #[m/s]  #horizontal wind speed at hub height
theta=0 #[deg] #relative wind direction 
phi=0
alpha=np.array([-0.05,0.05,0.1,0.2,0.3,0.4,1])
DZ1=np.arange(-40,41)
DZ2=np.arange(-40,41)
Z_H=100 #[m] hub height

#%% Functions
def cos(x):
    return np.cos(np.radians(x))
def sin(x):
    return np.sin(np.radians(x))

#%% Initialization
Data=xr.Dataset()
Data['DZ1']=xr.DataArray(DZ1,coords={'DZ1':DZ1})
Data['DZ2']=xr.DataArray(DZ2,coords={'DZ2':DZ2})
Data['alpha']=xr.DataArray(alpha,coords={'alpha':alpha})

#%% Main

#samples V_H a the two beams
V_r1=V_H*(((Data.DZ1+Z_H)/Z_H)**Data.alpha)*cos(gamma-theta)*cos(phi)
V_r2=V_H*(((Data.DZ2+Z_H)/Z_H)**Data.alpha)*cos(gamma+theta)*cos(phi)

#wind flow reconstruction
V_x_rec=(V_r1+V_r2)/(2*cos(phi)*cos(gamma))
V_y_rec=(V_r1-V_r2)/(2*cos(phi)*sin(gamma))
V_H_rec=(V_x_rec**2+V_y_rec**2)**0.5
theta_rec=np.degrees(np.arctan(V_y_rec/V_x_rec))

#%% Plots
fig=plt.figure(figsize=(18,8))
gs = gridspec.GridSpec(2, 4,
    width_ratios=[1, 1,1, 0.05], 
    wspace=0.3, hspace=0.3)

ctr=0
for i in range(2):
    for j in range(3):
        ax = fig.add_subplot(gs[i, j])
        pc=plt.pcolor(Data.DZ1,Data.DZ2,(V_H_rec.isel(alpha=ctr).T-V_H)/V_H*100,cmap='seismic',vmin=-10,vmax=10)
        plt.title(r'$\alpha='+str(alpha[ctr])+'$')
        if i==1:
            plt.xlabel(r'$\Delta Z_1$')
        else:
            ax.set_xticklabels([])
        if j==0:
            plt.ylabel(r'$\Delta Z_2$')
        else:
            ax.set_yticklabels([])
        ctr+=1
cax = fig.add_subplot(gs[:, 3])
fig.colorbar(pc, cax=cax,label=r'Bias in $V_H$ [$\%$]')

fig=plt.figure(figsize=(18,8))
gs = gridspec.GridSpec(2, 4,
    width_ratios=[1, 1,1, 0.05], 
    wspace=0.3, hspace=0.3)

ctr=0
for i in range(2):
    for j in range(3):
        ax = fig.add_subplot(gs[i, j])
        pc=plt.pcolor(Data.DZ1,Data.DZ2,theta_rec.isel(alpha=ctr).T-theta,cmap='seismic',vmin=-10,vmax=10)
        plt.title(r'$\alpha='+str(alpha[ctr])+'$')
        if i==1:
            plt.xlabel(r'$\Delta Z_1$')
        else:
            ax.set_xticklabels([])
        if j==0:
            plt.ylabel(r'$\Delta Z_2$')
        else:
            ax.set_yticklabels([])
        ctr+=1
cax = fig.add_subplot(gs[:, 3])
fig.colorbar(pc, cax=cax,label=r'Bias in $\theta$ [$^\circ$]')