#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 16:31:38 2020

@author: huiskamp
"""

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap 
from scipy.io import loadmat
import numpy as np

data = loadmat('DataFiles/ENSO_corrs.mat')
lat = data['lat']; lon = data['lon']
xgrid,ygrid = np.meshgrid(lon,lat)
n34_SAM_SAT = data['n34_SAM_SAT']; n34_SAM_pre = data['n34_SAM_pre']
SAT_corr_sig = data['SAT_corr_sig']; pre_corr_sig = data['pre_corr_sig']

minval = -1; maxval = 1
rng = 21
v = np.round(np.linspace(minval, maxval, rng),2)

for dat in range(2):
    if dat == 0:
        d = n34_SAM_SAT
        d_sig = SAT_corr_sig
    elif dat == 1:
        d = n34_SAM_pre
        d_sig = pre_corr_sig
    # Set up the map
    fig = plt.figure(figsize=(12,6))
    ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=9)
    ax2 = plt.subplot2grid((3, 3), (1, 0), colspan=9)
    ax3 = plt.subplot2grid((3, 3), (2, 0), colspan=9)
            
    m1=Basemap(llcrnrlon=0, llcrnrlat=-90,urcrnrlon=360,urcrnrlat=0, projection='cyl',resolution='c', ax=ax1 )
    m2=Basemap(llcrnrlon=0, llcrnrlat=-90,urcrnrlon=360,urcrnrlat=0, projection='cyl',resolution='c', ax=ax2 )
    m3=Basemap(llcrnrlon=0, llcrnrlat=-90,urcrnrlon=360,urcrnrlat=0, projection='cyl',resolution='c', ax=ax3 )
    m1.drawcoastlines(linewidth=1.5); m2.drawcoastlines(linewidth=1.5); m3.drawcoastlines(linewidth=1.5)
    lons, lats = m1(xgrid, ygrid)
    m1.drawparallels(np.arange(-90.,91.,15.),labels=[1,0,0,1],linewidth=0.5,fontsize=10)
    m1.drawmeridians(np.arange(-180.,181.,30.),labels=[1,0,0,1],linewidth=0.5,fontsize=10)
    m2.drawparallels(np.arange(-90.,91.,15.),labels=[1,0,0,1],linewidth=0.5,fontsize=10)
    m2.drawmeridians(np.arange(-180.,181.,30.),labels=[1,0,0,1],linewidth=0.5,fontsize=10)
    m3.drawparallels(np.arange(-90.,91.,15.),labels=[1,0,0,1],linewidth=0.5,fontsize=10)
    m3.drawmeridians(np.arange(-180.,181.,30.),labels=[1,0,0,1],linewidth=0.5,fontsize=10)
    
    # Plot maps and hatching of significant regions
    pl1 = m1.contourf(lons[0:45,:],lats[0:45,:],d[0,:,:],v,cmap='bwr',extend='both')
    pl1s = m1.pcolor(lons[0:45,:],lats[0:45,:],d_sig[0,:,:],hatch='///',alpha=0)
    pl2 = m2.contourf(lons[0:45,:],lats[0:45,:],d[1,:,:],v,cmap='bwr',extend='both')
    pl2s = m2.pcolor(lons[0:45,:],lats[0:45,:],d_sig[1,:,:],hatch='///',alpha=0)
    pl3 = m3.contourf(lons[0:45,:],lats[0:45,:],d[2,:,:],v,cmap='bwr',extend='both')
    pl3s = m3.pcolor(lons[0:45,:],lats[0:45,:],d_sig[2,:,:],hatch='///',alpha=0)
    
    # colourbar
    cax = plt.axes([0.75, 0.11, 0.019, 0.76]) 
    cbarl=plt.colorbar(pl1,cax=cax,orientation='vertical')
    cbarl.set_ticks(v)
    cbarl.ax.tick_params(labelsize=10)
    
    # Include a) b) c) labels....
    plt.text(-0.04, 2.15, 'a)', fontsize=10, transform=ax2.transAxes, horizontalalignment='left')
    plt.text(-0.04, 0.97, 'b)', fontsize=10, transform=ax2.transAxes, horizontalalignment='left')
    plt.text(-0.04, -0.25, 'c)', fontsize=10, transform=ax2.transAxes, horizontalalignment='left')
    
    if dat == 0:
        plt.savefig("figures/n34_SAM_SAT_corr.pdf",bbox_inches='tight')
        plt.savefig("figures/n34_SAM_SAT_corr.png",bbox_inches='tight',dpi=300)
    elif dat == 1:
        plt.savefig("figures/n34_SAM_pre_corr.pdf",bbox_inches='tight')
        plt.savefig("figures/n34_SAM_pre_corr.png",bbox_inches='tight',dpi=300)
    
    
    
    
    




















