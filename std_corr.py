#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This script produces a figure depicting the first standard deviation of 
# correlation coefficient variance over the 10 calibration windows for both
# surface air temperature and precipitation.
"""
Created on Mon Feb 10 16:31:38 2020

@author: huiskamp
"""

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap 
import numpy as np
from netCDF4 import Dataset as CDF
import copy as cp


SAM_corrs = CDF('SAM_corrs.nc','r')
std = CDF('corr_std.nc','r')

sam_sat_corr = SAM_corrs.variables['SAM_SAT_land_model']
sam_pre_corr = SAM_corrs.variables['SAM_pre_land_model']
std_sat = std.variables['std_sat_wdw']
std_pre = std.variables['std_precip_wdw']

lat = std.variables['latitude']; lon = std.variables['longitude']

xgrid,ygrid = np.meshgrid(lon,lat)

minval = 0; maxval = 0.2
rng = 21
v = np.round(np.linspace(minval, maxval, rng),2)

# Set up color pallette
palette = cp.copy(plt.cm.YlOrRd)
palette.set_over('r') # extreme values will be red

# Set up the map
fig = plt.figure(figsize=(12,6))
 
ax1 = plt.subplot2grid((6, 3), (0, 0), colspan=9)
ax2 = plt.subplot2grid((6, 3), (1, 0), colspan=9)
ax3 = plt.subplot2grid((6, 3), (2, 0), colspan=9)
ax4 = plt.subplot2grid((6, 3), (3, 0), colspan=9)
ax5 = plt.subplot2grid((6, 3), (4, 0), colspan=9)
ax6 = plt.subplot2grid((6, 3), (5, 0), colspan=9)
            
m1=Basemap(llcrnrlon=0, llcrnrlat=-90,urcrnrlon=360,urcrnrlat=0, projection='cyl',resolution='c', ax=ax1 )
m2=Basemap(llcrnrlon=0, llcrnrlat=-90,urcrnrlon=360,urcrnrlat=0, projection='cyl',resolution='c', ax=ax2 )
m3=Basemap(llcrnrlon=0, llcrnrlat=-90,urcrnrlon=360,urcrnrlat=0, projection='cyl',resolution='c', ax=ax3 )
m4=Basemap(llcrnrlon=0, llcrnrlat=-90,urcrnrlon=360,urcrnrlat=0, projection='cyl',resolution='c', ax=ax4 )
m5=Basemap(llcrnrlon=0, llcrnrlat=-90,urcrnrlon=360,urcrnrlat=0, projection='cyl',resolution='c', ax=ax5 )
m6=Basemap(llcrnrlon=0, llcrnrlat=-90,urcrnrlon=360,urcrnrlat=0, projection='cyl',resolution='c', ax=ax6 )


m1.drawcoastlines(linewidth=0.2,zorder=3); m1.fillcontinents(color='grey',zorder=1) 
m2.drawcoastlines(linewidth=0.2,zorder=3); m2.fillcontinents(color='grey',zorder=1)
m3.drawcoastlines(linewidth=0.2,zorder=3); m3.fillcontinents(color='grey',zorder=1) 
m4.drawcoastlines(linewidth=0.2,zorder=3); m4.fillcontinents(color='grey',zorder=1)
m5.drawcoastlines(linewidth=0.2,zorder=3); m5.fillcontinents(color='grey',zorder=1)
m6.drawcoastlines(linewidth=0.2,zorder=3); m6.fillcontinents(color='grey',zorder=1)
lons, lats = m1(xgrid, ygrid)
m1.drawparallels(np.arange(-90.,91.,15.),labels=[1,0,0,1],linewidth=0.3,fontsize=5)
m1.drawmeridians(np.arange(-180.,181.,30.),labels=[False,False,False,False],linewidth=0.3,fontsize=5)
m2.drawparallels(np.arange(-90.,91.,15.),labels=[1,0,0,1],linewidth=0.3,fontsize=5)
m2.drawmeridians(np.arange(-180.,181.,30.),labels=[False,False,False,False],linewidth=0.3,fontsize=5)
m3.drawparallels(np.arange(-90.,91.,15.),labels=[1,0,0,1],linewidth=0.3,fontsize=5)
m3.drawmeridians(np.arange(-180.,181.,30.),labels=[False,False,False,False],linewidth=0.3,fontsize=5)
m4.drawparallels(np.arange(-90.,91.,15.),labels=[1,0,0,1],linewidth=0.3,fontsize=5)
m4.drawmeridians(np.arange(-180.,181.,30.),labels=[False,False,False,False],linewidth=0.3,fontsize=5)
m5.drawparallels(np.arange(-90.,91.,15.),labels=[1,0,0,1],linewidth=0.3,fontsize=5)
m5.drawmeridians(np.arange(-180.,181.,30.),labels=[False,False,False,False],linewidth=0.3,fontsize=5)
m6.drawparallels(np.arange(-90.,91.,15.),labels=[1,0,0,1],linewidth=0.3,fontsize=5)
m6.drawmeridians(np.arange(-180.,181.,30.),labels=[1,0,0,1],linewidth=0.3,fontsize=5)
    
# Plot maps and hatching of significant regions
# Panel 1 - SAM-SAT (31yr wdw)
pl1 = m1.pcolor(lons[0:45,:],lats[0:45,:],std_sat[0,:,:].T,vmin=0,vmax=0.2,zorder=2,cmap=palette) # .T transposes the matrix (must be numpy array)
pl1s = m1.contour(lons[0:45,:],lats[0:45,:],sam_sat_corr[:,:].T,
	11, # number of levels
	linewidths=0.3,
	colors = 'k')
ax1.clabel(pl1s, fontsize=2, inline=1)
# Panel 2 - SAM-SAT (61yr wdw)
pl2 = m2.pcolor(lons[0:45,:],lats[0:45,:],std_sat[1,:,:].T,vmin=0,vmax=0.2,zorder=2,cmap=palette) 
pl2s = m2.contour(lons[0:45,:],lats[0:45,:],sam_sat_corr[:,:].T,
	11, # number of levels
	linewidths=0.3,
	colors = 'k')
ax2.clabel(pl2s, fontsize=2, inline=1)
# Panel 3 - SAM-SAT (91yr wdw)
pl3 = m3.pcolor(lons[0:45,:],lats[0:45,:],std_sat[2,:,:].T,vmin=0,vmax=0.2,zorder=2,cmap=palette)
pl3s = m3.contour(lons[0:45,:],lats[0:45,:],sam_sat_corr[:,:].T,
	11, # number of levels
	linewidths=0.3,
	colors = 'k')
ax3.clabel(pl3s, fontsize=2, inline=1)
# Panel 4 - SAM-Precip (31yr wdw)
pl4 = m4.pcolor(lons[0:45,:],lats[0:45,:],std_pre[0,:,:].T,vmin=0,vmax=0.2,zorder=2,cmap=palette)
pl4s = m4.contour(lons[0:45,:],lats[0:45,:],sam_pre_corr[:,:].T,
	11, # number of levels
	linewidths=0.3,
	colors = 'k')
ax4.clabel(pl4s, fontsize=2, inline=1)
# Panel 5 - SAM-Precip (61yr wdw)
pl5 = m5.pcolor(lons[0:45,:],lats[0:45,:],std_pre[1,:,:].T,vmin=0,vmax=0.2,zorder=2,cmap=palette) 
pl5s = m5.contour(lons[0:45,:],lats[0:45,:],sam_pre_corr[:,:].T,
	11, # number of levels
	linewidths=0.3,
	colors = 'k')
ax5.clabel(pl5s, fontsize=2, inline=1)
# Panel 6 - SAM-precip (91yr wdw)
pl6 = m6.pcolor(lons[0:45,:],lats[0:45,:],std_pre[2,:,:].T,vmin=0,vmax=0.2,zorder=2,cmap=palette)
pl6s = m6.contour(lons[0:45,:],lats[0:45,:],sam_pre_corr[:,:].T,
	11, # number of levels
	linewidths=0.3,
	colors = 'k')
ax6.clabel(pl6s, fontsize=2, inline=1)


# colourbar
cax = plt.axes([0.75, 0.11, 0.019, 0.76]) 
cbarl=plt.colorbar(pl1,cax=cax,orientation='vertical',extend='max',shrink=0.5)
cbarl.set_ticks(v)
cbarl.ax.tick_params(labelsize=10)

# Include a) b) c) labels....
plt.text(-0.04, 2.16, 'a)', fontsize=6, transform=ax2.transAxes, horizontalalignment='left')
plt.text(-0.04, 0.98, 'b)', fontsize=6, transform=ax2.transAxes, horizontalalignment='left')
plt.text(-0.04, -0.24, 'c)', fontsize=6, transform=ax2.transAxes, horizontalalignment='left')
plt.text(-0.04, -1.46, 'd)', fontsize=6, transform=ax2.transAxes, horizontalalignment='left')
plt.text(-0.04, -2.64, 'e)', fontsize=6, transform=ax2.transAxes, horizontalalignment='left')
plt.text(-0.04, -3.82, 'f)', fontsize=6, transform=ax2.transAxes, horizontalalignment='left')

plt.savefig("figures/std_corr_python.pdf",bbox_inches='tight')