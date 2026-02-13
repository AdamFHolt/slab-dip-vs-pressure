#!/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.gridspec import GridSpec
import sys, os, subprocess
from scipy.signal import savgol_filter
from scipy.interpolate import splrep, splev

analysis_depth = float(sys.argv[1])             # m (depth for DP extraction and central point of shear stress derivative)
analysis_depth_dz = float(sys.argv[2])  # m (depth interval for shear stress derivative)
ds = float(sys.argv[3])                                 # m (distance from slab to pull out DP)
dz = float(sys.argv[4])                                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)

model_bothfree = "2D_compositional_subd_lower-res_new"
model_fixedSP  = "2D_compositional_subd_FixedSP_lower-res_new"
model_fixedOP  = "2D_compositional_subd_FixedOP_lower-res_new"

text_bothfree = ''.join(['text_files/',model_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_fixedSP  = ''.join(['text_files/',model_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_fixedOP  = ''.join(['text_files/',model_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])

bothfree = np.loadtxt((text_bothfree)) # format: # time [Myr], DP_mod_shall, DP_mod_deep, DP_mod, DP_anal, dip, stress term, stress term B, stress term C
fixedSP  = np.loadtxt((text_fixedSP))
fixedOP  = np.loadtxt((text_fixedOP))
 
plot_name_png = ''.join(['plots/DP-comparisons/z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.png'])

fig=plt.figure()
gs=GridSpec(2,2) 

#### SCATTER PLOTS ####
def fixed_aspect_ratio(ratio):
    '''
    Set a fixed aspect ratio on matplotlib plots 
    regardless of axis units
    '''
    xvals,yvals = plt.gca().axes.get_xlim(),plt.gca().axes.get_ylim()

    xrange = xvals[1]-xvals[0]
    yrange = yvals[1]-yvals[0]
    plt.gca().set_aspect(ratio*(xrange/yrange), adjustable='box')


# plot viscosity field
ax=fig.add_subplot(gs[0,0])
# free plates
plt.scatter(bothfree[1:,4]/1.e6,bothfree[1:,3]/1.e6,color='slategray',s=10,edgecolor='black',lw=0,zorder=3,alpha=0.5) # uncorrected
plt.scatter(bothfree[1:,4]/1.e6,(bothfree[1:,3]-bothfree[1:,6])/1.e6,color='slategray',s=10,edgecolor='black',label='both free',lw=0.25,zorder=4) # corrected
for i in range(1,len(bothfree)):
   plt.plot([bothfree[i,4]/1.e6, bothfree[i,4]/1.e6], [bothfree[i,3]/1.e6, (bothfree[i,3]-bothfree[i,6])/1.e6], color='slategray', linewidth=1, zorder=2, alpha=0.5)
# fixed SP
plt.scatter(fixedSP[1:,4]/1.e6,fixedSP[1:,3]/1.e6,color='slateblue',s=10,edgecolor='black',lw=0,zorder=3,alpha=0.5) # uncorrected
plt.scatter(fixedSP[1:,4]/1.e6,(fixedSP[1:,3]-fixedSP[1:,6])/1.e6,color='slateblue',s=10,edgecolor='black',label='fixed SP',lw=0.25,zorder=4) # corrected
for i in range(1,len(fixedSP)):
    plt.plot([fixedSP[i,4]/1.e6, fixedSP[i,4]/1.e6], [fixedSP[i,3]/1.e6, (fixedSP[i,3]-fixedSP[i,6])/1.e6], color='slateblue', linewidth=1, zorder=2, alpha=0.5)
# fixed SP
plt.scatter(fixedOP[1:,4]/1.e6,fixedOP[1:,3]/1.e6,color='darkseagreen',s=10,edgecolor='black',lw=0,zorder=3,alpha=0.5) # uncorrected
plt.scatter(fixedOP[1:,4]/1.e6,(fixedOP[1:,3]-fixedOP[1:,6])/1.e6,color='darkseagreen',s=10,edgecolor='black',label='fixed OP',lw=0.25,zorder=4) # corrected
for i in range(1,len(fixedOP)):
    plt.plot([fixedOP[i,4]/1.e6, fixedOP[i,4]/1.e6], [fixedOP[i,3]/1.e6, (fixedOP[i,3]-fixedOP[i,6])/1.e6], color='darkseagreen', linewidth=1, zorder=2, alpha=0.5)
# axis stuff
plt.xlim(-20,  40); plt.ylim(-20,  40)
plt.plot([-20, 40], [-20, 40], color='bisque', linewidth=2, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel(r"$\Delta$P$\mathregular{_{analytical}}$  [MPa]",size=7.5)
plt.ylabel(r"$\Delta$P$\mathregular{_{model}}$  [MPa]",size=7.5)
plt.legend(fontsize=5.5,loc='lower right',facecolor='white',fancybox=True, framealpha=1, frameon=False)

fixed_aspect_ratio(1)
# coeff = np.corrcoef(dips_lall[:,2],slab2_dips_LallLocations[:,2])[1,0]
# coeff_string = ''.join(['$\mathregular{R_{Pearson}}$ = ',str(round(coeff, 3)),' (n = ',str(len(dips_lall)),')'])
# ax.text(0.05,0.91,coeff_string,size=7, color="black",transform = ax.transAxes)

plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)


