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

model_string = str(sys.argv[1])         # options: ref, weak, strong, core
analysis_depth = float(sys.argv[2])     # m (depth for DP extraction and central point of shear stress derivative)
analysis_depth_dz = float(sys.argv[3])  # m (depth interval for shear stress derivative)
ds = float(sys.argv[4])                 # m (distance from slab to pull out DP)
dz = float(sys.argv[5])                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)

if model_string == "ref":
    model_bothfree = "2D_compositional_subd_lower-res_new"
    model_fixedSP  = "2D_compositional_subd_FixedSP_lower-res_new"
    model_fixedOP  = "2D_compositional_subd_FixedOP_lower-res_new"
elif model_string == "weak":
    model_bothfree="2D_compositional_subd_lower-res_new_WeakPlates"
    model_fixedSP="2D_compositional_subd_lower-res_new_FixedSP_WeakPlates"
    model_fixedOP="2D_compositional_subd_lower-res_new_FixedOP_WeakPlates"
elif model_string == "strong":
    model_bothfree="2D_compositional_subd_lower-res_new_StiffPlates"
    model_fixedSP="2D_compositional_subd_lower-res_new_FixedSP_StiffPlates"
    model_fixedOP="2D_compositional_subd_lower-res_new_FixedOP_StiffPlates"
else:
    exit()

text_bothfree = ''.join(['text_files/',model_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_fixedSP  = ''.join(['text_files/',model_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_fixedOP  = ''.join(['text_files/',model_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])

bothfree = np.loadtxt((text_bothfree)) # format: # time [Myr], DP_mod_shall, DP_mod_deep, DP_mod, DP_anal, dip, stress term, stress term B, stress term C
fixedSP  = np.loadtxt((text_fixedSP))
fixedOP  = np.loadtxt((text_fixedOP))
 
plot_name_png = ''.join(['plots/DP-comparisons/compilations/',model_string,'-mods_z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.w-curvature.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/',model_string,'-mods_z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.w-curvature.pdf'])

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
curvature_thresh = 0.05
# free plates
for i in range(3,len(bothfree)):
    if np.abs(bothfree[i,12]) < curvature_thresh: 
        plt.scatter(bothfree[i,4]/1.e6,bothfree[i,3]/1.e6,color='slategray',s=8,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
        plt.scatter(bothfree[i,4]/1.e6,(bothfree[i,3]-bothfree[i,6])/1.e6,color='slategray',s=10,edgecolor='black',label='both free',lw=0.25,zorder=4) # corrected
        plt.plot([bothfree[i,4]/1.e6, bothfree[i,4]/1.e6], [bothfree[i,3]/1.e6, (bothfree[i,3]-bothfree[i,6])/1.e6], color='slategray', linewidth=0.7, zorder=2, alpha=0.25)
# fixed SP
for i in range(3,len(fixedSP)):
    if np.abs(fixedSP[i,12]) < curvature_thresh:     
        plt.scatter(fixedSP[i,4]/1.e6,fixedSP[i,3]/1.e6,color='slateblue',s=8,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
        plt.scatter(fixedSP[i,4]/1.e6,(fixedSP[i,3]-fixedSP[i,6])/1.e6,color='slateblue',s=10,edgecolor='black',label='fixed SP',lw=0.25,zorder=4) # corrected
        plt.plot([fixedSP[i,4]/1.e6, fixedSP[i,4]/1.e6], [fixedSP[i,3]/1.e6, (fixedSP[i,3]-fixedSP[i,6])/1.e6], color='slateblue', linewidth=0.7, zorder=2, alpha=0.25)
# fixed OP
for i in range(3,len(fixedOP)):
    if np.abs(fixedOP[i,12]) < curvature_thresh: 
        plt.scatter(fixedOP[i,4]/1.e6,fixedOP[i,3]/1.e6,color='darkseagreen',s=8,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
        plt.scatter(fixedOP[i,4]/1.e6,(fixedOP[i,3]-fixedOP[i,6])/1.e6,color='darkseagreen',s=10,edgecolor='black',label='fixed OP',lw=0.25,zorder=4) # corrected
        plt.plot([fixedOP[i,4]/1.e6, fixedOP[i,4]/1.e6], [fixedOP[i,3]/1.e6, (fixedOP[i,3]-fixedOP[i,6])/1.e6], color='darkseagreen', linewidth=0.7, zorder=2, alpha=0.25)

# axis stuff
plt.xlim(-5,  30); plt.ylim(-5,  30)
plt.plot([-5, 30], [-5, 30], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
# plt.xlabel("slab buoyancy  [MPa]",size=7.5)
plt.ylabel("extracted stress  [MPa]",size=7.5)
# plt.legend(fontsize=5.5,loc='best',facecolor='white',fancybox=True, framealpha=1, frameon=False)
ax.annotate('low curvature', xy=(0.0125,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=9,color='k') 
fixed_aspect_ratio(1)

anal_curve = np.zeros((100,2))
for j in range(len(anal_curve)):
	anal_curve[j,0] = j # dip [deg]
	drho = 50. # kg/m3
	slabnorm_thick = 80.e3 # m  -- CHECK!!
	anal_curve[j,1] = DP_anal = (drho * 9.81 * slabnorm_thick * np.cos(np.deg2rad(j)))/1.e6 # DP  [MPa]

ax2=fig.add_subplot(gs[0,1])
# free plates
for i in range(3,len(bothfree)):
    if np.abs(bothfree[i,12]) < curvature_thresh: 
        plt.scatter(bothfree[i,5],bothfree[i,3]/1.e6,color='slategray',s=8,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
        plt.scatter(bothfree[i,5],(bothfree[i,3]-bothfree[i,6])/1.e6,color='slategray',s=10,edgecolor='black',label='both free',lw=0.25,zorder=4) # corrected
        plt.plot([bothfree[i,5], bothfree[i,5]], [bothfree[i,3]/1.e6, (bothfree[i,3]-bothfree[i,6])/1.e6], color='slategray', linewidth=0.7, zorder=2, alpha=0.25)
# fixed SP
for i in range(3,len(fixedSP)):
    if np.abs(fixedSP[i,12]) < curvature_thresh: 
        plt.scatter(fixedSP[i,5],fixedSP[i,3]/1.e6,color='slateblue',s=8,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
        plt.scatter(fixedSP[i,5],(fixedSP[i,3]-fixedSP[i,6])/1.e6,color='slateblue',s=10,edgecolor='black',label='fixed SP',lw=0.25,zorder=4) # corrected
        plt.plot([fixedSP[i,5], fixedSP[i,5]], [fixedSP[i,3]/1.e6, (fixedSP[i,3]-fixedSP[i,6])/1.e6], color='slateblue', linewidth=0.7, zorder=2, alpha=0.25)
# fixed SP
for i in range(3,len(fixedOP)):
    if np.abs(fixedOP[i,12]) < curvature_thresh: 
        plt.scatter(fixedOP[i,5],fixedOP[i,3]/1.e6,color='darkseagreen',s=8,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
        plt.scatter(fixedOP[i,5],(fixedOP[i,3]-fixedOP[i,6])/1.e6,color='darkseagreen',s=10,edgecolor='black',label='fixed OP',lw=0.25,zorder=4) # corrected
        plt.plot([fixedOP[i,5], fixedOP[i,5]], [fixedOP[i,3]/1.e6, (fixedOP[i,3]-fixedOP[i,6])/1.e6], color='darkseagreen', linewidth=0.7, zorder=2, alpha=0.25)
# axis stuff
plt.xlim(45,  100); plt.ylim(-5,  30)
# analytical curve
plt.plot(anal_curve[:,0], anal_curve[:,1], color='black', linewidth=1, zorder=1)

# plt.plot([-20, 40], [-20, 40], color='bisque', linewidth=2, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=90, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax2.tick_params(axis='x', labelsize=6)
ax2.tick_params(axis='y', labelsize=6)
# plt.xlabel('$\Theta$  [$^\circ$]',size=7.5)
plt.ylabel("extracted stress  [MPa]",size=7.5)
fixed_aspect_ratio(1)


# plot viscosity field
ax3=fig.add_subplot(gs[1,0])
# free plates
for i in range(3,len(bothfree)):
    if np.abs(bothfree[i,12]) > curvature_thresh: 
        plt.scatter(bothfree[i,4]/1.e6,bothfree[i,3]/1.e6,color='slategray',s=8,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
        plt.scatter(bothfree[i,4]/1.e6,(bothfree[i,3]-bothfree[i,6])/1.e6,color='slategray',s=10,edgecolor='black',label='both free',lw=0.25,zorder=4) # corrected
        plt.plot([bothfree[i,4]/1.e6, bothfree[i,4]/1.e6], [bothfree[i,3]/1.e6, (bothfree[i,3]-bothfree[i,6])/1.e6], color='slategray', linewidth=0.7, zorder=2, alpha=0.25)
# fixed SP
for i in range(3,len(fixedSP)):
    if np.abs(fixedSP[i,12]) > curvature_thresh:     
        plt.scatter(fixedSP[i,4]/1.e6,fixedSP[i,3]/1.e6,color='slateblue',s=8,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
        plt.scatter(fixedSP[i,4]/1.e6,(fixedSP[i,3]-fixedSP[i,6])/1.e6,color='slateblue',s=10,edgecolor='black',label='fixed SP',lw=0.25,zorder=4) # corrected
        plt.plot([fixedSP[i,4]/1.e6, fixedSP[i,4]/1.e6], [fixedSP[i,3]/1.e6, (fixedSP[i,3]-fixedSP[i,6])/1.e6], color='slateblue', linewidth=0.7, zorder=2, alpha=0.25)
# fixed OP
for i in range(3,len(fixedOP)):
    if np.abs(fixedOP[i,12]) > curvature_thresh: 
        plt.scatter(fixedOP[i,4]/1.e6,fixedOP[i,3]/1.e6,color='darkseagreen',s=8,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
        plt.scatter(fixedOP[i,4]/1.e6,(fixedOP[i,3]-fixedOP[i,6])/1.e6,color='darkseagreen',s=10,edgecolor='black',label='fixed OP',lw=0.25,zorder=4) # corrected
        plt.plot([fixedOP[i,4]/1.e6, fixedOP[i,4]/1.e6], [fixedOP[i,3]/1.e6, (fixedOP[i,3]-fixedOP[i,6])/1.e6], color='darkseagreen', linewidth=0.7, zorder=2, alpha=0.25)

# axis stuff
plt.xlim(-5,  30); plt.ylim(-5,  30)
plt.plot([-5, 30], [-5, 30], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax3.tick_params(axis='x', labelsize=6)
ax3.tick_params(axis='y', labelsize=6)
plt.xlabel("slab buoyancy  [MPa]",size=7.5)
plt.ylabel("extracted stress  [MPa]",size=7.5)
# plt.legend(fontsize=5.5,loc='best',facecolor='white',fancybox=True, framealpha=1, frameon=False)
ax3.annotate('high curvature', xy=(0.0125,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=9,color='k') 
fixed_aspect_ratio(1)

anal_curve = np.zeros((100,2))
for j in range(len(anal_curve)):
    anal_curve[j,0] = j # dip [deg]
    drho = 50. # kg/m3
    slabnorm_thick = 80.e3 # m  -- CHECK!!
    anal_curve[j,1] = DP_anal = (drho * 9.81 * slabnorm_thick * np.cos(np.deg2rad(j)))/1.e6 # DP  [MPa]

ax4=fig.add_subplot(gs[1,1])
# free plates
for i in range(3,len(bothfree)):
    if np.abs(bothfree[i,12]) > curvature_thresh: 
        plt.scatter(bothfree[i,5],bothfree[i,3]/1.e6,color='slategray',s=8,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
        plt.scatter(bothfree[i,5],(bothfree[i,3]-bothfree[i,6])/1.e6,color='slategray',s=10,edgecolor='black',label='both free',lw=0.25,zorder=4) # corrected
        plt.plot([bothfree[i,5], bothfree[i,5]], [bothfree[i,3]/1.e6, (bothfree[i,3]-bothfree[i,6])/1.e6], color='slategray', linewidth=0.7, zorder=2, alpha=0.25)
# fixed SP
for i in range(3,len(fixedSP)):
    if np.abs(fixedSP[i,12]) > curvature_thresh: 
        plt.scatter(fixedSP[i,5],fixedSP[i,3]/1.e6,color='slateblue',s=8,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
        plt.scatter(fixedSP[i,5],(fixedSP[i,3]-fixedSP[i,6])/1.e6,color='slateblue',s=10,edgecolor='black',label='fixed SP',lw=0.25,zorder=4) # corrected
        plt.plot([fixedSP[i,5], fixedSP[i,5]], [fixedSP[i,3]/1.e6, (fixedSP[i,3]-fixedSP[i,6])/1.e6], color='slateblue', linewidth=0.7, zorder=2, alpha=0.25)
# fixed SP
for i in range(3,len(fixedOP)):
    if np.abs(fixedOP[i,12]) > curvature_thresh: 
        plt.scatter(fixedOP[i,5],fixedOP[i,3]/1.e6,color='darkseagreen',s=8,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
        plt.scatter(fixedOP[i,5],(fixedOP[i,3]-fixedOP[i,6])/1.e6,color='darkseagreen',s=10,edgecolor='black',label='fixed OP',lw=0.25,zorder=4) # corrected
        plt.plot([fixedOP[i,5], fixedOP[i,5]], [fixedOP[i,3]/1.e6, (fixedOP[i,3]-fixedOP[i,6])/1.e6], color='darkseagreen', linewidth=0.7, zorder=2, alpha=0.25)
# axis stuff
plt.xlim(45,  100); plt.ylim(-5,  30)
# analytical curve
plt.plot(anal_curve[:,0], anal_curve[:,1], color='black', linewidth=1, zorder=1)

# plt.plot([-20, 40], [-20, 40], color='bisque', linewidth=2, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=90, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax4.tick_params(axis='x', labelsize=6)
ax4.tick_params(axis='y', labelsize=6)
plt.xlabel('$\Theta$  [$^\circ$]',size=7.5)
plt.ylabel("extracted stress  [MPa]",size=7.5)
fixed_aspect_ratio(1)


plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)
# plt.savefig(plot_name_pdf, bbox_inches='tight', format='pdf')


