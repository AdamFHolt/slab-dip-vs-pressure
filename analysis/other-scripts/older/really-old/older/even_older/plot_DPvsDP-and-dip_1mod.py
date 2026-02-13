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

model_name = str(sys.argv[1])
analysis_depth_dz = float(sys.argv[2])  				# m (depth interval for shear stress derivative)
ds = float(sys.argv[3])                                 # m (distance from slab to pull out DP)
dz = float(sys.argv[4])                                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)

anal_depth1 = 200.0e3
anal_depth2 = 300.0e3
anal_depth3 = 400.0e3
anal_depth4 = 500.0e3

text_model1 = ''.join(['text_files/lo-res/',model_name,'.z',str(anal_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_model2 = ''.join(['text_files/lo-res/',model_name,'.z',str(anal_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_model3 = ''.join(['text_files/lo-res/',model_name,'.z',str(anal_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_model4 = ''.join(['text_files/lo-res/',model_name,'.z',str(anal_depth4/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])

#bothfree = np.loadtxt((text_bothfree)) # format: # time [Myr], DP_mod_shall, DP_mod_deep, DP_mod, DP_anal, dip, stress term, stress term B, stress term C
model1  = np.loadtxt((text_model1))
model2  = np.loadtxt((text_model2))
model3  = np.loadtxt((text_model3))
model4  = np.loadtxt((text_model4))

plot_loc   = ''.join(['plots/DP-comparisons/lo-res/single-models/',model_name])
if not os.path.exists(plot_loc):
    os.mkdir(plot_loc)

plot_name_png = ''.join([plot_loc,'/shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.png'])

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

# 200 km
plt.scatter(model1[1:,4]/1.e6,model1[1:,3]/1.e6,color='wheat',s=10,lw=0.0,zorder=3,alpha=0.5) # uncorrected
plt.scatter(model1[1:,4]/1.e6,(model1[1:,3]-model1[1:,6])/1.e6,color='wheat',s=10,edgecolor='black',label='200 km',lw=0.25,zorder=4) # corrected
for i in range(1,len(model1)):
    plt.plot([model1[i,4]/1.e6, model1[i,4]/1.e6], [model1[i,3]/1.e6, (model1[i,3]-model1[i,6])/1.e6], color='wheat', linewidth=1, zorder=2, alpha=0.5)
# 300 km
plt.scatter(model2[1:,4]/1.e6,model2[1:,3]/1.e6,color='gold',s=10,lw=0.0,zorder=3,alpha=0.5) # uncorrected
plt.scatter(model2[1:,4]/1.e6,(model2[1:,3]-model2[1:,6])/1.e6,color='gold',s=10,edgecolor='black',label='300 km',lw=0.25,zorder=4) # corrected
for i in range(1,len(model2)):
    plt.plot([model2[i,4]/1.e6, model2[i,4]/1.e6], [model2[i,3]/1.e6, (model2[i,3]-model2[i,6])/1.e6], color='gold', linewidth=1, zorder=2, alpha=0.5)
# 400 km
plt.scatter(model3[1:,4]/1.e6,model3[1:,3]/1.e6,color='darkorange',s=10,lw=0.0,zorder=3,alpha=0.5) # uncorrected
plt.scatter(model3[1:,4]/1.e6,(model3[1:,3]-model3[1:,6])/1.e6,color='darkorange',s=10,edgecolor='black',label='400 km',lw=0.25,zorder=4) # corrected
for i in range(1,len(model3)):
    plt.plot([model3[i,4]/1.e6, model3[i,4]/1.e6], [model3[i,3]/1.e6, (model3[i,3]-model3[i,6])/1.e6], color='darkorange', linewidth=1, zorder=2, alpha=0.5)
# 500 km
plt.scatter(model4[1:,4]/1.e6,model4[1:,3]/1.e6,color='brown',s=10,lw=0.0,zorder=3,alpha=0.5) # uncorrected
plt.scatter(model4[1:,4]/1.e6,(model4[1:,3]-model4[1:,6])/1.e6,color='brown',s=10,edgecolor='black',label='500 km',lw=0.25,zorder=4) # corrected
for i in range(1,len(model4)):
    plt.plot([model4[i,4]/1.e6, model4[i,4]/1.e6], [model4[i,3]/1.e6, (model4[i,3]-model4[i,6])/1.e6], color='brown', linewidth=1, zorder=2, alpha=0.5)

# axis stuff
plt.xlim(-30,  50); plt.ylim(-30,  50)
plt.plot([-30, 50], [-30, 50], color='bisque', linewidth=2, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel(r"$\Delta$P$\mathregular{_{analytical}}$  [MPa]",size=7.5)
plt.ylabel(r"$\Delta$P$\mathregular{_{model}}$  [MPa]",size=7.5)
plt.legend(fontsize=5.5,loc='lower right',facecolor='white',fancybox=True, framealpha=1, frameon=False)
fixed_aspect_ratio(1)

anal_curve = np.zeros((140,2))
for j in range(len(anal_curve)):
    anal_curve[j,0] = j # dip [deg]
    drho = 50. # kg/m3
    slabnorm_thick = 80.e3 # m  -- CHECK!!
    anal_curve[j,1] = DP_anal = (drho * 9.81 * slabnorm_thick * np.cos(np.deg2rad(j)))/1.e6 # DP  [MPa]



ax2=fig.add_subplot(gs[0,1])
# 200 km 
plt.scatter(model1[1:,5],model1[1:,3]/1.e6,color='wheat',s=8,edgecolor='black',lw=0,zorder=3,alpha=0.5) # uncorrected
plt.scatter(model1[1:,5],(model1[1:,3]-model1[1:,6])/1.e6,color='wheat',s=10,edgecolor='black',label='both free',lw=0.25,zorder=4) # corrected
for i in range(1,len(model1)):
   plt.plot([model1[i,5], model1[i,5]], [model1[i,3]/1.e6, (model1[i,3]-model1[i,6])/1.e6], color='wheat', linewidth=0.5, zorder=2, alpha=0.5)
# 300 km 
plt.scatter(model2[1:,5],model2[1:,3]/1.e6,color='gold',s=8,edgecolor='black',lw=0,zorder=3,alpha=0.5) # uncorrected
plt.scatter(model2[1:,5],(model2[1:,3]-model2[1:,6])/1.e6,color='gold',s=10,edgecolor='black',label='both free',lw=0.25,zorder=4) # corrected
for i in range(1,len(model2)):
   plt.plot([model2[i,5], model2[i,5]], [model2[i,3]/1.e6, (model2[i,3]-model2[i,6])/1.e6], color='gold', linewidth=0.5, zorder=2, alpha=0.5)
# 400 km 
plt.scatter(model3[1:,5],model3[1:,3]/1.e6,color='darkorange',s=8,edgecolor='black',lw=0,zorder=3,alpha=0.5) # uncorrected
plt.scatter(model3[1:,5],(model3[1:,3]-model3[1:,6])/1.e6,color='darkorange',s=10,edgecolor='black',label='both free',lw=0.25,zorder=4) # corrected
for i in range(1,len(model3)):
   plt.plot([model3[i,5], model3[i,5]], [model3[i,3]/1.e6, (model3[i,3]-model3[i,6])/1.e6], color='darkorange', linewidth=0.5, zorder=2, alpha=0.5)
# 500 km 
plt.scatter(model4[1:,5],model4[1:,3]/1.e6,color='brown',s=8,edgecolor='black',lw=0,zorder=3,alpha=0.5) # uncorrected
plt.scatter(model4[1:,5],(model4[1:,3]-model4[1:,6])/1.e6,color='brown',s=10,edgecolor='black',label='both free',lw=0.25,zorder=4) # corrected
for i in range(1,len(model4)):
   plt.plot([model4[i,5], model4[i,5]], [model4[i,3]/1.e6, (model4[i,3]-model4[i,6])/1.e6], color='brown', linewidth=0.5, zorder=2, alpha=0.5)

# axis stuff
plt.xlim(40,  130); plt.ylim(-30,  50)
# analytical curve
plt.plot(anal_curve[:,0], anal_curve[:,1], color='black', linewidth=1, zorder=1)

# plt.plot([-20, 40], [-20, 40], color='bisque', linewidth=2, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=90, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax2.tick_params(axis='x', labelsize=6)
ax2.tick_params(axis='y', labelsize=6)
plt.xlabel('$\Theta$  [$^\circ$]',size=7.5)
plt.ylabel(r"$\Delta$P$\mathregular{_{extracted}}$  [MPa]",size=7.5)
fixed_aspect_ratio(1)


plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)


