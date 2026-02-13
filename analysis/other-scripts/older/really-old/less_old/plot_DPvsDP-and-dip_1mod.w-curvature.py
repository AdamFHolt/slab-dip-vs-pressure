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

anal_depth1 = 230.0e3
anal_depth2 = 330.0e3
anal_depth3 = 430.0e3

text_model1 = ''.join(['text_files/',model_name,'.z',str(anal_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_model2 = ''.join(['text_files/',model_name,'.z',str(anal_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_model3 = ''.join(['text_files/',model_name,'.z',str(anal_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])

#bothfree = np.loadtxt((text_bothfree)) # format: # time [Myr], DP_mod_shall, DP_mod_deep, DP_mod, DP_anal, dip, stress term, stress term B, stress term C
model1  = np.loadtxt((text_model1))
model2  = np.loadtxt((text_model2))
model3  = np.loadtxt((text_model3))

plot_loc   = 'plots/DP-comparisons/single-models/'
plot_name_png = ''.join([plot_loc,'/',model_name,'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.w-curvature.png'])
plot_name_pdf = ''.join([plot_loc,'/',model_name,'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.w-curvature.pdf'])

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

ax=fig.add_subplot(gs[0,0])
curvature_thresh = 0.05
for i in range(2,len(model1)):
	if np.abs(model1[i,12]) < curvature_thresh:	
		plt.scatter(model1[i,4]/1.e6,model1[i,3]/1.e6,color='blue',s=15,lw=0.0,zorder=3,alpha=0.6) # uncorrected
		plt.scatter(model1[i,4]/1.e6,(model1[i,3]-model1[i,6])/1.e6,color='blue',s=15,edgecolor='black',label='230 km',lw=0.25,zorder=4) # corrected		
		plt.plot([model1[i,4]/1.e6, model1[i,4]/1.e6], [model1[i,3]/1.e6, (model1[i,3]-model1[i,6])/1.e6], color='blue', linewidth=1, zorder=2, alpha=0.5)
	else:
		pass

# intermediate
for i in range(2,len(model2)):
	if np.abs(model2[i,12]) < curvature_thresh:	
		plt.scatter(model2[i,4]/1.e6,model2[i,3]/1.e6,color='black',s=15,lw=0.0,zorder=3,alpha=0.6) # uncorrected
		plt.scatter(model2[i,4]/1.e6,(model2[i,3]-model2[i,6])/1.e6,color='black',s=15,edgecolor='black',label='330 km',lw=0.25,zorder=4) # corrected		
		plt.plot([model2[i,4]/1.e6, model2[i,4]/1.e6], [model2[i,3]/1.e6, (model2[i,3]-model2[i,6])/1.e6], color='black', linewidth=1, zorder=2, alpha=0.5)
	else:
		pass

# deep
for i in range(2,len(model3)):
	if np.abs(model3[i,12]) < curvature_thresh:	
		plt.scatter(model3[i,4]/1.e6,model3[i,3]/1.e6,color='red',s=15,lw=0.0,zorder=3,alpha=0.6) # uncorrected
		plt.scatter(model3[i,4]/1.e6,(model3[i,3]-model3[i,6])/1.e6,color='red',s=15,edgecolor='black',label='430 km',lw=0.25,zorder=4) # corrected		
		plt.plot([model3[i,4]/1.e6, model3[i,4]/1.e6], [model3[i,3]/1.e6, (model3[i,3]-model3[i,6])/1.e6], color='red', linewidth=1, zorder=2, alpha=0.5)
	else:
		pass

# axis stuff
plt.xlim(-17.5,  40); plt.ylim(-17.5,  40)
plt.plot([-17.5, 40], [-17.5, 40], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
# plt.xlabel("slab buoyancy  [MPa]",size=7.5)
plt.ylabel("extracted stress  [MPa]",size=7.5)
# plt.legend(fontsize=5.5,loc='lower right',facecolor='white',fancybox=True, framealpha=1, frameon=False)
fixed_aspect_ratio(1)
ax.annotate('low curvature', xy=(0.0125,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=9,color='k') 

anal_curve = np.zeros((140,2))
for j in range(len(anal_curve)):
	anal_curve[j,0] = j # dip [deg]
	drho = 50. # kg/m3
	slabnorm_thick = 80.e3 # m  -- CHECK!!
	anal_curve[j,1] = DP_anal = (drho * 9.81 * slabnorm_thick * np.cos(np.deg2rad(j)))/1.e6 # DP  [MPa]

ax2=fig.add_subplot(gs[0,1])
# shallow
for i in range(2,len(model1)):
	if np.abs(model1[i,12]) < curvature_thresh:
	    plt.scatter(model1[i,5],model1[i,3]/1.e6,color='blue',s=15,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
	    plt.scatter(model1[i,5],(model1[i,3]-model1[i,6])/1.e6,color='blue',s=15,edgecolor='black',lw=0.25,zorder=4) # corrected
	    plt.plot([model1[i,5], model1[i,5]], [model1[i,3]/1.e6, (model1[i,3]-model1[i,6])/1.e6], color='blue', linewidth=0.5, zorder=2, alpha=0.5)
	else:
		pass

# intermediate
for i in range(2,len(model2)):
	if np.abs(model2[i,12]) < curvature_thresh:
	    plt.scatter(model2[i,5],model2[i,3]/1.e6,color='black',s=15,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
	    plt.scatter(model2[i,5],(model2[i,3]-model2[i,6])/1.e6,color='black',s=15,edgecolor='black',lw=0.25,zorder=4) # corrected
	    plt.plot([model2[i,5], model2[i,5]], [model2[i,3]/1.e6, (model2[i,3]-model2[i,6])/1.e6], color='black', linewidth=0.5, zorder=2, alpha=0.5)
	else:
		pass

# deep
for i in range(2,len(model3)):
	if np.abs(model3[i,12]) < curvature_thresh:
	    plt.scatter(model3[i,5],model3[i,3]/1.e6,color='red',s=15,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
	    plt.scatter(model3[i,5],(model3[i,3]-model3[i,6])/1.e6,color='red',s=15,edgecolor='black',lw=0.25,zorder=4) # corrected
	    plt.plot([model3[i,5], model3[i,5]], [model3[i,3]/1.e6, (model3[i,3]-model3[i,6])/1.e6], color='red', linewidth=0.5, zorder=2, alpha=0.5)
	else:
		pass

# axis stuff
plt.xlim(40,  102.5); plt.ylim(-17.5,  50)
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


ax3=fig.add_subplot(gs[1,0])
for i in range(2,len(model1)):
	if np.abs(model1[i,12]) > curvature_thresh:	
		plt.scatter(model1[i,4]/1.e6,model1[i,3]/1.e6,color='blue',s=15,lw=0.0,zorder=3,alpha=0.6) # uncorrected
		plt.scatter(model1[i,4]/1.e6,(model1[i,3]-model1[i,6])/1.e6,color='blue',s=15,edgecolor='black',label='230 km',lw=0.25,zorder=4) # corrected		
		plt.plot([model1[i,4]/1.e6, model1[i,4]/1.e6], [model1[i,3]/1.e6, (model1[i,3]-model1[i,6])/1.e6], color='blue', linewidth=1, zorder=2, alpha=0.5)
	else:
		pass

# intermediate
for i in range(2,len(model2)):
	if np.abs(model2[i,12]) > curvature_thresh:	
		plt.scatter(model2[i,4]/1.e6,model2[i,3]/1.e6,color='black',s=15,lw=0.0,zorder=3,alpha=0.6) # uncorrected
		plt.scatter(model2[i,4]/1.e6,(model2[i,3]-model2[i,6])/1.e6,color='black',s=15,edgecolor='black',label='330 km',lw=0.25,zorder=4) # corrected		
		plt.plot([model2[i,4]/1.e6, model2[i,4]/1.e6], [model2[i,3]/1.e6, (model2[i,3]-model2[i,6])/1.e6], color='black', linewidth=1, zorder=2, alpha=0.5)
	else:
		pass

# deep
for i in range(2,len(model3)):
	if np.abs(model3[i,12]) > curvature_thresh:	
		plt.scatter(model3[i,4]/1.e6,model3[i,3]/1.e6,color='red',s=15,lw=0.0,zorder=3,alpha=0.6) # uncorrected
		plt.scatter(model3[i,4]/1.e6,(model3[i,3]-model3[i,6])/1.e6,color='red',s=15,edgecolor='black',label='430 km',lw=0.25,zorder=4) # corrected		
		plt.plot([model3[i,4]/1.e6, model3[i,4]/1.e6], [model3[i,3]/1.e6, (model3[i,3]-model3[i,6])/1.e6], color='red', linewidth=1, zorder=2, alpha=0.5)
	else:
		pass


# axis stuff
plt.xlim(-17.5,  40); plt.ylim(-17.5,  40)
plt.plot([-17.5, 40], [-17.5, 40], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax3.tick_params(axis='x', labelsize=6)
ax3.tick_params(axis='y', labelsize=6)
plt.xlabel("slab buoyancy  [MPa]",size=7.5)
plt.ylabel("extracted stress  [MPa]",size=7.5)
# plt.legend(fontsize=5.5,loc='lower right',facecolor='white',fancybox=True, framealpha=1, frameon=False)
fixed_aspect_ratio(1)
ax3.annotate('high curvature', xy=(0.0125,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=9,color='k') 


anal_curve = np.zeros((140,2))
for j in range(len(anal_curve)):
	anal_curve[j,0] = j # dip [deg]
	drho = 50. # kg/m3
	slabnorm_thick = 80.e3 # m  -- CHECK!!
	anal_curve[j,1] = DP_anal = (drho * 9.81 * slabnorm_thick * np.cos(np.deg2rad(j)))/1.e6 # DP  [MPa]

ax4=fig.add_subplot(gs[1,1])
# shallow
for i in range(2,len(model1)):
	if np.abs(model1[i,12]) > curvature_thresh:
	    plt.scatter(model1[i,5],model1[i,3]/1.e6,color='blue',s=15,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
	    plt.scatter(model1[i,5],(model1[i,3]-model1[i,6])/1.e6,color='blue',s=15,edgecolor='black',lw=0.25,zorder=4) # corrected
	    plt.plot([model1[i,5], model1[i,5]], [model1[i,3]/1.e6, (model1[i,3]-model1[i,6])/1.e6], color='blue', linewidth=0.5, zorder=2, alpha=0.5)
	else:
		pass

# intermediate
for i in range(2,len(model2)):
	if np.abs(model2[i,12]) > curvature_thresh:
	    plt.scatter(model2[i,5],model2[i,3]/1.e6,color='black',s=15,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
	    plt.scatter(model2[i,5],(model2[i,3]-model2[i,6])/1.e6,color='black',s=15,edgecolor='black',lw=0.25,zorder=4) # corrected
	    plt.plot([model2[i,5], model2[i,5]], [model2[i,3]/1.e6, (model2[i,3]-model2[i,6])/1.e6], color='black', linewidth=0.5, zorder=2, alpha=0.5)
	else:
		pass

# deep
for i in range(2,len(model3)):
	if np.abs(model3[i,12]) > curvature_thresh:
	    plt.scatter(model3[i,5],model3[i,3]/1.e6,color='red',s=15,edgecolor='black',lw=0,zorder=3,alpha=0.6) # uncorrected
	    plt.scatter(model3[i,5],(model3[i,3]-model3[i,6])/1.e6,color='red',s=15,edgecolor='black',lw=0.25,zorder=4) # corrected
	    plt.plot([model3[i,5], model3[i,5]], [model3[i,3]/1.e6, (model3[i,3]-model3[i,6])/1.e6], color='red', linewidth=0.5, zorder=2, alpha=0.5)
	else:
		pass

# axis stuff
plt.xlim(40,  102.5); plt.ylim(-17.5,  50)
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


