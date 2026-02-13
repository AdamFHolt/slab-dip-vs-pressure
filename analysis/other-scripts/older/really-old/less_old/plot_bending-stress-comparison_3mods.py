#!/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.gridspec import GridSpec
import sys, os, math, statistics
from scipy.signal import savgol_filter
from scipy.interpolate import splrep, splev
from functions import get_average_stress_and_dips

analysis_depth_dz = float(sys.argv[1])  # m (depth interval for shear stress derivative)
ds = float(sys.argv[2])                 # m (distance from slab to pull out DP)
dz = float(sys.argv[3])                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)

model_bothfree = "2D_compositional_subd_lower-res_new"
model_fixedSP  = "2D_compositional_subd_FixedSP_lower-res_new"
model_fixedOP  = "2D_compositional_subd_FixedOP_lower-res_new"

zshall = 230.e3 	# m
zmed   = 330.e3		# m
zdeep  = 430.e3		# m
drho   = 50. 		# kg/m3

text_bothfree_zshall = ''.join(['text_files/',model_bothfree,'.z',str(zshall/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_bothfree_zmed   = ''.join(['text_files/',model_bothfree,'.z',str(zmed/1.e3),  '.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_bothfree_zdeep  = ''.join(['text_files/',model_bothfree,'.z',str(zdeep/1.e3), '.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])

text_fixedSP_zshall  = ''.join(['text_files/',model_fixedSP,'.z',str(zshall/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_fixedSP_zmed    = ''.join(['text_files/',model_fixedSP,'.z',str(zmed/1.e3),  '.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_fixedSP_zdeep   = ''.join(['text_files/',model_fixedSP,'.z',str(zdeep/1.e3), '.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])

text_fixedOP_zshall  = ''.join(['text_files/',model_fixedOP,'.z',str(zshall/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_fixedOP_zmed    = ''.join(['text_files/',model_fixedOP,'.z',str(zmed/1.e3),  '.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_fixedOP_zdeep   = ''.join(['text_files/',model_fixedOP,'.z',str(zdeep/1.e3), '.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])

bothfree_zshall = np.loadtxt((text_bothfree_zshall)) # format: # time [Myr], DP_mod_shall, DP_mod_deep, DP_mod, DP_anal, dip, stress term, stress term B, stress term C
bothfree_zmed   = np.loadtxt((text_bothfree_zmed)) 
bothfree_zdeep  = np.loadtxt((text_bothfree_zdeep))

fixedSP_zshall  = np.loadtxt((text_fixedSP_zshall))
fixedSP_zmed    = np.loadtxt((text_fixedSP_zmed))
fixedSP_zdeep   = np.loadtxt((text_fixedSP_zdeep))

fixedOP_zshall  = np.loadtxt((text_fixedOP_zshall))
fixedOP_zmed    = np.loadtxt((text_fixedOP_zmed))
fixedOP_zdeep   = np.loadtxt((text_fixedOP_zdeep))

plot_name_png = ''.join(['plots/DP-comparisons/compilations/bending-stress.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.png'])

# both free
mean_bothfree_bending_shall, mean_bothfree_dip_shall, mean_bothfree_DP_shall, mean_bothfree_DPanal1_shall, mean_bothfree_DP2_shall = \
	get_average_stress_and_dips(bothfree_zshall,drho)
mean_bothfree_bending_med,   mean_bothfree_dip_med,   mean_bothfree_DP_med,   mean_bothfree_DPanal1_med,   mean_bothfree_DP2_med = \
	get_average_stress_and_dips(bothfree_zmed,drho)
mean_bothfree_bending_deep,  mean_bothfree_dip_deep,  mean_bothfree_DP_deep,  mean_bothfree_DPanal1_deep,  mean_bothfree_DP2_deep = \
	get_average_stress_and_dips(bothfree_zdeep,drho)
sd_bending_bothfree_shall 	= statistics.stdev(bothfree_zshall[2:,6]/1.e6)
sd_bending_bothfree_med   	= statistics.stdev(bothfree_zmed[2:,6]/1.e6)
sd_bending_bothfree_deep  	= statistics.stdev(bothfree_zdeep[2:,6]/1.e6)
sd_dip_bothfree_shall 		= statistics.stdev(bothfree_zshall[2:,5])
sd_dip_bothfree_med   		= statistics.stdev(bothfree_zmed[2:,5])
sd_dip_bothfree_deep  		= statistics.stdev(bothfree_zdeep[2:,5])

# fixed SP
mean_fixedSP_bending_shall, mean_fixedSP_dip_shall, mean_fixedSP_DP_shall, mean_fixedSP_DPanal1_shall, mean_fixedSP_DP2_shall = \
	get_average_stress_and_dips(fixedSP_zshall,drho)
mean_fixedSP_bending_med,   mean_fixedSP_dip_med,   mean_fixedSP_DP_med,   mean_fixedSP_DPanal1_med,   mean_fixedSP_DP2_med = \
	get_average_stress_and_dips(fixedSP_zmed,drho)
mean_fixedSP_bending_deep,  mean_fixedSP_dip_deep,  mean_fixedSP_DP_deep,  mean_fixedSP_DPanal1_deep,  mean_fixedSP_DP2_deep = \
	get_average_stress_and_dips(fixedSP_zdeep,drho)
sd_bending_fixedSP_shall 	= statistics.stdev(fixedSP_zshall[2:,6]/1.e6)
sd_bending_fixedSP_med   	= statistics.stdev(fixedSP_zmed[2:,6]/1.e6)
sd_bending_fixedSP_deep  	= statistics.stdev(fixedSP_zdeep[2:,6]/1.e6)
sd_dip_fixedSP_shall 		= statistics.stdev(fixedSP_zshall[2:,5])
sd_dip_fixedSP_med   		= statistics.stdev(fixedSP_zmed[2:,5])
sd_dip_fixedSP_deep  		= statistics.stdev(fixedSP_zdeep[2:,5])

# fixed OP
mean_fixedOP_bending_shall, mean_fixedOP_dip_shall, mean_fixedOP_DP_shall, mean_fixedOP_DPanal1_shall, mean_fixedOP_DP2_shall = \
	get_average_stress_and_dips(fixedOP_zshall,drho)
mean_fixedOP_bending_med,   mean_fixedOP_dip_med,   mean_fixedOP_DP_med,   mean_fixedOP_DPanal1_med,   mean_fixedOP_DP2_med = \
	get_average_stress_and_dips(fixedOP_zmed,drho)
mean_fixedOP_bending_deep,  mean_fixedOP_dip_deep,  mean_fixedOP_DP_deep,  mean_fixedOP_DPanal1_deep,  mean_fixedOP_DP2_deep = \
	get_average_stress_and_dips(fixedOP_zdeep,drho)
sd_bending_fixedOP_shall 	= statistics.stdev(fixedOP_zshall[2:,6]/1.e6)
sd_bending_fixedOP_med   	= statistics.stdev(fixedOP_zmed[2:,6]/1.e6)
sd_bending_fixedOP_deep  	= statistics.stdev(fixedOP_zdeep[2:,6]/1.e6)
sd_dip_fixedOP_shall 		= statistics.stdev(fixedOP_zshall[2:,5])
sd_dip_fixedOP_med   		= statistics.stdev(fixedOP_zmed[2:,5])
sd_dip_fixedOP_deep  		= statistics.stdev(fixedOP_zdeep[2:,5])

# plotting
fig=plt.figure()
gs=GridSpec(3,1) 

####### DIPS ###################

ax1=fig.add_subplot(gs[0,0])
# fixed SP
pshall =plt.scatter(0.9,mean_fixedSP_dip_shall,color='blue', s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
pmed   =plt.scatter(1,  mean_fixedSP_dip_med,  color='black',s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
pdeep  =plt.scatter(1.1,mean_fixedSP_dip_deep, color='red',  s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
ax1.plot([0.9,0.9], [mean_fixedSP_dip_shall-sd_dip_fixedSP_shall,mean_fixedSP_dip_shall+sd_dip_fixedSP_shall],color='blue',linewidth=5,zorder=2,alpha=0.4)
ax1.plot([1.0,1.0], [mean_fixedSP_dip_med  -sd_dip_fixedSP_med,  mean_fixedSP_dip_med  +sd_dip_fixedSP_med],  color='black',linewidth=5,zorder=2,alpha=0.4)
ax1.plot([1.1,1.1], [mean_fixedSP_dip_deep -sd_dip_fixedSP_deep, mean_fixedSP_dip_deep +sd_dip_fixedSP_deep], color='red', linewidth=5,zorder=2,alpha=0.4)

# legend
leg_string1 = ''.join([str(int(zshall/1.e3)),' km'])
leg_string2 = ''.join([str(int(zmed/1.e3)),' km'])
leg_string3 = ''.join([str(int(zdeep/1.e3)),' km'])
legend = ax1.legend((pshall,pmed,pdeep), (leg_string1,leg_string2,leg_string3),loc='upper left',fontsize='6.5')
legend.get_frame().set_facecolor('white')
legend.get_frame().set_edgecolor('white')

# both free
plt.scatter(1.9,mean_bothfree_dip_shall,color='blue', s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
plt.scatter(2,  mean_bothfree_dip_med,  color='black',s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
plt.scatter(2.1,mean_bothfree_dip_deep, color='red',  s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
ax1.plot([1.9,1.9], [mean_bothfree_dip_shall-sd_dip_bothfree_shall,mean_bothfree_dip_shall+sd_dip_bothfree_shall],color='blue',linewidth=5,zorder=2,alpha=0.4)
ax1.plot([2.0,2.0], [mean_bothfree_dip_med  -sd_dip_bothfree_med,  mean_bothfree_dip_med  +sd_dip_bothfree_med],  color='black',linewidth=5,zorder=2,alpha=0.4)
ax1.plot([2.1,2.1], [mean_bothfree_dip_deep -sd_dip_bothfree_deep, mean_bothfree_dip_deep +sd_dip_bothfree_deep], color='red', linewidth=5,zorder=2,alpha=0.4)

# fixed SP
plt.scatter(2.9,mean_fixedOP_dip_shall,color='blue', s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
plt.scatter(3,  mean_fixedOP_dip_med,  color='black',s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
plt.scatter(3.1,mean_fixedOP_dip_deep, color='red',  s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
ax1.plot([2.9,2.9], [mean_fixedOP_dip_shall-sd_dip_fixedOP_shall,mean_fixedOP_dip_shall+sd_dip_fixedOP_shall],color='blue',linewidth=5,zorder=2,alpha=0.4)
ax1.plot([3.0,3.0], [mean_fixedOP_dip_med  -sd_dip_fixedOP_med,  mean_fixedOP_dip_med  +sd_dip_fixedOP_med],  color='black',linewidth=5,zorder=2,alpha=0.4)
ax1.plot([3.1,3.1], [mean_fixedOP_dip_deep -sd_dip_fixedOP_deep, mean_fixedOP_dip_deep +sd_dip_fixedOP_deep], color='red', linewidth=5,zorder=2,alpha=0.4)

# axis stuff
plt.xlim(0.3,  3.7); plt.ylim(40,  100)
plt.axhline(y=90, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax1.tick_params(axis='x', labelsize=7.5)
ax1.tick_params(axis='y', labelsize=6)
plt.ylabel("dip [$^\circ$]",size=7.5)
ax1.set_xticks( [1,2,3] )
ax1.set_xticklabels( ['','',''] )

####### bending stress ###################

ax2=fig.add_subplot(gs[1,0])
# fixed SP
pshall =plt.scatter(0.9,mean_fixedSP_bending_shall,color='blue', s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
pmed   =plt.scatter(1,  mean_fixedSP_bending_med,  color='black',s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
pdeep  =plt.scatter(1.1,mean_fixedSP_bending_deep, color='red',  s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
ax2.plot([0.9,0.9], [mean_fixedSP_bending_shall-sd_bending_fixedSP_shall,mean_fixedSP_bending_shall+sd_bending_fixedSP_shall],color='blue',linewidth=5,zorder=2,alpha=0.4)
ax2.plot([1.0,1.0], [mean_fixedSP_bending_med  -sd_bending_fixedSP_med,  mean_fixedSP_bending_med  +sd_bending_fixedSP_med],  color='black',linewidth=5,zorder=2,alpha=0.4)
ax2.plot([1.1,1.1], [mean_fixedSP_bending_deep -sd_bending_fixedSP_deep, mean_fixedSP_bending_deep +sd_bending_fixedSP_deep], color='red', linewidth=5,zorder=2,alpha=0.4)

# both free
plt.scatter(1.9,mean_bothfree_bending_shall,color='blue', s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
plt.scatter(2,  mean_bothfree_bending_med,  color='black',s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
plt.scatter(2.1,mean_bothfree_bending_deep, color='red',  s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
ax2.plot([1.9,1.9], [mean_bothfree_bending_shall-sd_bending_bothfree_shall,mean_bothfree_bending_shall+sd_bending_bothfree_shall],color='blue',linewidth=5,zorder=2,alpha=0.4)
ax2.plot([2.0,2.0], [mean_bothfree_bending_med  -sd_bending_bothfree_med,  mean_bothfree_bending_med  +sd_bending_bothfree_med],  color='black',linewidth=5,zorder=2,alpha=0.4)
ax2.plot([2.1,2.1], [mean_bothfree_bending_deep -sd_bending_bothfree_deep, mean_bothfree_bending_deep +sd_bending_bothfree_deep], color='red', linewidth=5,zorder=2,alpha=0.4)

# fixed SP
plt.scatter(2.9,mean_fixedOP_bending_shall,color='blue', s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
plt.scatter(3,  mean_fixedOP_bending_med,  color='black',s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
plt.scatter(3.1,mean_fixedOP_bending_deep, color='red',  s=30,edgecolor='none',lw=0,zorder=3,alpha=1)
ax2.plot([2.9,2.9], [mean_fixedOP_bending_shall-sd_bending_fixedOP_shall,mean_fixedOP_bending_shall+sd_bending_fixedOP_shall],color='blue',linewidth=5,zorder=2,alpha=0.4)
ax2.plot([3.0,3.0], [mean_fixedOP_bending_med  -sd_bending_fixedOP_med,  mean_fixedOP_bending_med  +sd_bending_fixedOP_med],  color='black',linewidth=5,zorder=2,alpha=0.4)
ax2.plot([3.1,3.1], [mean_fixedOP_bending_deep -sd_bending_fixedOP_deep, mean_fixedOP_bending_deep +sd_bending_fixedOP_deep], color='red', linewidth=5,zorder=2,alpha=0.4)

# axis stuff
plt.xlim(0.3,  3.7); plt.ylim(-35,  10)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax2.tick_params(axis='x', labelsize=7.5)
ax2.tick_params(axis='y', labelsize=6)
plt.ylabel("bending stress [MPa]",size=7.5)
ax2.set_xticks( [1,2,3] )
ax2.set_xticklabels( ['','',''] )


############# DP ################################

ax3=fig.add_subplot(gs[2,0])

# fixed SP:
# analytic DP
plt.scatter(0.9,mean_fixedSP_DPanal1_shall,color='blue', s=40,edgecolor='none',lw=0,zorder=3)
panal=plt.scatter(1,  mean_fixedSP_DPanal1_med,  color='black',s=40,edgecolor='none',lw=0,zorder=3)
plt.scatter(1.1,mean_fixedSP_DPanal1_deep, color='red',  s=40,edgecolor='none',lw=0,zorder=3)
# extracted DP
plt.scatter(0.9,mean_fixedSP_DP_shall,color='black', s=30,edgecolor='none',lw=0,zorder=2,alpha=0.4)
pext=plt.scatter(1,  mean_fixedSP_DP_med,  color='black',s=30,edgecolor='none',lw=0,zorder=2,alpha=0.4)
plt.scatter(1.1,mean_fixedSP_DP_deep, color='black',  s=30,edgecolor='none',lw=0,zorder=2,alpha=0.4)
# extracted DP + bending stress
plt.scatter(0.9,mean_fixedSP_DP2_shall,color='blue', s=20,edgecolor='gray',lw=0.75,zorder=5,alpha=1)
pext2=plt.scatter(1,  mean_fixedSP_DP2_med,  color='black',s=20,edgecolor='gray',lw=0.75,zorder=5,alpha=1)
plt.scatter(1.1,mean_fixedSP_DP2_deep, color='red',  s=20,edgecolor='gray',lw=0.75,zorder=5,alpha=1)

# legend
legend = ax3.legend((panal,pext,pext2), ('analytical','extracted','extracted w/bending'),loc='lower left',fontsize='6.5')
legend.get_frame().set_facecolor('white')
legend.get_frame().set_edgecolor('white')

# both free:
# analytic DP
plt.scatter(1.9,mean_bothfree_DPanal1_shall,color='blue', s=40,edgecolor='none',lw=0,zorder=3)
plt.scatter(2,  mean_bothfree_DPanal1_med,  color='black',s=40,edgecolor='none',lw=0,zorder=3)
plt.scatter(2.1,mean_bothfree_DPanal1_deep, color='red',  s=40,edgecolor='none',lw=0,zorder=3)
# extracted DP
plt.scatter(1.9,mean_bothfree_DP_shall,color='black', s=30,edgecolor='none',lw=0,zorder=2,alpha=0.4)
plt.scatter(2,  mean_bothfree_DP_med,  color='black',s=30,edgecolor='none',lw=0,zorder=2,alpha=0.4)
plt.scatter(2.1,mean_bothfree_DP_deep, color='black',  s=30,edgecolor='none',lw=0,zorder=2,alpha=0.4)
# extracted DP + bending stress
plt.scatter(1.9,mean_bothfree_DP2_shall,color='blue', s=20,edgecolor='gray',lw=0.75,zorder=5,alpha=1)
plt.scatter(2,  mean_bothfree_DP2_med,  color='black',s=20,edgecolor='gray',lw=0.75,zorder=5,alpha=1)
plt.scatter(2.1,mean_bothfree_DP2_deep, color='red',  s=20,edgecolor='gray',lw=0.75,zorder=5,alpha=1)

# fixed SP:
# analytic DP
plt.scatter(2.9,mean_fixedOP_DPanal1_shall,color='blue', s=40,edgecolor='none',lw=0,zorder=3)
plt.scatter(3,  mean_fixedOP_DPanal1_med,  color='black',s=40,edgecolor='none',lw=0,zorder=3)
plt.scatter(3.1,mean_fixedOP_DPanal1_deep, color='red',  s=40,edgecolor='none',lw=0,zorder=3)
# extracted DP
plt.scatter(2.9,mean_fixedOP_DP_shall,color='black', s=30,edgecolor='none',lw=0,zorder=2,alpha=0.4)
plt.scatter(3,  mean_fixedOP_DP_med,  color='black',s=30,edgecolor='none',lw=0,zorder=2,alpha=0.4)
plt.scatter(3.1,mean_fixedOP_DP_deep, color='black',  s=30,edgecolor='none',lw=0,zorder=2,alpha=0.4)
# extracted DP + bending stress
plt.scatter(2.9,mean_fixedOP_DP2_shall,color='blue', s=20,edgecolor='gray',lw=0.75,zorder=5,alpha=1)
plt.scatter(3,  mean_fixedOP_DP2_med,  color='black',s=20,edgecolor='gray',lw=0.75,zorder=5,alpha=1)
plt.scatter(3.1,mean_fixedOP_DP2_deep, color='red',  s=20,edgecolor='gray',lw=0.75,zorder=5,alpha=1)

# axis stuff
plt.xlim(0.3,  3.7); plt.ylim(-10,  35)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax3.tick_params(axis='x', labelsize=7.5)
ax3.tick_params(axis='y', labelsize=6)
plt.ylabel(r"$\Delta$P [MPa]",size=7.5)

ax3.set_xticks( [1,2,3] )
ax3.set_xticklabels( ['fixed SP','free plates','fixed OP'] )

plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)


