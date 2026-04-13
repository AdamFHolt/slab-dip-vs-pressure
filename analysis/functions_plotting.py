#!/bin/python
import numpy as np
import scipy.special
import scipy.integrate
from scipy.interpolate import griddata
from scipy.signal import savgol_filter
import sys, os, math, statistics, subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def plot_BvsFullForce_wKthresh(tmin,model,curvature_thresh,center_color,edge_color,marker):

	for i in range(tmin,len(model)):
		if np.abs(model[i,11])*1000 < curvature_thresh: 
			plt.scatter(model[i,4]/1.e6,(model[i,3]-model[i,6]+model[i,17])/1.e6,color=center_color,s=10,edgecolor=edge_color,marker=marker,lw=0.1,zorder=3) 
		else:
			plt.scatter(model[i,4]/1.e6,(model[i,3]-model[i,6]+model[i,17])/1.e6,color=center_color,s=5,edgecolor=edge_color,marker=marker,lw=0.1,zorder=3,alpha=0.2) 

def plot_BvsFullForce_wKthresh_overturned(tmin,model,curvature_thresh,edge_color,marker):

	for i in range(tmin,len(model)):
		if np.abs(model[i,11])*1000 < curvature_thresh: 
			plt.scatter(model[i,4]/1.e6,(model[i,3]-model[i,6]+model[i,17])/1.e6,facecolors='none',s=10,edgecolor=edge_color,marker=marker,lw=0.2,zorder=2,alpha=0.7) 
		else:
			plt.scatter(model[i,4]/1.e6,(model[i,3]-model[i,6]+model[i,17])/1.e6,facecolors='none',s=5,edgecolor=edge_color,marker=marker,lw=0.2,zorder=2,alpha=0.7) 

def plot_BvsFullForce_Kcolored(tmin,model,edge_color,marker,color_map,norm,zorder):

	for i in range(tmin,len(model)):
		plt.scatter(model[i,4]/1.e6,(model[i,3]-model[i,6]+model[i,17])/1.e6,s=10,c=np.abs(model[i,11])*1e3,cmap=color_map,norm=norm,edgecolor=edge_color,linewidth=0.1,zorder=zorder,marker=marker)

def plot_BvsFullForce_Kcolored_overturned(
    tmin, model, marker, color_map, norm, zorder
):

    for i in range(tmin, len(model)):

        x = model[i,4] / 1.e6
        y = (model[i,3] - model[i,6] + model[i,17]) / 1.e6

        val = np.abs(model[i,11]) * 1e3
        edge_col = color_map(norm(val))  # ← map to RGBA

        plt.scatter(
            x, y,
            s=10,
            facecolors='none',        # hollow marker
            edgecolors=edge_col,      # colored outline
            linewidth=0.4,
            zorder=zorder,
            marker=marker
        )


def plot_BvsDP_wKthresh(tmin,model,curvature_thresh,center_color,edge_color,marker):

	for i in range(tmin,len(model)):
		if np.abs(model[i,11])*1000 < curvature_thresh: 
			plt.scatter(model[i,4]/1.e6,model[i,3]/1.e6,color=center_color,s=10,edgecolor=edge_color,marker=marker,lw=0.1,zorder=3) 
		else:
			plt.scatter(model[i,4]/1.e6,model[i,3]/1.e6,color=center_color,s=5,edgecolor=edge_color,marker=marker,lw=0.1,zorder=3,alpha=0.2)

def plot_BvsDP_wKthresh_overturned(tmin,model,curvature_thresh,edge_color,marker):

	for i in range(tmin,len(model)):
		if np.abs(model[i,11])*1000 < curvature_thresh: 
			plt.scatter(model[i,4]/1.e6,model[i,3]/1.e6,facecolors='none',s=10,edgecolor=edge_color,marker=marker,lw=0.2,zorder=2,alpha=0.7) 
		else:
			plt.scatter(model[i,4]/1.e6,model[i,3]/1.e6,facecolors='none',s=5,edgecolor=edge_color,marker=marker,lw=0.2,zorder=2,alpha=0.7)


def plot_BvsDP_scalingcolored(tmin,model,edge_color,marker,color_map,norm,mant_visc,zorder,viscosity):

	cmyr_to_ms = 0.01/(365.25*24*3600)

	for i in range(tmin,len(model)):

		plt.scatter(model[i,4]/1.e6,model[i,3]/1.e6,s=10,c=(100./1624.)*model[i,11]*viscosity*mant_visc*model[i,19]*cmyr_to_ms*1e-6,cmap=color_map,norm=norm,edgecolor=edge_color,linewidth=0.1,zorder=zorder,marker=marker)

def plot_BvsDP_scalingcolored_overturned(
    tmin, model, marker, color_map, norm, mant_visc, zorder, viscosity
):
    cmyr_to_ms = 0.01 / (365.25 * 24 * 3600)

    for i in range(tmin, len(model)):
        x = model[i, 4] / 1.e6
        y = model[i, 3] / 1.e6

        # value you want to color by (same expression you had in c=...)
        val = (
            (100./1624.) * model[i, 11] * viscosity * mant_visc * model[i, 19] * cmyr_to_ms * 1e-6
        )

        edge_col = color_map(norm(val))  # map -> RGBA using norm + cmap

        plt.scatter(
            x, y,
            s=10,
            facecolors='none',        # transparent center
            edgecolors=edge_col,      # colored edge from colormap
            linewidth=0.4,
            zorder=zorder,
            marker=marker
        )


def plot_forcecomponent_fullstressmisfit(tmin,model,curve_thresh,center_color,marker,misfit_color):

	K_ind    = 11; ss_ind   = 17; sn_ind   = 6
	anal_ind = 4; DP_ind   = 3

	for i in range(tmin,len(model)):
		if model[i,K_ind]*1e3 >= curve_thresh:
			plt.scatter(model[i,K_ind]*1e3, (model[i,DP_ind]-model[i,anal_ind]+model[i,ss_ind]-model[i,sn_ind])/1.e6,s=10,color=center_color,edgecolor=misfit_color,linewidth=0.25,zorder=3,marker=marker)
		else:
			plt.scatter(model[i,K_ind]*1e3, (model[i,DP_ind]-model[i,anal_ind]+model[i,ss_ind]-model[i,sn_ind])/1.e6,s=10,color=center_color,edgecolor='black',linewidth=0.25,zorder=3,marker=marker)




def plot_forcecomponent_fullstressmisfit_overturned(tmin,model,center_color,marker):

	K_ind    = 11; ss_ind   = 17; sn_ind   = 6
	anal_ind = 4; DP_ind   = 3

	for i in range(tmin,len(model)):
		plt.scatter(model[i,K_ind]*1e3, (model[i,DP_ind]-model[i,anal_ind]+model[i,ss_ind]-model[i,sn_ind])/1.e6,s=10,facecolors='none',edgecolor=center_color,linewidth=0.25,zorder=3,marker=marker,alpha=0.7)



def append_DPcomponents(tmin,model,model_pref,DP_array):

	K_ind    = 11; DP_ind   = 3; Psp_ind = 15; Pop_ind = 16;
	
	for i in range(tmin,min(len(model),len(model_pref))):
		a = model[i,DP_ind]
		b = (model[i,Psp_ind]-model_pref[i,1])/model[i,DP_ind]
		DP_array = np.vstack([DP_array, np.array([[a, b]])])

	return DP_array


def plot_forcecomponent_dpmisfit(tmin,model,curve_thresh,x_ind,center_color,marker,misfit_color):

	K_ind    = 11; anal_ind = 4; DP_ind   = 3
	for i in range(tmin,len(model)):
		if model[i,K_ind]*1e3 >= curve_thresh:
			plt.scatter(model[i,x_ind], -1.0*(model[i,DP_ind]-model[i,anal_ind])/1.e6,s=10,color=center_color,edgecolor=misfit_color,linewidth=0.25,zorder=3,marker=marker)
		else:
			plt.scatter(model[i,x_ind], -1.0*(model[i,DP_ind]-model[i,anal_ind])/1.e6,s=10,color=center_color,edgecolor='black',linewidth=0.25,zorder=3,marker=marker)



def plot_forcecomponent_dpmisfit_overturned(tmin,model,x_ind,center_color,marker):

	K_ind    = 11; anal_ind = 4; DP_ind   = 3
	for i in range(tmin,len(model)):
		plt.scatter(model[i,x_ind], -1.0*(model[i,DP_ind]-model[i,anal_ind])/1.e6,s=10,facecolors='none',edgecolor=center_color,linewidth=0.25,zorder=3,marker=marker,alpha=0.7)


def plot_forcecomponent_dqds(tmin,model,curve_thresh,x_ind,center_color,marker,misfit_color):

	K_ind    = 11; anal_ind = 4; DP_ind   = 3
	ss_ind   = 17; sn_ind   = 6

	for i in range(tmin,len(model)):
		if model[i,K_ind]*1e3 >= curve_thresh:
			plt.scatter(model[i,x_ind], -1.0*model[i,sn_ind]/1.e6,s=10,color=center_color,edgecolor=misfit_color,linewidth=0.25,zorder=3,marker=marker)
		else:
			plt.scatter(model[i,x_ind], -1.0*model[i,sn_ind]/1.e6,s=10,color=center_color,edgecolor='black',linewidth=0.25,zorder=3,marker=marker)


def plot_forcecomponent_dqds_overturned(tmin,model,x_ind,center_color,marker):

	K_ind    = 11; anal_ind = 4; DP_ind   = 3
	ss_ind   = 17; sn_ind   = 6

	for i in range(tmin,len(model)):
		plt.scatter(model[i,x_ind], -1.0*model[i,sn_ind]/1.e6,s=10,facecolors='none',edgecolor=center_color,linewidth=0.25,zorder=3,marker=marker,alpha=0.7)


def plot_forcecomponent_dqds_vsK(tmin,model,curve_thresh,x_ind,center_color,marker,misfit_color):

	K_ind    = 11; anal_ind = 4; DP_ind   = 3
	ss_ind   = 17; sn_ind   = 6

	for i in range(tmin,len(model)):
		if model[i,K_ind]*1e3 >= curve_thresh:
			plt.scatter(model[i,K_ind]*1e3, -1.0*model[i,sn_ind]/1.e6,s=10,color=center_color,edgecolor=misfit_color,linewidth=0.25,zorder=3,marker=marker)
		else:
			plt.scatter(model[i,K_ind]*1e3, -1.0*model[i,sn_ind]/1.e6,s=10,color=center_color,edgecolor='black',linewidth=0.25,zorder=3,marker=marker)


def plot_forcecomponent_dqds_vsK_overturned(tmin,model,x_ind,center_color,marker):

	K_ind    = 11; anal_ind = 4; DP_ind   = 3
	ss_ind   = 17; sn_ind   = 6

	for i in range(tmin,len(model)):
		plt.scatter(model[i,K_ind]*1e3, -1.0*model[i,sn_ind]/1.e6,s=10,facecolors='none',edgecolor=center_color,linewidth=0.25,zorder=3,marker=marker,alpha=0.7)


def plot_forcecomponent_dqds_vsVc(tmin,model,curve_thresh,x_ind,center_color,marker,misfit_color):

	K_ind    = 11; anal_ind = 4; DP_ind   = 3
	ss_ind   = 17; sn_ind   = 6; vc_ind   = 19

	for i in range(tmin,len(model)):
		if model[i,K_ind]*1e3 >= curve_thresh:
			plt.scatter(model[i,vc_ind], -1.0*model[i,sn_ind]/1.e6,s=10,color=center_color,edgecolor=misfit_color,linewidth=0.25,zorder=3,marker=marker)
		else:
			plt.scatter(model[i,vc_ind], -1.0*model[i,sn_ind]/1.e6,s=10,color=center_color,edgecolor='black',linewidth=0.25,zorder=3,marker=marker)


def plot_forcecomponent_dqds_vsVc_overturned(tmin,model,x_ind,center_color,marker):

	K_ind    = 11; anal_ind = 4; DP_ind   = 3
	ss_ind   = 17; sn_ind   = 6; vc_ind   = 19

	for i in range(tmin,len(model)):
		plt.scatter(model[i,vc_ind], -1.0*model[i,sn_ind]/1.e6,s=10,facecolors='none',edgecolor=center_color,linewidth=0.25,zorder=3,marker=marker,alpha=0.7)



def plot_forcecomponent_dqds_vsVisc(tmin,model,curve_thresh,x_ind,center_color,marker,misfit_color,viscosity):

	K_ind    = 11; anal_ind = 4; DP_ind   = 3
	ss_ind   = 17; sn_ind   = 6; vc_ind   = 19

	for i in range(tmin,len(model)):
		if model[i,K_ind]*1e3 >= curve_thresh:
			plt.scatter(viscosity, -1.0*model[i,sn_ind]/1.e6,s=10,color=center_color,edgecolor=misfit_color,linewidth=0.25,zorder=3,marker=marker)
		else:
			plt.scatter(viscosity, -1.0*model[i,sn_ind]/1.e6,s=10,color=center_color,edgecolor='black',linewidth=0.25,zorder=3,marker=marker)


def plot_forcecomponent_dqds_vsVisc_overturned(tmin,model,x_ind,center_color,marker,viscosity):

	K_ind    = 11; anal_ind = 4; DP_ind   = 3
	ss_ind   = 17; sn_ind   = 6; vc_ind   = 19

	for i in range(tmin,len(model)):
		plt.scatter(viscosity, -1.0*model[i,sn_ind]/1.e6,s=10,facecolors='none',edgecolor=center_color,linewidth=0.25,zorder=3,marker=marker,alpha=0.7)


def plot_forcecomponent_KN(tmin,model,curve_thresh,center_color,marker,misfit_color):

	K_ind    = 11; anal_ind = 4; DP_ind   = 3
	ss_ind   = 17; sn_ind   = 6

	for i in range(tmin,len(model)):
		if model[i,K_ind]*1e3 >= curve_thresh:
			plt.scatter(model[i,K_ind]*1e3, model[i,ss_ind]/1.e6,s=10,color=center_color,edgecolor=misfit_color,linewidth=0.25,zorder=3,marker=marker)
		else:
			plt.scatter(model[i,K_ind]*1e3, model[i,ss_ind]/1.e6,s=10,color=center_color,edgecolor='black',linewidth=0.25,zorder=3,marker=marker)


def plot_dK(tmin,model,curve_thresh,center_color,marker,misfit_color):

	K_ind    = 11; anal_ind = 4; DP_ind   = 3
	dK_ind   = 12; sn_ind   = 6

	for i in range(tmin,len(model)):
		if model[i,K_ind]*1e3 >= curve_thresh:
			plt.scatter(model[i,K_ind]*1e3, model[i,dK_ind]*1e6,s=10,color=center_color,edgecolor=misfit_color,linewidth=0.25,zorder=3,marker=marker)
		else:
			plt.scatter(model[i,K_ind]*1e3, model[i,dK_ind]*1.e6,s=10,color=center_color,edgecolor='black',linewidth=0.25,zorder=3,marker=marker)


def plot_forcecomponent_N(tmin,model,curve_thresh,center_color,marker,misfit_color):

	K_ind    = 11; anal_ind = 4; DP_ind   = 3
	ss_ind   = 17; sn_ind   = 6

	for i in range(tmin,len(model)):
		if model[i,K_ind]*1e3 >= curve_thresh:
			plt.scatter(model[i,K_ind]*1e3, (model[i,ss_ind]/model[i,K_ind])/1.e6,s=10,color=center_color,edgecolor=misfit_color,linewidth=0.25,zorder=3,marker=marker)
		else:
			plt.scatter(model[i,K_ind]*1e3, (model[i,ss_ind]/model[i,K_ind])/1.e6,s=10,color=center_color,edgecolor='black',linewidth=0.25,zorder=3,marker=marker)

