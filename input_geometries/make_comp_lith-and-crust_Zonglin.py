#!/usr/bin/env python 
import sys
import numpy
import matplotlib.pyplot as plt

## to do:
## - add a straight portion of the slab (70 deg dip, to 500 km depth)
## - crust also on the top of the overriding palte
## - same thickness SP and OP
## - extend to 3-D
##		- box size = 6000 * 6000 * 1000 (depth)
##      - SP size  = 2000 * 2000
##      - OP size  = 2000 * 1000
## - make a Python script to plot

ofile="outputs/comp_lith-and-crust_Zonglin.txt"

# box dimensions (i.e. "extent" in ASPECT input)
xmin=0;xmax=5800.e3;
ymin=0;ymax=1450.e3;
# number of cells, (i.e. "number of repetitions" in ASPECT input)
xnum= 580*4
ynum= 145*4

x_gap = 500.e3; 
x_SP  = 3000.e3; 
y_crust = 8e3;
depth_notch  = 200e3;
radius_outer = 245e3;
slab_dip = 70.;
OPthick = 50.e3;
SPthick = 100.e3;

No_nodes= (xnum + 1) * (ynum + 1)
C=numpy.zeros([No_nodes,4],float)
 
ind=0

for j in range(ynum + 1): 
	for i in range(xnum + 1):

		x = xmin + i * ((xmax - xmin)/xnum)
		y = ymin + j * ((ymax - ymin)/ynum) 
  
		C[ind,0] = x
		C[ind,1] = y

		# crust along top of flat portion of SP
		if x > (x_gap) and x <= (x_gap + x_SP - radius_outer):
			if y >= (ymax - y_crust):
				C[ind,2]=1
			# top half of the lithosphere
			elif y < (ymax - y_crust) and y > (ymax - (SPthick)):
				C[ind,3]=1

		# curved portion of lithosphere/crust ("notch")
		elif x > (x_gap C[ind,3]=1

		# overriding plate (OP) above notch
		if x > (x_gap + x_SP - radius_outer) and x < (x_gap + x_SP):
			x1 = x_gap + x_SP - radius_outer; 
			y1 = ymax - radius_outer;
			if ((x-x1)**2 + (y-y1)**2) >= radius_outer**2 and y > (ymax - OPthick): 
				C[ind,3]= 1
				
		# rest of the OP
		if  x >= (x_gap + x_SP) and x < (xmax - x_gap) and y > (ymax - OPthick): 
			C[ind,3]= 1

		ind=ind+1;
 

# write to file
f= open(ofile,"w+")
f.write("# POINTS: %s %s\n" % (str(xnum+1),str(ynum+1)))
f.write("# Columns: x y composition1 composition2\n")
for k in range(0,ind):
	f.write("%.6f %.6f %.2f %.2f\n" % (C[k,0],C[k,1],C[k,2],C[k,3]))
f.close() 

