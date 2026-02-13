#!/usr/bin/env python3 
#cp'd from make_tempSPandOP_halfspace80Ma25Ma_rad200km_BigBox_SmallerOP.py

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.ioff()

ofile="outputs/temp_constant.txt"

Tmant = 1573.

# box dimensions (i.e. "extent" in ASPECT input)
xmin=0;xmax=5800.e3;
ymin=0;ymax=1450.e3;
# number of input field points
xnum=5800
ynum=1450

No_nodes= (xnum + 1) * (ynum + 1)
T=np.zeros([No_nodes,3],float)
 
ind=0
print("writting text file...")
for j in range(ynum + 1): 
	for i in range(xnum + 1):

		x = xmin + i * ((xmax - xmin)/xnum)
		y = ymin + j * ((ymax - ymin)/ynum) 

		T[ind,0] = x
		T[ind,1] = y
		T[ind,2] = Tmant

		ind=ind+1;
 
# write to file
f= open(ofile,"w+")
f.write("# POINTS: %s %s\n" % (str(xnum+1),str(ynum+1)))
f.write("# Columns: x y temperature\n")
for k in range(0,ind):
	f.write("%.6f %.6f %.6f\n" % (T[k,0],T[k,1],T[k,2]))
f.close() 