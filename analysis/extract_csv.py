#!/home/holt/software/ParaView-5.11.1-MPI-Linux-Python3.9-x86_64/bin/pvpython
#
# extract csv files for all timesteps of a .pvd file
#
from paraview.simple import *
import os,sys
import numpy as np
mod_name=str(sys.argv[1])
max_time=int(sys.argv[2])
models_loc = 'raw_outputs/'
output_loc = 'csv_outputs/'
solution = ''.join([models_loc,str(mod_name),'/solution.pvd'])
reader = OpenDataFile(solution)

# directory for csv files
if not os.path.exists(''.join([output_loc,str(mod_name)])):
        os.mkdir(''.join([output_loc,str(mod_name)]))

times = reader.TimestepValues

#for i in np.array([max_time]):
for i in range(max_time):

    time = times[i]
    ofile = 'full.%d.csv' % i
    ofiled=''.join([output_loc,str(mod_name),'/',ofile])
    print('extracting',i,time,ofiled)
    writer = CreateWriter( ofiled, reader)
    writer.FieldAssociation = "Point Data" # or "Cells"
    writer.UpdatePipeline(time)
    del writer
  

