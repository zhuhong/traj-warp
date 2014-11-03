#!/usr/bin/env python
"""
Trj_modify is used to fit the solute group in trajectory to the center.
@change:
    - 2011.12.26.
        - Add those words.
        - Modified it.
    - 2014.06.23.
        - Update the notes.
"""

from MDAnalysis import Universe, Writer
from MDAnalysis.core.util import greedy_splitext

from MDPackage import Simple_atom
from MDPackage import usage
from MDPackage import Index
# from MDPackage.mymath import least_squares_fitting
import os.path
import numpy as np
import sys
import math
import time as Time

#import psyco ; psyco.jit() 
#from psyco.classes import *


def Usage():
    print "Usage: Trj_modify.py topfile trjfile index.ndx output.xtc"

def Move_2_center(top_file,trj_file,index_file,trjout_file):
    u=Universe(top_file,trj_file)
    TRAJ_FRAMES=u.trajectory.numframes
    w = Writer(trjout_file, u.trajectory.numatoms)

    
    atoms=Simple_atom.Get_Simple_atom_list(top_file)
    index = Index.Read_index_to_Inclass(index_file)
    Index.Print_Index(index)
    while True:
        try:
            solute_index=int(raw_input("Choosing the group for centering:"))
            # solvent_index = int(raw_input("Choosing the solvent group:"))
            break
        except:
            print "You should input a number."
            continue
    solute_group  = index[solute_index].group_list
    # solvent_group = index[solvent_index].group_list
    solute_atoms  =len(solute_group)
    NUM_ATOMS = u.trajectory.numatoms
    # print "\t Reading %d frames from trajectory file: %s" %(nframes,traj_file)    
    START_TIME=Time.time()
    for ts in u.trajectory:
        ref_com =np.zeros((3),dtype=np.float32)    

        # sys.stderr.write("\t Reading frame %8d\r" %ts.frame)
        # sys.stderr.flush()

        for i,ai in list(enumerate(solute_group)):
            ref_com[0] += ts._x[ai-1]
            ref_com[1] += ts._y[ai-1]
            ref_com[2] += ts._z[ai-1]
        ref_com = ref_com/solute_atoms
        dimensions=ts.dimensions
        ref_com = ref_com - np.array([dimensions[0]/2, dimensions[1]/2, dimensions[2]/2])  

        
        # ref_com = np.array([sum(traj_data[:,0])/solute_atoms - dimensions[0]/2,\
            # sum(traj_data[:,1])/solute_atoms - dimensions[1]/2,\
            # sum(traj_data[:,2])/solute_atoms - dimensions[2]/2])  

        for i in range(NUM_ATOMS):
            ts._x[i] = ts._x[i] - ref_com[0] 
            ts._y[i] = ts._y[i] - ref_com[1] 
            ts._z[i] = ts._z[i] - ref_com[2] 

            if (i+1) in solute_group:
                continue

            if ts._x[i] > dimensions[0] or ts._x[i] <0:
                ts._x[i]=ts._x[i]%dimensions[0]

            if ts._y[i] > dimensions[1] or ts._y[i] <0:
                ts._y[i]=ts._y[i]%dimensions[1]

            if ts._z[i] > dimensions[2] or ts._z[i] < 0:
                ts._z[i]=ts._z[i]%dimensions[2]

        NOW_TIME=Time.time()
        BIN_TIME=NOW_TIME-START_TIME
        # if ts.frame % 10 == 0:
#            usage.echo("%8.4f   %8.4f   %8.4f\r" %(dimensions[0],dimensions[1],dimensions[2]))
        usage.echo(" "*40+"Converted frame %d, time used: %8.2f s, time left: %8.2f s \r" \
            % (ts.frame,BIN_TIME,BIN_TIME*(float(TRAJ_FRAMES)/ts.frame-1) ))
#    for ts in u.trajectory:
        w.write(ts)
        # del traj_data
#        usage.echo("Writting frame %d\r"  %ts.frame)
    w.close_trajectory()
    print "Converted %r --> %r" % (intrj, outtrj)



if __name__=="__main__":
    if len(sys.argv)!=5:
        Usage()
        sys.exit()

    topol =sys.argv[1] #PRMpbc
    intrj =sys.argv[2] #TRJpbc_bz2
    index = sys.argv[3]
    outtrj = sys.argv[4]
    Move_2_center(topol,intrj,index,outtrj)    




