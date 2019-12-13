#
# Copyright (C) 2019 Lei Liu & Changbong Hyeon
# This file is part of HLM-epidom.
#
import os
import sys
import h5py
from numpy import *
from itertools import islice

def l2f(List): return map(float, List)

# Main Func. #################################
# change the coordinates format from blockfile to HDF5
if not len(sys.argv) == 5: 
    print("usage:: python gnm.bk2h5.py xxx.equ ts te NFRAME") 
    sys.exit() 
# 
bf = str(sys.argv[1])[:-7]
ts = int(sys.argv[2])
te = int(sys.argv[3])
nf = int(sys.argv[4])
fx = bf+".t%d-%d.h5" % (ts, te)

# check sys file 
if not os.path.isfile("%s.t%d.sys"%(bf,ts)):
    print("cannot find %s.t%d.sys" % (bf, ts))
    sys.exit()

# read sys file 
fr = open(bf+".t%d.sys"%(ts), 'r') 
for lt in fr.readlines(): 
    if '{box_l' in lt:
        lt = lt.strip() 
        lt = lt.split(' ') 
        dimx = float(lt[1]) 
        dimy = float(lt[2]) 
        dimz = float(lt[3][:-1])
    if 'n_part ' in lt: 
        lt= lt.strip() 
        lt= lt.split(' ') 
        num= int(lt[1][:-1]) 
fr.close() 
dims = array([dimx, dimy, dimz])

# read blockfile output by ESPResSo 
hf = h5py.File(fx, 'w')
trj= hf.create_dataset('trj', (nf, num, 3), dtype='f')
#
fdx = 0 
for t in range(ts, te+1): 
    # 
    with open(bf+".t%d.equ"%(t)) as fr: 
        # 
        while True: 
            chunk = list(islice(fr, num+2)) 
            if not chunk: 
                break   
            # 
            chunk = chunk[1:-1] 
            if not len(chunk) == num: 
                print("insistent particle number.") 
                sys.exit() 
            #
            chunk = char.strip(chunk)
            chunk = char.strip(chunk, "{")
            chunk = char.strip(chunk, "}")
            chunk = char.split(chunk)
            l2as = array( map(l2f, chunk) )
            uxyz = l2as[:, 1:4]
            #
            if fdx == 0:
                hf.create_dataset('box', data=dims, dtype='f') # simulation box
                tys = l2as[:,4].astype(int)
                hf.create_dataset('tys', data=tys, dtype='i4') # particle types
            # save
            trj[fdx,:,:] = uxyz
            # proc
            if fdx%100 == 0: print("t%2d %7d" % (t, fdx))
            # 
            fdx += 1
#
hf.close()

