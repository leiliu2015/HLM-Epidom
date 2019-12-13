#
# Copyright (C) 2019 Lei Liu & Changbong Hyeon
# This file is part of HLM-epidom.
#
import sys
import h5py
from numpy import *
from numba import jit

#
@jit(nopython=True)
def calCM(ups, N, rc):
    cm = zeros((N, N))
    #
    for i in range(0, N):
        cm[i,i] = 1.0
        for j in range(i+1, N):
            r2 = 0.0
            for k in range(0, 3):
                r2 += (ups[j,k] - ups[i,k])**2
            #
            rij= sqrt(r2)
            #
            if rij <= rc:
                cm[i,j] = 1.0
                cm[j,i] = 1.0
    #
    return cm

# Main Func. #################################
# calculate {p_ij} based on structures
if not len(sys.argv) == 6:
    print('usage:: python gnm.h5cm.py xxx.h5 Rc fs fe fq')
    sys.exit()
#
fx = str(sys.argv[1])
rc=float(sys.argv[2]) # r_{contact}, cutoff
fs = int(sys.argv[3])
fe = int(sys.argv[4])
fq = int(sys.argv[5])

# system paras
f5 = h5py.File(fx, 'r')
fmx= f5.get('trj').shape[0]
fi = range(fs, min(fmx, fe), fq)
nf = len(fi)
N  = f5.get('tys').shape[0]
f5.close()

#
cm = zeros((N, N))
fr = h5py.File(fx, 'r')
for f in range(0, nf):
    ups = fr.get('trj')[fi[f]]
    cm += calCM(ups, N, rc)
fr.close()

# frequency -> probability
cm /= float(nf)

# save {p_ij}
fw = open(fx[:-3]+".rc%3.1f.cm"%(rc), 'w')
fw.write("#shape: %d\n"%(N))
#
for i in range(0, N):
    lt = ''
    for j in range(0, N):
        if isnan(cm[i,j]) or cm[i,j]==0:
            lt += "%11s " % (NaN)
        else:
            lt += "%11.5e " % (cm[i,j])
    fw.write(lt+'\n')
#
fw.close()

