#
# Copyright (C) 2019 Lei Liu & Changbong Hyeon
# This file is part of HLM-epidom.
#
import os
import sys
import h5py
from numpy import *
from numba import jit

#
@jit(nopython=True)
def caldrgs(xyz, xs, xe):
    ns = xe - xs + 1
    rgs= zeros((ns, 2))
    #
    for i in range(xs, xe+1):
        for j in range(i+1, xe+1):
            s = j - i
            r2= 0.0
            for k in range(0, 3):
                r2 += var( xyz[i:(j+1), k] )
            #
            rgs[s,0] += sqrt(r2)
            rgs[s,1] += r2
    #
    for k in range(0, 2):
        rgs[:,k] = rgs[:,k]/arange(ns,0,-1)
    return rgs

# Main Func. #################################
# intra-chain r_{g}(s) in experimentally measured regions
if not len(sys.argv) == 6:
    print('usage:: python 2016zhuang418.rgs.py xxx.h5 A-23/I-14/R-11 fs fe fq')
    sys.exit()
#
fx = str(sys.argv[1])
dq = ['A-23','I-14','R-11'].index(str(sys.argv[2]))
fs = int(sys.argv[3])
fe = int(sys.argv[4])
fq = int(sys.argv[5])
rs = 5000.0 # 5 kb per monomer
gs = array([19600000,15700000,2450000]).astype(int)[dq]
ge = gs + 500000

# readin genomic regions measured in exp.
if True:
    if dq == 0:
        ds = ['A-23']
        gz = array([[19726615,20092780]]).astype(int)
    elif dq==1:
        ds = ['I-14']
        gz = array([[15700000,16128463]]).astype(int) # missing: 15602615-15700000
    else:
        ds = ['R-11']
        gz = array([[2487143,2889707]]).astype(int)
    #
    nz = len(gz)
    bz = zeros((nz,2)).astype(int) # the range of monomer index
    nb = zeros(nz).astype(int) # the length of domain in units of monomers
    #
    for z in range(0, nz):
        g = array(gz[z])
        bz[z] = map(int, (g-gs)/rs)
        nb[z] = bz[z,1] - bz[z,0] + 1

# system paras
f5 = h5py.File(fx, 'r')
fmx= f5.get('trj').shape[0]
fi = range(fs, min(fmx, fe), fq)
nf = len(fi)
N  = len(f5.get('tys'))
f5.close()

# r_{g}(s)
rgs= [zeros((nb[z], 2)) for z in range(0, nz)]
fr = h5py.File(fx, 'r')
#
for f in range(0, nf):
    ups = fr.get('trj')[fi[f]]
    for z in range(0, nz):
        rgs[z] += caldrgs(ups, bz[z,0], bz[z,1])
fr.close()
for z in range(0, nz):
    rgs[z] /= float(nf)

#
fw = open(fx[:-3]+'.2016zhuang418.rgs', 'w')
for z in range(0, nz):
    ct = "#%s\n" % (ds[z])
    ct+= "#%8d %8d\n" % (gz[z,0], gz[z,1])
    ct+= "#%8d %8d\n" % (bz[z,0]*5000+gs, (bz[z,1]+1)*5000+gs)
    ct+= "#%2d %2d\n" % (bz[z,0], bz[z,1])
    ct+= '#s[5kb] r_{g}(s) std\n'
    fw.write(ct)
    #
    for s in range(0, nb[z]):
        lt  = "%3d" % (s)
        lt += "%11.4e %11.4e " % ( rgs[z][s,0], sqrt(rgs[z][s,1] - rgs[z][s,0]**2) )
        fw.write(lt+'\n')
    fw.write('\n\n')
fw.close()

