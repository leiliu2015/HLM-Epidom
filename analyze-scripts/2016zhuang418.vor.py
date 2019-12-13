#
# Copyright (C) 2019 Lei Liu & Changbong Hyeon
# This file is part of HLM-epidom.
#
import os
import sys
import h5py
from numpy import *
from numba import jit

# asphericity: sphere 0; nonsphere (0,1]
@jit
def calAsphericity(xyz):
    # centre of domain
    dc = mean(xyz, axis=0)
    xyz= xyz - dc
    # gyratio tensor
    gt = zeros((3, 3))
    for i in range(3):
        for j in range(3):
            gt[i,j] = mean(xyz[:,i]*xyz[:,j])
    #
    vs = linalg.eigvals(gt)
    mv = mean(vs)
    #
    asp = 1.5*sum((vs-mv)**2)/(sum(vs)**2)
    #
    return asp

# R_{g}
@jit
def calRg(xyz):
    r2 = 0.0
    for k in range(0, 3):
        r2 += var(xyz[:,k])
    rg = sqrt(r2)
    #
    return rg

# Main Func. #################################
# Density, Surface Roughness, Asphericity, L/R_{g}^{3} of experimentally measured regions
if not len(sys.argv) == 6:
    print('usage:: python 2016zhuang418.vor.py xxx.h5 A-23/I-14/R-11 fs fe fq')
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

# default box size and particle size were set in L11 & L24 of 'irregular.cc', respectively.
dv = '../voro/'
bx = 50.0
cb = ones(3)*bx/2.0
vq = 1 if os.path.isfile(dv+'irregular') else 0 # check whether Voro++ is installed.

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

#
cc = (36.0*pi)**(1.0/3)
vds= zeros((nf, nz)) # volume of domain
sds= zeros((nf, nz)) # surface-roughness of domain
ads= zeros((nf, nz)) # aspherisity of domain
rds= zeros((nf, nz)) # r_g of domain
fr = h5py.File(fx, 'r')
#
for f in range(0, nf):
    xyz = fr.get('trj')[fi[f]]
    #
    if vq:
        # -- prepare xyz for voro
        cog= mean(xyz, axis=0)
        xyz= xyz + (cb-cog)
        fw = open('test.xyz', 'w')
        for i in range(0, len(xyz)):
            lt = "%d %.3f %.3f %.3f" % (i+1, xyz[i,0], xyz[i,1], xyz[i,2])
            fw.write(lt+'\n')
        fw.close()
        # -- run voro
        os.system(dv+'irregular')
        # -- read voro outputs
        vf = zeros(N)
        sf = zeros(N) # exposed surface area (with vaccum)
        sm = zeros((N, N)) # surface between neighboring monomers 
        fv = open('test.vol', 'r')
        for line in fv.readlines():
            lt = line.strip()
            lt = lt.split()
            ip = int(lt[0])-1
            vi=float(lt[1])
            #
            vf[ip] = vi
            #
            ni = int(lt[2])
            if not len(lt) == (2*ni+3):
                print('check'+line)
                sys.exit()
            else:
                js = map(int, lt[-ni:])
                ss = map(float, lt[3:(3+ni)])
                #
                for k in range(0, ni):
                    j = js[k]
                    s = ss[k]
                    #
                    if j>0:
                        sm[ip, j-1] = s
                        sm[j-1, ip] = s
                    else:
                        sf[ip] += s
        fv.close()
        #
        for z in range(0, nz):
            bs = bz[z,0]
            be = bz[z,1]+1
            vd = vf[bs:be].sum() # volume of domain
            sd = sf[bs:be].sum() + sm[bs:be,:].sum() - sm[bs:be, bs:be].sum() # surface of domain
            #
            vds[f,z] = vd
            sds[f,z] = sd / vd**(2.0/3) / cc
        # -- remove voro outputs
        os.remove('test.xyz')
        os.remove('test.vol')
        os.remove('test.gnu')
    #
    for z in range(0, nz):
        bs = bz[z,0]
        be = bz[z,1]+1
        ads[f,z] = calAsphericity(xyz[bs:be])
        rds[f,z] = calRg(xyz[bs:be])
fr.close()


# -- save density of domains
if vq:
    fw = open(fx[:-3]+'.2016zhuang418.rho', 'w')
    for z in range(0, nz):
        ys = float(nb[z]) / vds[:,z] # number density
        lt = "%s %s [ %3d %3d ] " % (ds[z],ds[z][0],bz[z,0],bz[z,1])
        lt+= "rho: %11.5e %11.5e qt/ " % (mean(ys), std(ys))
        for q in range(0, 101, 25):
            lt += "%11.5e " % (percentile(ys, q))
        fw.write(lt+'\n')
    fw.close()

# -- save surface roughness of domains
if vq:
    fw = open(fx[:-3]+'.2016zhuang418.srf', 'w')
    for z in range(0, nz):
        ys = sds[:,z]
        lt = "%s %s [ %3d %3d ] " % (ds[z],ds[z][0],bz[z,0],bz[z,1])
        lt+= "S/S_{0}: %11.5e %11.5e qt/ " % (mean(ys), std(ys))
        for q in range(0, 101, 25):
            lt += "%11.5e " % (percentile(ys, q))
        fw.write(lt+'\n')
    fw.close()

# -- save asphericity of domains
fw = open(fx[:-3]+'.2016zhuang418.asp', 'w')
for z in range(0, nz):
    ys = ads[:,z]
    lt = "%s %s [ %3d %3d ] " % (ds[z],ds[z][0],bz[z,0],bz[z,1])
    lt+= "asphericity: %11.5e %11.5e qt/ " % (mean(ys), std(ys))
    for q in range(0, 101, 25):
        lt += "%11.5e " % (percentile(ys, q))
    fw.write(lt+'\n')
fw.close()

# -- save L/R_{g}^{3} of domains
fw = open(fx[:-3]+'.2016zhuang418.irg', 'w')
for z in range(0, nz):
    ys = float(nb[z]) / (rds[:,z]**3)
    lt = "%s %s [ %3d %3d ] " % (ds[z],ds[z][0],bz[z,0],bz[z,1])
    lt+= "L/R_g^3: %11.5e %11.5e qt/ " % (mean(ys), std(ys))
    for q in range(0, 101, 25):
        lt += "%11.5e " % (percentile(ys, q))
    fw.write(lt+'\n')
fw.close()

