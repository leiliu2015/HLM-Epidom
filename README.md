## HLM-Epidom

As an application of [HLM](https://github.com/leiliu2015/HLM), we used this polymer-based approach to generate three-dimensional (3D) structures of different epigenetic domains in *Drosophila*, and compared their structural properties. This repository provides necessary files to reproduce most results shown in our recent [study]() (e.g., Fig. 2A-C, Fig. 3B and C), which includes mainly
1. Input files of Molecular Dynamics (MD) simulations to generate 3D structures.
2. Python scripts to analyze the structures.

### System Requirements
The code was tested on ubuntu 14.04/16.04 LTS. 3D structures are generated by running MD simulations with [ESPResSo 3.3.1 package](http://espressomd.org/wordpress/) (***optional***). We recommand [Anaconda](https://www.anaconda.com/distribution/) to manage the Python environment (Python 2.7) and other packages (H5py, NumPy, Numba) required for the subsequent analysis. The final results are visualized with [Gnuplot](gnuplot.sourceforge.net).  Please refer to their websites for the instructions of installation. 

In addition, we used Voronoi tessellation (***optional***) to define the volume and surface of the space occupied by chromatin chain. To apply this method in the analysis, please first download and install [Voro++ 0.4.6](math.lbl.gov/voro++/) on your computer. Then, compile the C++ file in the directory [voro](voro/) by typing the following line at the root of the repository in your terminal, 
```
g++ -Wall -ansi -pedantic -O3 -I VoroPath/src -L VoroPath/src -o voro/irregular voro/irregular.cc -lvoro++
```
where `VoroPath` should be the path to the directory where Voro++ 0.4.6 has been installed. 

### File Description
- Hi-C/
  - Kc167-5kb.chrxxxxx.cm (Contact frequency matrices, {*p<sub>ij</sub>*}, measured by [Eagen *et al.*](https://www.pnas.org/content/114/33/8764) with Hi-C, based on which we inferred the values of model parameters {*k<sub>ij</sub>*} and $\chi$. See also Table S1 in our [preprint]() for their genomic positions.)
 - A-23/ I-14/ R-11/ Fig-S6/ Fig-S7/ Fig-S8/ Fig-S9/
   - MD/
     - [hlm_espresso.tcl](A-23/MD/hlm_espresso.tcl) (A TCL script used as the input of ESPResSo)
     - [hlm.iniCfg](A-23/MD/hlm.iniCfg) (An initial input configuration read by the TCL script)
     - [hlm.km](A-23/MD/hlm.km) ({*k<sub>ij</sub>*} derived from Hi-C, which will be read by the TCL script)
   - gnm.t0-0.h5 (A HDF5 file which stores 10<sup>4</sup> configurations generated by MD simulations)
- analyze-scripts/
  - [gnm.bk2h5.py](analyze-scripts/gnm.bk2h5.py) (A Python script to change the coordinates format to HDF5)
  - [gnm.h5cm.py](analyze-scripts/gnm.h5cm.py) (A Python script to calculate {*p<sub>ij</sub>*} based on structures)
  - [2016zhuang418.rgs.py](analyze-scripts/2016zhuang418.rgs.py) (A Python script to calculate *r<sub>g</sub>(s)* in A-23, I-14, and R-11 domains based on structures)
  - [2016zhuang418.vor.py](analyze-scripts/2016zhuang418.vor.py) (A Python script to calculate chromatin density, surface roughness, and asphericity of A-23, I-14, and R-11 domains based on structures)
- gnuplot-scripts/
  - [pij.gnu](gnuplot-scripts/pij.gnu) (A Gnuplot script to compare {*p<sub>ij</sub>*} based on Hi-C and HLM)
  - [rgs.gnu](gnuplot-scripts/rgs.gnu) (A Gnuplot script to compare *r<sub>g</sub>(s)* in different domains)
  - [vor.gnu](gnuplot-scripts/vor.gnu) (A Gnuplot script to compare chromatin density and other properties of different domains)
- voro/
  - [irregular.cc](voro/irregular.cc) (A C++ file which tessellates the space occupied by the domain with dodecahedrons)

### User Guide
> ***Note:*** Since we have provided the structures of A-23, I-14, and R-11 domain, if you don't want to spend hours in MD simulations, you may skip the first part. In that case, the whole analysis takes less than 10 minutes.

**F**irst, we perform MD simulations using ESPResSo to generate an ensemble of 3D structures, which usually takes a few hours.
```
$ cd ./A-23/MD
$ Espresso hlm_espresso.tcl
$ python ../../analyze-scripts/gnm.bk2h5.py gnm.t0.equ 0 0 10000
$ mv ./gnm.t0-0.h5 ../
```
The coordinates of the chromatin chain will be stored in a text file called `gnm.t0.equ`. To facilitate the following analysis, we change the coordinates format to a binary HDF5 format with the Python script named `gnm.bk2h5.py` , and move the output to a upper level directory.

**N**ext, we calculate the contact probability matrix {*p<sub>ij</sub>*}, the mean radius of gyration of subchains *r<sub>g</sub>(s)*, the asphericity *Asp.*, and other properties of the domain based on the structures.
```
$ cd ./A-23
$ python ../analyze-scripts/gnm.h5cm.py gnm.t0-0.h5 2.2 0 10000 2
$ python ../analyze-scripts/2016zhuang418.rgs.py gnm.t0-0.h5 A-23 0 10000 10
$ python ../analyze-scripts/2016zhuang418.vor.py gnm.t0-0.h5 A-23 0 10000 10
```
The modeling and analyzing of other two domains follow similar steps, with `A-23` being replaced by `I-14` and `R-11`. After finishing the analysis for those three domains, type the following lines at the root of the repository to compare the results.
```
$ cd ./gnuplot-scripts
$ gnuplot -persist pij.gnu
$ gnuplot -persist rgs.gnu
$ gnuplot -persist vor.gnu
```
These will reproduce Fig. 2A-C, Fig. 2B and C in our [preprint](), respectively. We have deposited input files for MD simulations in the directories Fig-Sx/MD as well. However, since each genomic region contains more than one epigenetic domains (see Table S1 and S2 in the [preprint]()), analysis on the structures needs additional care. All the output files can be deleted by typing `$ bash ./clearAll.sh` at the repository root. If you are interested in other possible applications of HLM, or have further questions about it, please contact Lei Liu (leiliu2015@163.com).

