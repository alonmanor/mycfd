#!/usr/bin/env python



import os
import sys
import shutil
from glob import *
import pdb

run_name = sys.argv[1]
wd = os.getcwd()
print(wd)
print('run name:',run_name)
run_dir = os.path.join(wd,run_name)
if os.path.exists(run_dir):
    print('run already exists')
    sys.exit()
os.mkdir(run_dir)
os.mkdir(os.path.join(run_dir, 'code'))
v = glob(os.path.join(wd,'*'))
#pdb.set_trace()
for nme in v:
    if nme[-4:] in ['.f90']:
        print (nme, run_dir)
        shutil.copy(nme, './'+run_name+'/code')
    if nme[-4:] in ['.prm']:
        print (nme, run_dir)
        shutil.copy(nme, './'+run_name)
    if nme.split('/')[-1] == 'incompact3d':
        shutil.copy(nme, './'+run_name)
os.system('gfortran -o paraview_incompact3d paraview_incompact3d.f90')
shutil.copy(os.path.join(wd,'paraview_incompact3d'), './'+run_name)
os.chdir(run_dir)
os.system('mpirun -np 4 ./incompact3d')