#!/usr/bin/env python



import os
import sys
import shutil
from glob import *
import pdb
import re


def create_paraview_parameters_file():
    with open(os.path.join(wd,'module_param.f90'),'r') as fp:
        file_text = fp.readlines()
    #strip comments
    for n,l in enumerate(file_text):
        file_text[n] = l[:l.find('!')]
    with open(os.path.join(wd,'param_file.txt'),'w') as fp:
        for n,l in enumerate(file_text):
            if 'integer,parameter :: nx=' in l:
                for prm in re.findall(r'\d+', l):
                    fp.write(prm+'\n')
    with open(os.path.join(wd,'incompact3d.prm'),'r') as fp:
        file_text = fp.readlines()
    with open(os.path.join(wd,'param_file.txt'),'a') as fp:
        for n,l in enumerate(file_text):
            if '#xlx' in l:
                fp.write(l[:l.find('#')].strip()+'\n')
            if '#yly' in l:
                fp.write(l[:l.find('#')].strip()+'\n')
            if '#zlz' in l:
                fp.write(l[:l.find('#')].strip()+'\n')
            if '#nclx' in l:
                fp.write(l[:l.find('#')].strip()+'\n')
            if '#ncly' in l:
                fp.write(l[:l.find('#')].strip()+'\n')
            if '#nclz' in l:
                fp.write(l[:l.find('#')].strip()+'\n')
            if '#istret' in l:
                fp.write(l[:l.find('#')].strip()+'\n')
            if '#iscalar' in l:
                fp.write(l[:l.find('#')].strip()+'\n')

                
                    
            
            

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
#        print (nme, run_dir)
        shutil.copy(nme, './'+run_name+'/code')
    if nme[-4:] in ['.prm']:
#        print (nme, run_dir)
        shutil.copy(nme, './'+run_name)
    if nme.split('/')[-1] == 'incompact3d':
        shutil.copy(nme, './'+run_name)
create_paraview_parameters_file()
os.system('gfortran -o paraview_incompact3d paraview_incompact3d.f90')
shutil.copy(os.path.join(wd,'paraview_incompact3d'), './'+run_name)
shutil.copy(os.path.join(wd,'param_file.txt'), './'+run_name)
os.chdir(run_dir)
#sys.exit()

os.system('mpirun -np 4 ./incompact3d')