# !/usr/bin/env python

# written on 9 Apr 2026 by Charitra Jain
# updated on 10 Sep 2026
# to plot viscosity profiles at different times
# for models presented in the manuscript
# Jain & Sobolev (2026): Controls on the tectonic evolution of the early Earth: 1. Composite rheology
# submitted to JGR: Solid Earth

# This script uses StagPy library to post-process StagYY data.
# Please cite as:
# Adrien Morison, Martina Ulvrova, Stephane Labrosse, B4rsh, theofatou, & tfrass49. (2024).
# StagPython/StagPy: v0.20.1 (v0.20.1). Zenodo. https://doi.org/10.5281/zenodo.13684760

# Other python packages required: os, numpy, matplotlib

# Once required packages are installed, execute the script from terminal with command "python plot_viscosity_profiles.py"

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from stagpy import stagyydata


### beginning of input ###
# user only has to input figure number (2/6/9) as batch
batch=2
### end of input ###


if batch == 2:
   model_folder = 'Table4_models/'
   model_list = ["02_QDf01r2/", "10_QCf01r3/", "15_QPf01r1/"]
   profile_frame_list=[80,255,281,447,\
                    92,350,467,781,\
                    71,255,296,505]
   profile_time_list=[85,652,670,1250,\
                   85,480,705,1250,\
                   85,526,556,1244]
elif batch == 6:
   model_folder = 'Table6_models/'
   model_list = ["01_NoMeltQDf005r1/", "08_NoMeltQCf005r2/", "13_NoMeltQPf005r1/"]   
   profile_frame_list=[18,110,707,715,\
                    25,1857,1864,1873,\
                    19,140,425,432]
   profile_time_list=[49,300,1296,1299,\
                   48,1171,1172,1173,\
                   50,302,808,809]
elif batch == 9:
   model_folder = 'Table7_models/'   
   model_list = ["11_QPf005r2/", "14_QPf01r5/", "17_QPf02r2/"]
   profile_frame_list=[120,296,466,613,\
                    87,287,420,642,\
                    64,337,429,688]
   profile_time_list=[218,1405,3044,4500,\
                   235,1335,2454,4500,\
                   210,1481,2111,4500]
else:
   print('Wrong figure')

here = os.getcwd()

if os.path.isdir(here+'/../figures/'):
   print('Directory for figures exists')
else:
   print('Making directory for figures')
   os.mkdir(here+'/../figures')

total_runs = len(model_list)
print('# models being processed: '+str(total_runs))


# customising matplotlib and figure
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.sans-serif'] = "Helvetica"
mpl.rcParams['font.size'] = 8
mpl.rcParams['axes.titlesize'] = 8
mpl.rcParams['axes.titleweight'] = "bold"
mpl.rcParams['axes.labelweight'] = "normal"
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['legend.fontsize'] = 6
mpl.rcParams['figure.dpi'] = 600
mpl.rcParams['figure.figsize'] = [1.5,4]
mpl.rcParams['xtick.direction'] = "in"
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.direction'] = "in"
mpl.rcParams['ytick.minor.visible'] = False
mpl.rcParams["pdf.use14corefonts"] = True
box_ratio=0.3


seconds_in_year = 31556926
here = os.getcwd()

nz=128

run=int(0)
for path in model_list:
   
   current_model=path[0:len(path)-1]

   # from stagpy import stagyydata
   sdat = stagyydata.StagyyData(here+'/../'+str(model_folder)+path)   
   
   # reading in the timesteps when snapshots are saved
   snaps = []
   file = (here+'/../'+str(model_folder)+path+'+op/earthDIM_rprof.dat')
   f = open(file,'r')
   lines=f.readlines()[1:]
   f.close
   nlines=int(len(lines))
   for i in np.arange(0,nlines,nz+1):
      p = lines[i].split()
      snaps.append(int(p[1]))
   nframes = int(len(snaps))
   print('Number of frames: ' + str(nframes))
   # timesteps for initial and final states
   step=snaps[1]
   timef=sdat.tseries.time[len(sdat.tseries.time)-1]*1e-6/seconds_in_year
   print('Time at last frame: '+str('{:.2f}'.format(timef)))

   visc1=sdat.steps[snaps[profile_frame_list[(run*4)]]].rprofs['etalog'].values
   visc2=sdat.steps[snaps[profile_frame_list[(run*4)+1]]].rprofs['etalog'].values
   visc3=sdat.steps[snaps[profile_frame_list[(run*4)+2]]].rprofs['etalog'].values
   visc4=sdat.steps[snaps[profile_frame_list[(run*4)+3]]].rprofs['etalog'].values

   radius=sdat.steps[step].rprofs['r'].values
   depth_CMB_to_surface = abs(2890.e3 - radius)

   fig, ax = plt.subplots(1,1)
   ax.plot(visc1,depth_CMB_to_surface/1e3,color = '#C25B60',linestyle='-',alpha=1)
   ax.plot(visc2,depth_CMB_to_surface/1e3,color = '#EE9B00',linestyle='-',alpha=1)
   ax.plot(visc3,depth_CMB_to_surface/1e3,color = '#43B7C2',linestyle='-',alpha=1)
   ax.plot(visc4,depth_CMB_to_surface/1e3,color = '#00647A',linestyle='-',alpha=1)
   ax.set_xlim(1e18, 1e24)
   ax.set_xscale('log')
   ax.set_xticks([1e18,1e20,1e22,1e24])
   ax.set_xlabel('Viscosity (Pa s)')
   ax.set_ylim(660.,0.)
   ax.set_yscale('linear')
   ax.set_ylabel('Depth (km)')

   ax.set_box_aspect(1/box_ratio)
   ax.grid(color='k', linestyle='--', linewidth=0.3,axis='both',alpha=0.3)
   labels_run = [str(profile_time_list[(run*4)])+' Myr',str(profile_time_list[(run*4)+1])+' Myr',str(profile_time_list[(run*4)+2])+' Myr',str(profile_time_list[(run*4)+3])+' Myr']
   print(labels_run)
   location='center right'
   if current_model == '11_QPf005r2':
       location='lower left'
   ax.legend(ncol = 1, labels=labels_run, loc=location, \
   labelspacing = 0.5, edgecolor = 'none', facecolor = 'none', fancybox = False)
   
   fig.tight_layout()
   file_name='fig'+str(batch)+'_viscosity_'+current_model+'.pdf'
   plt.savefig(here+'/../figures/'+file_name,transparent=True)

   run = run + 1
