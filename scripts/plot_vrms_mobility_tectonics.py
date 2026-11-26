# !/usr/bin/env python

# written on 9 Apr 2026 by Charitra Jain
# updated on 10 Sep 2026
# added legends on 26 Nov 2026
# to plot rms surface velocities, mobility and inferred tectonic regimes
# for models presented in the manuscript
# Jain & Sobolev (2026): Controls on the tectonic evolution of the early Earth: 1. Composite rheology
# submitted to JGR: Solid Earth

# This script uses StagPy library to post-process StagYY data.
# Please cite as:
# Adrien Morison, Martina Ulvrova, Stephane Labrosse, B4rsh, theofatou, & tfrass49. (2024).
# StagPython/StagPy: v0.20.1 (v0.20.1). Zenodo. https://doi.org/10.5281/zenodo.13684760

# Other python packages required: os, math, numpy, h5py, pandas, matplotlib

# Once required packages are installed, execute the script from terminal with command "python plot_vrms_mobility_tectonics.py"

import os
import math
import numpy as np
import h5py
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from stagpy import stagyydata


### beginning of input ###
# user only has to input figure number (2/6/9) as batch
batch=2
### end of input ###


if batch == 2:
   model_folder = 'Table4_models/'
   model_list = ["02_QDf01r2/", "10_QCf01r3/", "15_QPf01r1/"]
   evolved_time = 1500.
   time_ticks=[0.,300.,600.,900.,1200.,1500.]
elif batch == 6:
   model_folder = 'Table6_models/'
   model_list = ["01_NoMeltQDf005r1/", "08_NoMeltQCf005r2/", "13_NoMeltQPf005r1/"]   
   evolved_time = 1500.
   time_ticks=[0.,300.,600.,900.,1200.,1500.]   
elif batch == 9:
   model_folder = 'Table7_models/'   
   model_list = ["11_QPf005r2/", "14_QPf01r5/", "17_QPf02r2/"]
   evolved_time = 4500.
   time_ticks=[0,500,1000,1500,2000,2500,3000,3500,4000,4500]   
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
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['figure.dpi'] = 600
mpl.rcParams['figure.figsize'] = [6.25,1] # (1 panel, width x height, 1=600px with dpi=300)
mpl.rcParams['xtick.direction'] = "in"
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.direction'] = "in"
mpl.rcParams['ytick.minor.visible'] = False
mpl.rcParams["pdf.use14corefonts"] = True

box_ratio=12
max_mobility = 3.0
mobility_ticks=[0.0,1.0,2.0,3.0]
min_vrms = 1e-3
max_vrms = 1e3
vrms_ticks=[min_vrms,1,max_vrms]

ML=['#024b79']
PSL=['#43b7c2']
RO=['#ffad48']
EL=['#c25b60']
SL=['#822b42']

velocity_rms_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
time_Myr_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
mobility_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
steps_mobility_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
mobility_time_avg = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
time_Myr_mobility_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
moving_average_mobility_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
plateness_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
moving_average_plateness_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
plateness_time_avg = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))

seconds_in_year = 31556926
seconds_in_Myear = seconds_in_year*1e6
earth_radius = 6371.e3
core_radius = 3481.e3
nz = 128

# define size of moving window for rolling average
window_size_mobility  = 10
window_size_plateness = 10

run = int(0)
for path in model_list:

   volume_mantle = (4 * math.pi*(earth_radius**3-core_radius**3))/3
   # quadrant analytical volume is 1/16 of full annulus, which in itself is pi/ny times actual volume

   nytot = 512
   analytical_volume_quad = (4 * math.pi*(earth_radius**3-core_radius**3))/3 * (math.pi)/(nytot*16)
   factor=volume_mantle/analytical_volume_quad
   if run == 0:
      print('Scaling factor for quadrant annulus: '+str(factor))

   print('Run # '+str(run+1))
   print('Model  '+path)

   # from stagpy import stagyydata
   sdat = stagyydata.StagyyData(here+'/../'+str(model_folder)+path)

   # extracting time from time_botT.h5
   # this is needed when plotting the time evolution of a field
   time_list = []
   f = h5py.File(here+'/../'+str(model_folder)+path+'+hdf5/time_botT.h5', 'r')

   ndatasets = len(f.keys())
   for key in f.keys():
      data = f[key]
      time_list.append(float(data[0]))
   f.close()

   time=np.array(time_list)
   time_Myr=time/seconds_in_Myear
   print('# HDF datasets: '+str(ndatasets))
   time_Myr_run[run]=time_Myr
   time_list.clear()



   # reading in the timesteps when snapshots are saved
   snaps = []
   file = (here+'/../'+str(model_folder)+path+'snaps.dat')
   f = open(file,'r')
   lines = f.readlines()
   f.close
   for line in lines:
      p = line.split()
      snaps.append(int(p[0]))
   nframes = int(len(snaps))
   print('# frames (rprof.dat): ' + str(nframes))
   print('Evolution time (Myr):    '+str('{:.2f}'.format(sdat.tseries.time[len(sdat.tseries.time)-1]*1e-6/seconds_in_year)))
   # timesteps for initial and final states
   step1=snaps[1]
   step2=snaps[nframes-1]


   # Using pandas to calculate moving average (with frames) of mobility and plateness
   # models presented in figures 3 and 4 have extra columns in their files.
   file = (here+'/../'+str(model_folder)+path+'+op/earthDIM_plates_analyse.dat')
   df = pd.read_table(file,sep='  ',engine='python', header=None, skiprows=1,usecols=range(6))                                                       
   df = df.set_axis(['istep','time', 'Rtp_vol', 'Rtp_surf', 'mobility', 'f80'], axis=1)

   df['mobility_10'] = df.mobility.rolling(window_size_mobility, min_periods=1).mean()
   moving_average_mobility=np.array(df['mobility_10'])
   moving_average_mobility_run[run] = moving_average_mobility
   df['plateness'] = 1.-(df['f80']/0.6)
   df['plateness_10'] = df.plateness.rolling(window_size_plateness, min_periods=1).mean()
   moving_average_plateness=np.array(df['plateness_10'])
   moving_average_plateness_run[run] = moving_average_plateness


   # loading mobility and plateness from _plates_analyse.dat
   mobility_list = []
   f80_list = []
   plateness_list = []
   time_list = []
   steps_mobility_list = []
   file = (here+'/../'+str(model_folder)+path+'+op/earthDIM_plates_analyse.dat')
   f = open(file,'r')
   lines=f.readlines()[1:]
   f.close
   nlines=int(len(lines))

   for i in np.arange(0,nlines):
      p = lines[i].split()
      steps_mobility_list.append(int(p[0]))
      time_list.append(float(p[1]))
      mobility_list.append(float(p[4]))
      f80_list.append(float(p[5]))
      plateness_list.append(1.-(float(p[5])/0.6))

   # calculating time-averaged mobility
   sum=0.
   for i in np.arange(0,nlines-1):
      sum = sum + (time_list[i+1]-time_list[i]) * mobility_list[i]

   mobility_time_avg[run] = sum/time_list[nlines-1]
   print('Time-averaged mobility:  '+str('{:.2f}'.format(mobility_time_avg[run])))

   # calculating time-averaged plateness
   sum=0.
   for i in np.arange(0,nlines-1):
      sum = sum + (time_list[i+1]-time_list[i]) * plateness_list[i]

   plateness_time_avg[run] = sum/time_list[nlines-1]
   # print('Time-averaged plateness:      '+str('{:.2f}'.format(plateness_time_avg[run])))

   mobility_run[run] = np.array(mobility_list)
   plateness_run[run] = np.array(plateness_list)
   steps_mobility_run[run] = np.array(steps_mobility_list)
   time_Myr_mobility_run[run]=np.array(time_list)/seconds_in_Myear

   mobility_list.clear()
   f80_list.clear()
   plateness_list.clear()
   steps_mobility_list.clear()
   time_list.clear()



   # loading velocities (in the top cell, nz-1) from radial profiles
   velocity_rms_list = [0] # this is vsurf_rms
   for i in np.arange(1,nframes):
      velocity_rms_list.append(float(sdat.steps[snaps[i]].rprofs['vrms'].values[nz-1]*100*seconds_in_year))

   velocity_rms_run[run] = np.array(velocity_rms_list)
   velocity_rms_list.clear()


   # plotting
   fig, ax = plt.subplots(1, 1, squeeze=False)
   ax[0][0].plot(time_Myr_run[run],velocity_rms_run[run], linestyle='-', linewidth=1, alpha=1, color = '#000000')
   
   current_model=path[0:len(path)-1]

   if current_model == '02_QDf01r2':
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,150.,    linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],150.,205.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],205.,215.,    linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],215.,320.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],320.,350.,    linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],350.,655.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],655.,680.,    linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],680.,765.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],765.,800.,    linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],800.,840.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],840.,860.,    linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],860.,890.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],890.,1500., linewidth=0,color = SL,alpha=1,edgecolor='none')
      ax[0][0].axvline(85., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(652., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(670., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(1250., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].set_xlim(0., evolved_time)
      loc1='upper right'
      loc2='lower right'
   elif current_model == '10_QCf01r3':
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,290.,     linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],290.,465.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],465.,530.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],530.,810.,   linewidth=0,color = SL,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],810.,1150.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1150.,1200.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1200.,1280.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1280.,1320.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1320.,1460.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1460.,1500.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].axvline(85., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(480., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(705., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(1250., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      loc1='upper center'
      loc2='lower center'
   elif current_model == '15_QPf01r1':
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,310.,     linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],310.,530.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],530.,620.,     linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],620.,1220.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1220.,1320., linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1320.,1375.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1375.,1500., linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].axvline(85., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(526., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(556., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(1244., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      loc1='upper center'
      loc2='lower center'
   elif current_model == '01_NoMeltQDf005r1':
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,90.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],90.,500.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],500.,990.,  linewidth=0,color = SL,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],990.,1025.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1025.,1296.,  linewidth=0,color = SL,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1296.,1340.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1340.,1500.,  linewidth=0,color = SL,alpha=1,edgecolor='none')
      ax[0][0].axvline(49., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(300., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(1296., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(1299., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      loc1='upper center'
      loc2='lower center'
   elif current_model == '08_NoMeltQCf005r2':
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,140.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],140.,480.,  linewidth=0,color = SL,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],480.,1500.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].axvline(48., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(1171., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(1172., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(1173., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      loc1='center'
      loc2='center right'
   elif current_model == '13_NoMeltQPf005r1':
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,200.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],200.,380.,  linewidth=0,color = SL,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],380.,600.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],600.,810.,  linewidth=0,color = SL,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],810.,1500.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].axvline(50., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(302., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(808., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(809., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      loc1='lower center'
      loc2='lower right'
   elif current_model == '11_QPf005r2':
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,4500.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].axvline(218., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(1405., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(3044., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      loc1='upper center'
      loc2='upper right'
   elif current_model == '14_QPf01r5':
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,275.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],275.,400.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],400.,750.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],750.,1200.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1200.,1400.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1400.,1800.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1800.,2150.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],2150.,2400.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],2400.,2500.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],2500.,2700.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],2700.,3250.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3250.,3500.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3500.,4500.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].axvline(235., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(1335., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(2454., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      loc1='upper center'
      loc2='upper right'
   elif current_model == '17_QPf02r2':
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,350.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],350.,480.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],480.,1470.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1470.,1500.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1500.,1600.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1600.,1690.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1690.,2100.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],2100.,2150.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],2150.,2200.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],2200.,2280.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],2280.,2930.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],2930.,3030.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3030.,3250.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3250.,3500.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3500.,4060.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],4060.,4180.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],4180.,4500.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].axvline(210., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(1481., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(2111., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      loc1='upper right'
      loc2='lower right'

   ax[0][0].set_xlim(0., evolved_time)
   ax[0][0].set_xscale('linear')
   ax[0][0].set_xticks(time_ticks)      
   ax[0][0].set_ylim(min_vrms,max_vrms)
   ax[0][0].set_yscale('log')
   ax[0][0].set_yticks(vrms_ticks)
   ax[0][0].set_ylabel('$v_{rms}$ (cm/yr)',labelpad=4)
   ax[0][0].set_xlabel('Time (Myr)')
   ax[0][0].set_box_aspect(1/box_ratio)
   ax[0][0].legend(ncol = 1, labels=['v$_{rms}$'], loc=loc1, \
   labelspacing = 0.5, edgecolor = 'none', facecolor = 'none', fancybox = False)

   ax2 = ax[0][0].twinx()
   ax2.plot(time_Myr_run[run],moving_average_mobility_run[run], linestyle='-', linewidth=1, alpha=0.7, color = '#ffffff')
   ax2.set_xlim(0., evolved_time)
   ax2.set_xscale('linear')
   ax2.set_xticks(time_ticks)      
   ax2.set_yscale('linear')
   ax2.set_ylim(0.,max_mobility)
   ax2.set_yticks(mobility_ticks)
   ax2.set_ylabel('Mobility',labelpad=4)
   ax2.set_box_aspect(1/box_ratio)    
   ax2.legend(ncol = 1, labels=['mobility'], loc=loc2, \
   labelspacing = 0.5, labelcolor = '#ffffff', edgecolor = 'none', facecolor = 'none', fancybox = False)  

   fig.tight_layout()

   file_name='fig'+str(batch)+'_vrms_'+current_model+'.pdf'
   plt.savefig(here+'/../figures/'+file_name,transparent=True)

   run = run + 1
