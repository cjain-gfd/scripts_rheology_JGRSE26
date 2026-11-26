# !/usr/bin/env python

# written on 9 Apr 2026 by Charitra Jain
# updated on 10 Sep 2026
# to make figures for the manuscript
# Jain & Sobolev (2026): Controls on the tectonic evolution of the early Earth: 1. Composite rheology
# submitted to JGR: Solid Earth

# This script uses StagPy library to post-process StagYY data.
# Please cite as:
# Adrien Morison, Martina Ulvrova, Stephane Labrosse, B4rsh, theofatou, & tfrass49. (2024).
# StagPython/StagPy: v0.20.1 (v0.20.1). Zenodo. https://doi.org/10.5281/zenodo.13684760

# Other python packages required: os, math, numpy, h5py, pandas, matplotlib

# Once required packages are installed, execute the script from terminal with command "python make_figures_rheology_NoMelt.py"

import os
import math
import numpy as np
import h5py
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from stagpy import stagyydata


### beginning of input ###
# user only has to input figure number (7 or 8) as batch
batch=7
### end of input ###


if batch == 7:
   model_folder = 'Table6_models/'
   model_list = ["01_NoMeltQDf005r1/", "02_NoMeltQDf005r2/", "03_NoMeltQDf005r3/", \
                 "07_NoMeltQCf005r1/", "08_NoMeltQCf005r2/", "09_NoMeltQCf005r3/", \
                 "13_NoMeltQPf005r1/", "14_NoMeltQPf005r2/", "15_NoMeltQPf005r3/"]
   palette = ['#C25B60','#FFBF70','#43B7C2','#C25B60','#FFBF70','#43B7C2']
   labels_a = ('diffusion (3x)','composite (3x)','composite proxy (3x)','NoMeltQDf005r1','NoMeltQCf005r2','NoMeltQPf005r1')
   figure_title = ('Quadrant annulus | $\\it{T}$$_{P0}$ = 1900 K | $\\mu$ = 0.05 | different rheology without melting')
   location='lower center'
   style = ['solid','solid','solid','dashed','dashed','dashed']
   symbol = ['^','o','s','^','o','s']
   legend_columns = 3
   legend_panel = 2
elif batch == 8:
   model_folder = 'Table6_models/'   
   model_list = ["04_NoMeltQDf01r1/", "05_NoMeltQDf01r2/", "06_NoMeltQDf01r3/", \
                 "10_NoMeltQCf01r1/", "11_NoMeltQCf01r2/", "12_NoMeltQCf01r3/", \
                 "16_NoMeltQPf01r1/", "17_NoMeltQPf01r2/", "18_NoMeltQPf01r3/"]
   palette = ['#C25B60','#FFBF70','#43B7C2']
   labels_a = ('diffusion (3x)','composite (3x)','composite proxy (3x)')
   figure_title = ('Quadrant annulus | $\\it{T}$$_{P0}$ = 1900 K | $\\mu$ = 0.1 | different rheology without melting')
   location='lower center'
   style = ['solid','solid','solid']
   symbol = ['^','o','s']
   legend_columns = 3
   legend_panel = 0
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



# nz values written per core
nz_core = int(128)


# customising matplotlib and figure
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.sans-serif'] = "Helvetica"
mpl.rcParams['font.size'] = 11
mpl.rcParams['axes.titlesize'] = 11
mpl.rcParams['axes.titleweight'] = "bold"
mpl.rcParams['axes.labelweight'] = "normal"
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['figure.dpi'] = 600
mpl.rcParams['figure.figsize'] = [7,6] # (3 panels, width x height, 1=600px with dpi=300)
mpl.rcParams['xtick.direction'] = "in"
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.direction'] = "in"
mpl.rcParams["pdf.use14corefonts"] = True
twin_axes=0
box_ratio=4.5
evolved_time = 1500.
time_ticks=[0.,300.,600.,900.,1200.,1500.]
minTin=1600.
maxTin=2400.
Tin_ticks=[1600,1800,2000,2200,2400]
max_mobility = 3.0
mobility_ticks=[0.0,1.0,2.0,3.0]
minTmean=2300.
maxTmean=3100.
meanT_ticks=[2300,2500,2700,2900,3100]
time_xaxis_interpolation=1500


time_Myr_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
mobility_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
steps_mobility_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
mobility_time_avg = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
time_Myr_mobility_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
moving_average_mobility_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
plateness_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
moving_average_plateness_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
plateness_time_avg = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
tmean_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))

tp_run      = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
temp200_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
SD200_run   = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
temp400_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
SD400_run   = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))

# for interpolated values
tp_run_interp = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
mobility_run_even_interp = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
tmean_run_interp = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
mobility_batch = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
time_Myr_mobility_batch = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
moving_average_mobility_batch = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
plateness_batch = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))

seconds_in_year = 31556926
seconds_in_Myear = seconds_in_year*1e6
earth_radius = 6371.e3
core_radius = 3481.e3

# define size of moving window for rolling average
window_size_mobility  = 10
window_size_plateness = 10

load_internalT = True

run = int(0)
for path in model_list:

   volume_mantle = (4 * math.pi*(earth_radius**3-core_radius**3))/3
   # quadrant analytical volume is 1/16 of full annulus, which in itself is pi/ny times actual volume
   nytot = 512
   numcores_files = 8
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
   coreT_list = []
   f = h5py.File(here+'/../'+str(model_folder)+path+'+hdf5/time_botT.h5', 'r')

   ndatasets = len(f.keys())
   for key in f.keys():
      data = f[key]
      time_list.append(float(data[0]))
      coreT_list.append(float(data[1]))
   f.close()

   time=np.array(time_list)
   coreT=np.array(coreT_list)
   time_Myr=time/seconds_in_Myear
   print('Final core T (K):        '+str('{:.0f}'.format(coreT[len(coreT)-1])))
   # print('# HDF datasets: '+str(ndatasets))

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
   # print('# frames (rprof.dat): ' + str(nframes))
   print('Evolution time (Myr):    '+str('{:.2f}'.format(sdat.tseries.time[len(sdat.tseries.time)-1]*1e-6/seconds_in_year)))
   # timesteps for initial and final states
   step1=snaps[1]
   step2=snaps[nframes-1]



   # loading Tmean from time series
   tmean_list = []

   for i in snaps:
      tmean_list.append(float(sdat.tseries['Tmean'].values[i]))
   tmean_run[run] = np.array(tmean_list)
   tmean_list.clear()
   print('Final mantle T (K):      '+str('{:.0f}'.format(tmean_run[run][len(tmean_run[run])-1])))   

   # interpolating Tmean on a much finer time axis with 1500/4500 elements
   time_interp = np.linspace(0,time_xaxis_interpolation,num=time_xaxis_interpolation+1)
   tmean_run_interp[run] = np.interp(time_interp,time_Myr_run[run],tmean_run[run])

   # Using pandas to calculate moving average (with frames) of mobility and plateness
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

   mobility_batch[run] = mobility_run[run]
   time_Myr_mobility_batch[run] = time_Myr_mobility_run[run]
   moving_average_mobility_batch[run] = moving_average_mobility_run[run]
   plateness_batch[run] = plateness_run[run]
        
   # interpolating mobility and plateness on a much finer time axis with 1500/4500 elements
   time_interp2 = np.linspace(0,time_xaxis_interpolation,num=time_xaxis_interpolation+1)
   mobility_run_even_interp[run] = np.interp(time_interp2,time_Myr_mobility_batch[run],moving_average_mobility_batch[run])


   if load_internalT:

      tp_list = []
      t200_list = []
      sd200_list = []
      t400_list = []
      sd400_list = []

      file = (here+'/../'+str(model_folder)+path+'potentialT.dat')
      f = open(file,'r')
      lines=f.readlines()
      f.close
      nlines=int(len(lines))

      for i in np.arange(0,nlines):
         p = lines[i].split()
         tp_list.append(float(p[1]))
         t200_list.append(float(p[2]))
         sd200_list.append(float(p[3]))
         t400_list.append(float(p[4]))
         sd400_list.append(float(p[5]))
      
      tp_run[run]      = np.array(tp_list)
      temp200_run[run] = np.array(t200_list)
      SD200_run[run]   = np.array(sd200_list)
      temp400_run[run] = np.array(t400_list)
      SD400_run[run]   = np.array(sd400_list)

      tp_list.clear()
      t200_list.clear()
      sd200_list.clear()
      t400_list.clear()
      sd400_list.clear()

      # interpolating Tp on a much finer time axis with 1500/4500 elements
      time_interp4 = np.linspace(0,time_xaxis_interpolation,num=time_xaxis_interpolation+1)
      tp_run_interp[run] = np.interp(time_interp4,time_Myr_run[run],(temp200_run[run]+SD200_run[run]))
      print('Final pot. T (deg C):    '+str('{:.0f}'.format(temp200_run[run][len(temp200_run[run])-1]+SD200_run[run][len(SD200_run[run])-1])))


   # plotting potential temperature, mobility, mantle temperature per model 
   fig, ax = plt.subplots(3, 1, figsize = (7,6), squeeze=False, dpi = 600)

   ax[0][0].plot(time_Myr_run[run],temp200_run[run]+SD200_run[run], linestyle='-', linewidth=1, alpha=1, color = '#C25B60')
   ax[0][0].set_xlim(0., evolved_time)
   ax[0][0].set_xscale('linear')
   ax[0][0].set_xticks(time_ticks)
   ax[0][0].set_ylim(minTin,maxTin)
   ax[0][0].set_yscale('linear')
   ax[0][0].set_yticks(Tin_ticks)   
   ax[0][0].set_ylabel('Potential T ($\\degree$C)')
   ax[0][0].set_box_aspect(1/box_ratio)

   ax[1][0].plot(time_Myr_run[run],moving_average_mobility_run[run], linestyle='-', linewidth=1, alpha=1, color = '#FFBF70')
   ax[1][0].set_xlim(0., evolved_time)
   ax[1][0].set_xscale('linear')
   ax[0][0].set_xticks(time_ticks)
   ax[1][0].set_ylim(0., max_mobility)
   ax[1][0].set_yscale('linear')
   ax[1][0].set_yticks(mobility_ticks)   
   ax[1][0].set_ylabel('Mobility',labelpad=23)
   ax[1][0].set_box_aspect(1/box_ratio)

   ax[2][0].plot(time_Myr_run[run],tmean_run[run], linestyle='-', linewidth=1, alpha=1, color = '#00647A')
   ax[2][0].set_xlim(0., evolved_time)
   ax[2][0].set_xscale('linear')
   ax[2][0].set_xticks(time_ticks)   
   ax[2][0].set_ylim(minTmean,maxTmean)
   ax[2][0].set_yscale('linear')
   ax[2][0].set_yticks(meanT_ticks)
   ax[2][0].set_ylabel('Mantle T (K)')
   ax[2][0].set_xlabel('Time (Myr)')
   ax[2][0].set_box_aspect(1/box_ratio) 

   # saving it in the model folder
   file_name='model_'+str(path[0:len(path)-1])+'_diagnostics.png'
   fig.tight_layout()
   plt.savefig(here+'/../'+str(model_folder)+path+file_name)


   run = run + 1


# calculating range and mean from all runs
if batch == 7 or batch == 8:

    # potential T
    maxTin_diff = np.maximum.reduce([tp_run_interp[0],tp_run_interp[1],tp_run_interp[2]])
    minTin_diff = np.minimum.reduce([tp_run_interp[0],tp_run_interp[1],tp_run_interp[2]])
    meanTin_diff = np.mean(np.array([maxTin_diff,minTin_diff]),axis=0)

    maxTin_disl = np.maximum.reduce([tp_run_interp[3],tp_run_interp[4],tp_run_interp[5]])
    minTin_disl = np.minimum.reduce([tp_run_interp[3],tp_run_interp[4],tp_run_interp[5]])
    meanTin_disl = np.mean(np.array([maxTin_disl,minTin_disl]),axis=0)

    maxTin_proxy = np.maximum.reduce([tp_run_interp[6],tp_run_interp[7],tp_run_interp[8]])
    minTin_proxy = np.minimum.reduce([tp_run_interp[6],tp_run_interp[7],tp_run_interp[8]])
    meanTin_proxy = np.mean(np.array([maxTin_proxy,minTin_proxy]),axis=0)


    # Mobility M
    maxM_diff = np.maximum.reduce([mobility_run_even_interp[0],mobility_run_even_interp[1],mobility_run_even_interp[2]])
    minM_diff = np.minimum.reduce([mobility_run_even_interp[0],mobility_run_even_interp[1],mobility_run_even_interp[2]])
    meanM_diff = np.mean(np.array([maxM_diff,minM_diff]),axis=0)

    maxM_disl = np.maximum.reduce([mobility_run_even_interp[3],mobility_run_even_interp[4],mobility_run_even_interp[5]])
    minM_disl = np.minimum.reduce([mobility_run_even_interp[3],mobility_run_even_interp[4],mobility_run_even_interp[5]])
    meanM_disl = np.mean(np.array([maxM_disl,minM_disl]),axis=0)

    maxM_proxy = np.maximum.reduce([mobility_run_even_interp[6],mobility_run_even_interp[7],mobility_run_even_interp[8]])
    minM_proxy = np.minimum.reduce([mobility_run_even_interp[6],mobility_run_even_interp[7],mobility_run_even_interp[8]])
    meanM_proxy = np.mean(np.array([maxM_proxy,minM_proxy]),axis=0)


    # mantle T
    maxT_diff   = np.maximum.reduce([tmean_run_interp[0],tmean_run_interp[1],tmean_run_interp[2]])
    minT_diff   = np.minimum.reduce([tmean_run_interp[0],tmean_run_interp[1],tmean_run_interp[2]])
    meanT_diff  = np.mean(np.array([maxT_diff,minT_diff]),axis=0)
    
    maxT_disl   = np.maximum.reduce([tmean_run_interp[3],tmean_run_interp[4],tmean_run_interp[5]])
    minT_disl   = np.minimum.reduce([tmean_run_interp[3],tmean_run_interp[4],tmean_run_interp[5]])
    meanT_disl  = np.mean(np.array([maxT_disl,minT_disl]),axis=0)
    
    maxT_proxy  = np.maximum.reduce([tmean_run_interp[6],tmean_run_interp[7],tmean_run_interp[8]])
    minT_proxy  = np.minimum.reduce([tmean_run_interp[6],tmean_run_interp[7],tmean_run_interp[8]])
    meanT_proxy = np.mean(np.array([maxT_proxy,minT_proxy]),axis=0)

else:
    print('No interpolation done')



# plotting
fig, ax = plt.subplots(3, 1, squeeze=False)

ax[0][0].plot(time_interp4,meanTin_diff, linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 200, mfc='#ffffff')
ax[0][0].plot(time_interp4,meanTin_disl, linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 200, mfc='#ffffff')
ax[0][0].plot(time_interp4,meanTin_proxy, linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 200, mfc='#ffffff') 

if batch == 7:
   ax[0][0].plot(time_Myr_run[0],temp200_run[0]+SD200_run[0], linestyle=style[3], linewidth=1, alpha=1, color = palette[3])
   ax[0][0].plot(time_Myr_run[4],temp200_run[4]+SD200_run[4], linestyle=style[4], linewidth=1, alpha=1, color = palette[4])
   ax[0][0].plot(time_Myr_run[6],temp200_run[6]+SD200_run[6], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])   

ax[0][0].fill_between(time_interp4,minTin_diff,maxTin_diff,linewidth=0,color = palette[0],alpha=0.2)
ax[0][0].fill_between(time_interp4,minTin_disl,maxTin_disl,linewidth=0,color = palette[1],alpha=0.2)
ax[0][0].fill_between(time_interp4,minTin_proxy,maxTin_proxy,linewidth=0,color = palette[2],alpha=0.2)

ax[0][0].set_xlim(0., evolved_time)
ax[0][0].set_xscale('linear')
ax[0][0].set_xticks(time_ticks)
ax[0][0].set_ylim(minTin,maxTin)
ax[0][0].set_yscale('linear')
ax[0][0].set_yticks(Tin_ticks)
ax[0][0].set_ylabel('Potential T ($\\degree$C)')
ax[0][0].set_box_aspect(1/box_ratio)

labels = labels_a
ax[legend_panel][0].legend(ncol = legend_columns, labels = (labels), loc=location, labelspacing = 0.5, facecolor = 'none', edgecolor = 'none', fancybox = False)



ax[1][0].plot(time_interp2,meanM_diff,  linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 200, mfc='#ffffff')
ax[1][0].plot(time_interp2,meanM_disl,  linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 200, mfc='#ffffff')
ax[1][0].plot(time_interp2,meanM_proxy, linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 200, mfc='#ffffff')

if batch == 7:
   ax[1][0].plot(time_Myr_run[0],moving_average_mobility_run[0], linestyle=style[3], linewidth=1, alpha=1, color = palette[3])
   ax[1][0].plot(time_Myr_run[4],moving_average_mobility_run[4], linestyle=style[4], linewidth=1, alpha=1, color = palette[4])
   ax[1][0].plot(time_Myr_run[6],moving_average_mobility_run[6], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])

ax[1][0].fill_between(time_interp2,minM_diff,maxM_diff,linewidth=0,color = palette[0],alpha=0.2)
ax[1][0].fill_between(time_interp2,minM_disl,maxM_disl,linewidth=0,color = palette[1],alpha=0.2)
ax[1][0].fill_between(time_interp2,minM_proxy,maxM_proxy,linewidth=0,color = palette[2],alpha=0.2)

ax[1][0].set_xlim(0., evolved_time)
ax[1][0].set_xscale('linear')
ax[1][0].set_xticks(time_ticks)  
ax[1][0].set_ylim(0., max_mobility)
ax[1][0].set_yscale('linear')
ax[1][0].set_yticks(mobility_ticks)
ax[1][0].set_ylabel('Mobility',labelpad=23)
ax[1][0].set_box_aspect(1/box_ratio)

labels = labels_a
ax[legend_panel][0].legend(ncol = legend_columns, labels = (labels), loc=location, labelspacing = 0.5, facecolor = 'none', edgecolor = 'none', fancybox = False)


ax[2][0].plot(time_interp,meanT_diff, linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 200, mfc='#ffffff')
ax[2][0].plot(time_interp,meanT_disl, linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 200, mfc='#ffffff')
ax[2][0].plot(time_interp,meanT_proxy, linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 200, mfc='#ffffff')

if batch == 7:
   ax[2][0].plot(time_Myr_run[0],tmean_run[0], linestyle=style[3], linewidth=1, alpha=1, color = palette[3])
   ax[2][0].plot(time_Myr_run[4],tmean_run[4], linestyle=style[4], linewidth=1, alpha=1, color = palette[4])
   ax[2][0].plot(time_Myr_run[6],tmean_run[6], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])

ax[2][0].fill_between(time_interp,minT_diff,maxT_diff,linewidth=0,color = palette[0],alpha=0.2)
ax[2][0].fill_between(time_interp,minT_disl,maxT_disl,linewidth=0,color = palette[1],alpha=0.2)
ax[2][0].fill_between(time_interp,minT_proxy,maxT_proxy,linewidth=0,color = palette[2],alpha=0.2)

ax[2][0].set_xlim(0., evolved_time)
ax[2][0].set_xscale('linear')
ax[2][0].set_xticks(time_ticks)   
ax[2][0].set_ylim(minTmean,maxTmean)
ax[2][0].set_yscale('linear')
ax[2][0].set_yticks(meanT_ticks)
ax[2][0].set_ylabel('Mantle T (K)')
ax[2][0].set_xlabel('Time (Myr)')
ax[2][0].set_box_aspect(1/box_ratio) 

labels = labels_a
ax[legend_panel][0].legend(ncol = legend_columns, labels = (labels), loc=location, labelspacing = 0.5, facecolor = 'none', edgecolor = 'none', fancybox = False)

ax[0][0].set_title(figure_title)

fig.text(0.01, 0.96, 'a', fontsize=15)
fig.text(0.01, 0.65, 'b', fontsize=15)
fig.text(0.01, 0.34, 'c', fontsize=15)

fig.tight_layout()
file_name='fig'+str('{:01d}'.format(batch))+'.pdf'
plt.savefig(here+'/../figures/'+file_name)


