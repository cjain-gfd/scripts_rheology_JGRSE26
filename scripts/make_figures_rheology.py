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

# Once required packages are installed, execute the script from terminal with command "python make_figures_rheology.py"

import os
import math
import numpy as np
import h5py
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from stagpy import stagyydata


### beginning of input ###
# user only has to input figure number (3,4,5,10,11) as batch
# figure 11 is supplementary figure S1
batch=3
### end of input ###


if batch == 3:
   model_folder = 'Table4_models/'
   model_list = ["01_QDf01r1/", "02_QDf01r2/", "03_QDf01r3/", \
                 "08_QCf01r1/", "09_QCf01r2/", "10_QCf01r3/", \
                 "15_QPf01r1/", "16_QPf01r2/", "17_QPf01r3/"]
   palette = ['#C25B60','#FFBF70','#43B7C2','#C25B60','#FFBF70','#43B7C2']
   labels_a = ('diffusion (3x)','composite (3x)','composite proxy (3x)','QDf01r2','QCf01r3','QPf01r1')
   figure_title = ('Quadrant annulus | $\\it{T}$$_{P0}$ = 1900 K | $\\mu$ = 0.1 | different rheology')
   location='upper center'
   style = ['solid','solid','solid','dashed','dashed','dashed']
   symbol = ['^','o','s','^','o','s']
   legend_columns = 3
   legend_panel=2
   evolved_time = 1500.
   time_ticks=[0.,300.,600.,900.,1200.,1500.]
   minTin=1450.
   maxTin=1950.
   Tin_ticks=[1450,1575,1700,1825,1950]
   max_mobility = 3.0
   mobility_ticks=[0.0,1.0,2.0,3.0]
   minTmean=2100.
   maxTmean=2600.
   meanT_ticks=[2100,2225,2350,2475,2600]
   time_xaxis_interpolation=1500
elif batch == 4:
   model_folder = 'Table4_models/'
   model_list = ["04_HDf01r1/", "05_HDf01r2/", \
                 "11_HCf01r1/", "12_HCf01r2/", \
                 "18_HPf01r1/", "19_HPf01r2/"]
   palette = ['#C25B60','#FFBF70','#43B7C2']
   labels_a = ('diffusion (2x)','composite (2x)','composite proxy (2x)')
   figure_title = ('Half annulus | $\\it{T}$$_{P0}$ = 1900 K | $\\mu$ = 0.1 | different rheology')
   location='upper center'
   style = ['solid','solid','solid']
   symbol = ['^','o','s']
   legend_columns = 3
   legend_panel=2
   evolved_time = 1500.
   time_ticks=[0.,300.,600.,900.,1200.,1500.]
   minTin=1450.
   maxTin=1950.
   Tin_ticks=[1450,1575,1700,1825,1950]
   max_mobility = 3.0
   mobility_ticks=[0.0,1.0,2.0,3.0]
   minTmean=2100.
   maxTmean=2600.
   meanT_ticks=[2100,2225,2350,2475,2600]
   time_xaxis_interpolation=1500
elif batch == 5:
   model_folder = 'Table4_models/'
   model_list = ["06_ADf01r1/", "07_ADf01r2/", \
                 "13_ACf01r1/", "14_ACf01r2/", \
                 "20_APf01r1/", "21_APf01r2/"]
   palette = ['#C25B60','#FFBF70','#43B7C2']
   labels_a = ('diffusion (2x)','composite (2x)','composite proxy (2x)')
   figure_title = ('Full annulus | $\\it{T}$$_{P0}$ = 1900 K | $\\mu$ = 0.1 | different rheology')
   location='upper center'
   style = ['solid','solid','solid']
   symbol = ['^','o','s']
   legend_columns = 3
   legend_panel=2
   evolved_time = 1500.
   time_ticks=[0.,300.,600.,900.,1200.,1500.]
   minTin=1450.
   maxTin=1950.
   Tin_ticks=[1450,1575,1700,1825,1950]
   max_mobility = 3.0
   mobility_ticks=[0.0,1.0,2.0,3.0]
   minTmean=2100.
   maxTmean=2600.
   meanT_ticks=[2100,2225,2350,2475,2600]
   time_xaxis_interpolation=1500
elif batch == 10:
   model_folder = 'Table7_models/'
   model_list = ["10_QPf005r1/", "11_QPf005r2/", "12_QPf005r3/", "13_QPf01r4/", "14_QPf01r5/", "15_QPf01r6/", "16_QPf02r1/", "17_QPf02r2/", "18_QPf02r3/"]
   palette = ['#C25B60','#FFBF70','#43B7C2','#C25B60','#FFBF70','#43B7C2','#00647A','#A800FF','#000000']
   labels_a = ('$\\mu$ = 0.05 (3x)','$\\mu$ = 0.10 (3x)','$\\mu$ = 0.20 (3x)',\
            'QPf005r2','QPf01r5','QPf02r2')
   figure_title = ('Quadrant annulus | $\\it{T}$$_{P0}$ = 1900 K | composite proxy rheology | different $\\mu$')
   location='lower center'
   style = ['solid','solid','solid','solid','solid','solid','dashed','dashed','dashed']
   symbol = ['^','^','^','o','o','o']
   legend_columns = 3
   legend_panel=2
   evolved_time = 4500.
   time_ticks=[0,500,1000,1500,2000,2500,3000,3500,4000,4500]
   minTin=1150.
   maxTin=1950.
   Tin_ticks=[1150,1350,1550,1750,1950]
   max_mobility = 3.0
   mobility_ticks=[0.0,1.0,2.0,3.0]
   minTmean=1800.
   maxTmean=2600.
   meanT_ticks=[1800,2000,2200,2400,2600]
   time_xaxis_interpolation=4500
elif batch == 11:
   model_folder = 'Table7_models/'
   model_list = ["01_QDf005r1/", "02_QDf005r2/", "03_QDf005r3/", "04_QDf01r4/", "05_QDf01r5/", "06_QDf01r6/", "07_QDf02r1/", "08_QDf02r2/", "09_QDf02r3/"]
   palette = ['#C25B60','#FFBF70','#43B7C2','#C25B60','#FFBF70','#43B7C2','#00647A','#A800FF','#000000']
   labels_a = ('$\\mu$ = 0.05 (3x)','$\\mu$ = 0.10 (3x)','$\\mu$ = 0.20 (3x)')
   figure_title = ('Quadrant annulus | $\\it{T}$$_{P0}$ = 1900 K | diffusion creep rheology | different $\\mu$')
   location='lower center'
   style = ['solid','solid','solid','dashed','dashed','dashed','dashed','dashed','dashed']
   symbol = ['^','^','^','o','o','o']
   legend_columns = 3
   legend_panel=2
   evolved_time = 4500.
   time_ticks=[0,500,1000,1500,2000,2500,3000,3500,4000,4500]
   minTin=1150.
   maxTin=1950.
   Tin_ticks=[1150,1350,1550,1750,1950]
   max_mobility = 3.0
   mobility_ticks=[0.0,1.0,2.0,3.0]
   minTmean=1800.
   maxTmean=2600.
   meanT_ticks=[1800,2000,2200,2400,2600]
   time_xaxis_interpolation=4500   

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
mpl.rcParams['font.size'] = 11
mpl.rcParams['axes.titlesize'] = 11
mpl.rcParams['axes.titleweight'] = "bold"
mpl.rcParams['axes.labelweight'] = "normal"
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['figure.dpi'] = 600
mpl.rcParams['figure.figsize'] = [7,8] # (4 panels, width x height, 1=600px with dpi=300)
mpl.rcParams['xtick.direction'] = "in"
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.direction'] = "in"
mpl.rcParams["pdf.use14corefonts"] = True
twin_axes=0
box_ratio=4.5

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

time_Myr_step_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
TotalTTGmass_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))

# for interpolated values
tp_run_interp = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
mobility_run_even_interp = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
TotalTTGmass_run_interp = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
tmean_run_interp = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
mobility_batch = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
time_Myr_mobility_batch = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
moving_average_mobility_batch = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
plateness_batch = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))

seconds_in_year = 31556926
seconds_in_Myear = seconds_in_year*1e6
earth_radius = 6371.e3
core_radius = 3481.e3
density_TTG = 2700.

# empty lists
timesteps = []
time_TTGmass   = []
TTG1   = []
TTG2   = []
TTG3   = []
frames = []
steps  = []
time = []

# define size of moving window for rolling average
window_size_mobility  = 10
window_size_plateness = 10

load_internalT = True
compute_TTGmass = True

run = int(0)
for path in model_list:

   volume_mantle = (4 * math.pi*(earth_radius**3-core_radius**3))/3
   # quadrant analytical volume is 1/16 of full annulus, which in itself is pi/ny times actual volume
   if batch == 4:
      nytot = 1024
      numcores_files = 16
      analytical_volume_half = (4 * math.pi*(earth_radius**3-core_radius**3))/3 * (math.pi)/(nytot*4)
      factor=volume_mantle/analytical_volume_half
      if run == 0:
         print('Scaling factor for half annulus: '+str(factor))      
   elif batch == 5:
      nytot = 2048
      numcores_files = 32
      analytical_volume_full = (4 * math.pi*(earth_radius**3-core_radius**3))/3 * (math.pi)/(nytot)
      factor=volume_mantle/analytical_volume_full
      if run == 0:
         print('Scaling factor for full annulus: '+str(factor))      
   else:
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
   # f = h5py.File(here+'/fig'+str(batch)+'/'+path+'time_botT.h5', 'r')
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
   # models presented in figures 3 and 4 have extra columns in their files.
   file = (here+'/../'+str(model_folder)+path+'+op/earthDIM_plates_analyse.dat')
   if batch == 4 or batch == 5:
      df = pd.read_table(file,sep='  ',engine='python', header=None, skiprows=1,usecols=range(8))
      df = df.set_axis(['istep','time', 'Rtp_vol', 'Rtp_surf', 'mobility', 'f80', 'Vrms', 'Vsurf_rms'], axis=1)
   else:
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


   # computing total mass of TTG produced
   if compute_TTGmass:
      file = (here+'/../'+str(model_folder)+path+'TTGmass.dat')
      f = open(file,'r')
      lines = f.readlines()
      f.close
      for line in lines:
         p = line.split()
         timesteps.append(int(p[0]))
         time_TTGmass.append(float(p[1]))
         TTG1.append(float(p[2]))
         TTG2.append(float(p[3]))
         TTG3.append(float(p[4]))

      finaltime=time_TTGmass[len(time_TTGmass)-1]*1e-6/seconds_in_year

      time_Myr_step = np.array(time_TTGmass)*1e-6/seconds_in_year

      # importing arrays
      TTG1v = np.array(TTG1)
      TTG2v = np.array(TTG2)
      TTG3v = np.array(TTG3)

      # declaring empty array
      TTGmass = np.array([0])

      # first sum med-P and high-P TTG
      # then add that to low-P TTG
      TTGmass = np.add(TTG1v,np.add(TTG2v,TTG3v))

      # getting cumulative mass
      TTGmass_cumulative = np.cumsum(TTGmass)
      lowP_cumulative = np.cumsum(TTG1v)
      medP_cumulative = np.cumsum(TTG2v)
      highP_cumulative = np.cumsum(TTG3v)

      print('TTG mass (scaled, kg):   '+str('{:.2e}'.format(TTGmass_cumulative[len(time_TTGmass)-1]*factor)))

      # getting cumulative volume
      TTGvolume = TTGmass_cumulative * factor * 1e-9 / density_TTG

      timesteps.clear()
      time_TTGmass.clear()
      TTG1.clear()
      TTG2.clear()
      TTG3.clear()

      time_Myr_step_run[run] = time_Myr_step
      TotalTTGmass_run[run] = TTGmass_cumulative * factor

      # interpolating TTG mass on a much finer time axis with 1500/4500 elements
      time_interp3 = np.linspace(0,time_xaxis_interpolation,num=(time_xaxis_interpolation+1)*10)
      TotalTTGmass_run_interp[run] = np.interp(time_interp3,time_Myr_step_run[run],TotalTTGmass_run[run])


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




   # plotting potential temperature, mobility, TTG mass, mantle temperature per model 
   fig, ax = plt.subplots(4, 1, figsize = (7,8), squeeze=False, dpi = 600)

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

   ax[2][0].plot(time_Myr_step_run[run],TotalTTGmass_run[run], linestyle='-', linewidth=1, alpha=1, color = '#43B7C2')
   ax[2][0].set_xlim(0., evolved_time)
   ax[2][0].set_xscale('linear')
   ax[2][0].set_xticks(time_ticks)
   ax[2][0].set_ylim(0., 8.e22)
   ax[2][0].set_yscale('linear')
   ax[2][0].set_ylabel('TTG mass (kg)',labelpad=23)   
   ax[2][0].set_box_aspect(1/box_ratio)

   ax[3][0].plot(time_Myr_run[run],tmean_run[run], linestyle='-', linewidth=1, alpha=1, color = '#00647A')
   ax[3][0].set_xlim(0., evolved_time)
   ax[3][0].set_xscale('linear')
   ax[3][0].set_xticks(time_ticks)   
   ax[3][0].set_ylim(minTmean,maxTmean)
   ax[3][0].set_yscale('linear')
   ax[3][0].set_yticks(meanT_ticks)
   ax[3][0].set_ylabel('Mantle T (K)')
   ax[3][0].set_xlabel('Time (Myr)')
   ax[3][0].set_box_aspect(1/box_ratio) 

   # saving it in the model folder
   file_name='model_'+str(path[0:len(path)-1])+'_diagnostics.png'
   fig.tight_layout()
   plt.savefig(here+'/../'+str(model_folder)+path+file_name)


   run = run + 1


# calculating range and mean from all runs
if batch == 3:

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


    # TTG mass
    maxTTG_diff = np.maximum.reduce([TotalTTGmass_run_interp[0],TotalTTGmass_run_interp[1],TotalTTGmass_run_interp[2]])
    minTTG_diff = np.minimum.reduce([TotalTTGmass_run_interp[0],TotalTTGmass_run_interp[1],TotalTTGmass_run_interp[2]])
    meanTTG_diff = np.mean(np.array([maxTTG_diff,minTTG_diff]),axis=0)

    maxTTG_disl = np.maximum.reduce([TotalTTGmass_run_interp[3],TotalTTGmass_run_interp[4],TotalTTGmass_run_interp[5]])
    minTTG_disl = np.minimum.reduce([TotalTTGmass_run_interp[3],TotalTTGmass_run_interp[4],TotalTTGmass_run_interp[5]])
    meanTTG_disl = np.mean(np.array([maxTTG_disl,minTTG_disl]),axis=0)

    maxTTG_proxy = np.maximum.reduce([TotalTTGmass_run_interp[6],TotalTTGmass_run_interp[7],TotalTTGmass_run_interp[8]])
    minTTG_proxy = np.minimum.reduce([TotalTTGmass_run_interp[6],TotalTTGmass_run_interp[7],TotalTTGmass_run_interp[8]])
    meanTTG_proxy = np.mean(np.array([maxTTG_proxy,minTTG_proxy]),axis=0)
    
    
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


elif batch == 4 or batch == 5:
    # potential T
    maxTin_diff = np.maximum.reduce([tp_run_interp[0],tp_run_interp[1]])
    minTin_diff = np.minimum.reduce([tp_run_interp[0],tp_run_interp[1]])
    meanTin_diff = np.mean(np.array([maxTin_diff,minTin_diff]),axis=0)

    maxTin_disl = np.maximum.reduce([tp_run_interp[2],tp_run_interp[3]])
    minTin_disl = np.minimum.reduce([tp_run_interp[2],tp_run_interp[3]])
    meanTin_disl = np.mean(np.array([maxTin_disl,minTin_disl]),axis=0)

    maxTin_proxy = np.maximum.reduce([tp_run_interp[4],tp_run_interp[5]])
    minTin_proxy = np.minimum.reduce([tp_run_interp[4],tp_run_interp[5]])
    meanTin_proxy = np.mean(np.array([maxTin_proxy,minTin_proxy]),axis=0)


    # Mobility M
    maxM_diff = np.maximum.reduce([mobility_run_even_interp[0],mobility_run_even_interp[1]])
    minM_diff = np.minimum.reduce([mobility_run_even_interp[0],mobility_run_even_interp[1]])
    meanM_diff = np.mean(np.array([maxM_diff,minM_diff]),axis=0)

    maxM_disl = np.maximum.reduce([mobility_run_even_interp[2],mobility_run_even_interp[3]])
    minM_disl = np.minimum.reduce([mobility_run_even_interp[2],mobility_run_even_interp[3]])
    meanM_disl = np.mean(np.array([maxM_disl,minM_disl]),axis=0)

    maxM_proxy = np.maximum.reduce([mobility_run_even_interp[4],mobility_run_even_interp[5]])
    minM_proxy = np.minimum.reduce([mobility_run_even_interp[4],mobility_run_even_interp[5]])
    meanM_proxy = np.mean(np.array([maxM_proxy,minM_proxy]),axis=0)


    # TTG mass
    maxTTG_diff = np.maximum.reduce([TotalTTGmass_run_interp[0],TotalTTGmass_run_interp[1]])
    minTTG_diff = np.minimum.reduce([TotalTTGmass_run_interp[0],TotalTTGmass_run_interp[1]])
    meanTTG_diff = np.mean(np.array([maxTTG_diff,minTTG_diff]),axis=0)

    maxTTG_disl = np.maximum.reduce([TotalTTGmass_run_interp[2],TotalTTGmass_run_interp[3]])
    minTTG_disl = np.minimum.reduce([TotalTTGmass_run_interp[2],TotalTTGmass_run_interp[3]])
    meanTTG_disl = np.mean(np.array([maxTTG_disl,minTTG_disl]),axis=0)

    maxTTG_proxy = np.maximum.reduce([TotalTTGmass_run_interp[4],TotalTTGmass_run_interp[5]])
    minTTG_proxy = np.minimum.reduce([TotalTTGmass_run_interp[4],TotalTTGmass_run_interp[5]])
    meanTTG_proxy = np.mean(np.array([maxTTG_proxy,minTTG_proxy]),axis=0)


    # mantle T
    maxT_diff   = np.maximum.reduce([tmean_run_interp[0],tmean_run_interp[1]])
    minT_diff   = np.minimum.reduce([tmean_run_interp[0],tmean_run_interp[1]])
    meanT_diff  = np.mean(np.array([maxT_diff,minT_diff]),axis=0)
    
    maxT_disl   = np.maximum.reduce([tmean_run_interp[2],tmean_run_interp[3]])
    minT_disl   = np.minimum.reduce([tmean_run_interp[2],tmean_run_interp[3]])
    meanT_disl  = np.mean(np.array([maxT_disl,minT_disl]),axis=0)
    
    maxT_proxy  = np.maximum.reduce([tmean_run_interp[4],tmean_run_interp[5]])
    minT_proxy  = np.minimum.reduce([tmean_run_interp[4],tmean_run_interp[5]])
    meanT_proxy = np.mean(np.array([maxT_proxy,minT_proxy]),axis=0)    


elif batch == 10:
   # potential T
   maxTin_disl_f005 = np.maximum.reduce([tp_run_interp[0],tp_run_interp[1],tp_run_interp[2]])
   minTin_disl_f005 = np.minimum.reduce([tp_run_interp[0],tp_run_interp[1],tp_run_interp[2]])
   meanTin_disl_f005 = np.mean(np.array([maxTin_disl_f005,minTin_disl_f005]),axis=0)

   maxTin_disl_f01 = np.maximum.reduce([tp_run_interp[3],tp_run_interp[4],tp_run_interp[5]])
   minTin_disl_f01 = np.minimum.reduce([tp_run_interp[3],tp_run_interp[4],tp_run_interp[5]])
   meanTin_disl_f01 = np.mean(np.array([maxTin_disl_f01,minTin_disl_f01]),axis=0)

   maxTin_disl_f02 = np.maximum.reduce([tp_run_interp[6],tp_run_interp[7],tp_run_interp[8]])
   minTin_disl_f02 = np.minimum.reduce([tp_run_interp[6],tp_run_interp[7],tp_run_interp[8]])
   meanTin_disl_f02 = np.mean(np.array([maxTin_disl_f02,minTin_disl_f02]),axis=0)

   # Mobility M
   maxM_disl_f005 = np.maximum.reduce([mobility_run_even_interp[0],mobility_run_even_interp[1],mobility_run_even_interp[2]])
   minM_disl_f005 = np.minimum.reduce([mobility_run_even_interp[0],mobility_run_even_interp[1],mobility_run_even_interp[2]])
   meanM_disl_f005 = np.mean(np.array([maxM_disl_f005,minM_disl_f005]),axis=0)

   maxM_disl_f01 = np.maximum.reduce([mobility_run_even_interp[3],mobility_run_even_interp[4],mobility_run_even_interp[5]])
   minM_disl_f01 = np.minimum.reduce([mobility_run_even_interp[3],mobility_run_even_interp[4],mobility_run_even_interp[5]])
   meanM_disl_f01 = np.mean(np.array([maxM_disl_f01,minM_disl_f01]),axis=0)

   maxM_disl_f02 = np.maximum.reduce([mobility_run_even_interp[6],mobility_run_even_interp[7],mobility_run_even_interp[8]])
   minM_disl_f02 = np.minimum.reduce([mobility_run_even_interp[6],mobility_run_even_interp[7],mobility_run_even_interp[8]])
   meanM_disl_f02 = np.mean(np.array([maxM_disl_f02,minM_disl_f02]),axis=0)

   # TTG mass
   maxTTG_disl_f005 = np.maximum.reduce([TotalTTGmass_run_interp[0],TotalTTGmass_run_interp[1],TotalTTGmass_run_interp[2]])
   minTTG_disl_f005 = np.minimum.reduce([TotalTTGmass_run_interp[0],TotalTTGmass_run_interp[1],TotalTTGmass_run_interp[2]])
   meanTTG_disl_f005 = np.mean(np.array([maxTTG_disl_f005,minTTG_disl_f005]),axis=0)

   maxTTG_disl_f01 = np.maximum.reduce([TotalTTGmass_run_interp[3],TotalTTGmass_run_interp[4],TotalTTGmass_run_interp[5]])
   minTTG_disl_f01 = np.minimum.reduce([TotalTTGmass_run_interp[3],TotalTTGmass_run_interp[4],TotalTTGmass_run_interp[5]])
   meanTTG_disl_f01 = np.mean(np.array([maxTTG_disl_f01,minTTG_disl_f01]),axis=0)

   maxTTG_disl_f02 = np.maximum.reduce([TotalTTGmass_run_interp[6],TotalTTGmass_run_interp[7],TotalTTGmass_run_interp[8]])
   minTTG_disl_f02 = np.minimum.reduce([TotalTTGmass_run_interp[6],TotalTTGmass_run_interp[7],TotalTTGmass_run_interp[8]])
   meanTTG_disl_f02 = np.mean(np.array([maxTTG_disl_f02,minTTG_disl_f02]),axis=0)

   # mantle T
   maxT_disl_f005 = np.maximum.reduce([tmean_run_interp[0],tmean_run_interp[1],tmean_run_interp[2]])
   minT_disl_f005 = np.minimum.reduce([tmean_run_interp[0],tmean_run_interp[1],tmean_run_interp[2]])
   meanT_disl_f005 = np.mean(np.array([maxT_disl_f005,minT_disl_f005]),axis=0)

   maxT_disl_f01 = np.maximum.reduce([tmean_run_interp[3],tmean_run_interp[4],tmean_run_interp[5]])
   minT_disl_f01 = np.minimum.reduce([tmean_run_interp[3],tmean_run_interp[4],tmean_run_interp[5]])
   meanT_disl_f01 = np.mean(np.array([maxT_disl_f01,minT_disl_f01]),axis=0)

   maxT_disl_f02 = np.maximum.reduce([tmean_run_interp[6],tmean_run_interp[7],tmean_run_interp[8]])
   minT_disl_f02 = np.minimum.reduce([tmean_run_interp[6],tmean_run_interp[7],tmean_run_interp[8]])
   meanT_disl_f02 = np.mean(np.array([maxT_disl_f02,minT_disl_f02]),axis=0)

elif batch == 11:
   # potential T
   maxTin_diff_f005 = np.maximum.reduce([tp_run_interp[0],tp_run_interp[1],tp_run_interp[2]])
   minTin_diff_f005 = np.minimum.reduce([tp_run_interp[0],tp_run_interp[1],tp_run_interp[2]])
   meanTin_diff_f005 = np.mean(np.array([maxTin_diff_f005,minTin_diff_f005]),axis=0)

   maxTin_diff_f01 = np.maximum.reduce([tp_run_interp[3],tp_run_interp[4],tp_run_interp[5]])
   minTin_diff_f01 = np.minimum.reduce([tp_run_interp[3],tp_run_interp[4],tp_run_interp[5]])
   meanTin_diff_f01 = np.mean(np.array([maxTin_diff_f01,minTin_diff_f01]),axis=0)

   maxTin_diff_f02 = np.maximum.reduce([tp_run_interp[6],tp_run_interp[7],tp_run_interp[8]])
   minTin_diff_f02 = np.minimum.reduce([tp_run_interp[6],tp_run_interp[7],tp_run_interp[8]])
   meanTin_diff_f02 = np.mean(np.array([maxTin_diff_f02,minTin_diff_f02]),axis=0)

   # Mobility M
   maxM_diff_f005 = np.maximum.reduce([mobility_run_even_interp[0],mobility_run_even_interp[1],mobility_run_even_interp[2]])
   minM_diff_f005 = np.minimum.reduce([mobility_run_even_interp[0],mobility_run_even_interp[1],mobility_run_even_interp[2]])
   meanM_diff_f005 = np.mean(np.array([maxM_diff_f005,minM_diff_f005]),axis=0)

   maxM_diff_f01 = np.maximum.reduce([mobility_run_even_interp[3],mobility_run_even_interp[4],mobility_run_even_interp[5]])
   minM_diff_f01 = np.minimum.reduce([mobility_run_even_interp[3],mobility_run_even_interp[4],mobility_run_even_interp[5]])
   meanM_diff_f01 = np.mean(np.array([maxM_diff_f01,minM_diff_f01]),axis=0)

   maxM_diff_f02 = np.maximum.reduce([mobility_run_even_interp[6],mobility_run_even_interp[7],mobility_run_even_interp[8]])
   minM_diff_f02 = np.minimum.reduce([mobility_run_even_interp[6],mobility_run_even_interp[7],mobility_run_even_interp[8]])
   meanM_diff_f02 = np.mean(np.array([maxM_diff_f02,minM_diff_f02]),axis=0)

   # TTG mass
   maxTTG_diff_f005 = np.maximum.reduce([TotalTTGmass_run_interp[0],TotalTTGmass_run_interp[1],TotalTTGmass_run_interp[2]])
   minTTG_diff_f005 = np.minimum.reduce([TotalTTGmass_run_interp[0],TotalTTGmass_run_interp[1],TotalTTGmass_run_interp[2]])
   meanTTG_diff_f005 = np.mean(np.array([maxTTG_diff_f005,minTTG_diff_f005]),axis=0)

   maxTTG_diff_f01 = np.maximum.reduce([TotalTTGmass_run_interp[3],TotalTTGmass_run_interp[4],TotalTTGmass_run_interp[5]])
   minTTG_diff_f01 = np.minimum.reduce([TotalTTGmass_run_interp[3],TotalTTGmass_run_interp[4],TotalTTGmass_run_interp[5]])
   meanTTG_diff_f01 = np.mean(np.array([maxTTG_diff_f01,minTTG_diff_f01]),axis=0)

   maxTTG_diff_f02 = np.maximum.reduce([TotalTTGmass_run_interp[6],TotalTTGmass_run_interp[7],TotalTTGmass_run_interp[8]])
   minTTG_diff_f02 = np.minimum.reduce([TotalTTGmass_run_interp[6],TotalTTGmass_run_interp[7],TotalTTGmass_run_interp[8]])
   meanTTG_diff_f02 = np.mean(np.array([maxTTG_diff_f02,minTTG_diff_f02]),axis=0)

   # mantle T
   maxT_diff_f005 = np.maximum.reduce([tmean_run_interp[0],tmean_run_interp[1],tmean_run_interp[2]])
   minT_diff_f005 = np.minimum.reduce([tmean_run_interp[0],tmean_run_interp[1],tmean_run_interp[2]])
   meanT_diff_f005 = np.mean(np.array([maxT_diff_f005,minT_diff_f005]),axis=0)

   maxT_diff_f01 = np.maximum.reduce([tmean_run_interp[3],tmean_run_interp[4],tmean_run_interp[5]])
   minT_diff_f01 = np.minimum.reduce([tmean_run_interp[3],tmean_run_interp[4],tmean_run_interp[5]])
   meanT_diff_f01 = np.mean(np.array([maxT_diff_f01,minT_diff_f01]),axis=0)

   maxT_diff_f02 = np.maximum.reduce([tmean_run_interp[6],tmean_run_interp[7],tmean_run_interp[8]])
   minT_diff_f02 = np.minimum.reduce([tmean_run_interp[6],tmean_run_interp[7],tmean_run_interp[8]])
   meanT_diff_f02 = np.mean(np.array([maxT_diff_f02,minT_diff_f02]),axis=0)

else:
    print('No interpolation done')


# plotting
fig, ax = plt.subplots(4, 1, squeeze=False)

if batch == 3:
   ax[0][0].plot(time_interp4,meanTin_diff, linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 200, mfc='#ffffff')
   ax[0][0].plot(time_interp4,meanTin_disl, linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 200, mfc='#ffffff')
   ax[0][0].plot(time_interp4,meanTin_proxy, linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 200, mfc='#ffffff') 
   ax[0][0].plot(time_Myr_run[1],temp200_run[1]+SD200_run[1], linestyle=style[3], linewidth=1, alpha=1, color = palette[3])
   ax[0][0].plot(time_Myr_run[5],temp200_run[5]+SD200_run[5], linestyle=style[4], linewidth=1, alpha=1, color = palette[4])
   ax[0][0].plot(time_Myr_run[6],temp200_run[6]+SD200_run[6], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])
   ax[0][0].fill_between(time_interp4,minTin_diff,maxTin_diff,linewidth=0,color = palette[0],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTin_disl,maxTin_disl,linewidth=0,color = palette[1],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTin_proxy,maxTin_proxy,linewidth=0,color = palette[2],alpha=0.2)   
elif batch == 4 or batch == 5:
   ax[0][0].plot(time_interp4,meanTin_diff, linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 200, mfc='#ffffff')
   ax[0][0].plot(time_interp4,meanTin_disl, linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 200, mfc='#ffffff')
   ax[0][0].plot(time_interp4,meanTin_proxy, linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 200, mfc='#ffffff')
   ax[0][0].fill_between(time_interp4,minTin_diff,maxTin_diff,linewidth=0,color = palette[0],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTin_disl,maxTin_disl,linewidth=0,color = palette[1],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTin_proxy,maxTin_proxy,linewidth=0,color = palette[2],alpha=0.2)
elif batch == 10:
   ax[0][0].plot(time_interp4,meanTin_disl_f005, linestyle=style[3], linewidth=1, alpha=1, color = palette[3],marker=symbol[3], markersize=4, mec=palette[3], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[0][0].plot(time_interp4,meanTin_disl_f01,  linestyle=style[4], linewidth=1, alpha=1, color = palette[4],marker=symbol[4], markersize=4, mec=palette[4], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[0][0].plot(time_interp4,meanTin_disl_f02,  linestyle=style[5], linewidth=1, alpha=1, color = palette[5],marker=symbol[5], markersize=4, mec=palette[5], mew='0.5', markevery = 500, mfc='#ffffff')   
   ax[0][0].plot(time_Myr_run[1],temp200_run[1]+SD200_run[1], linestyle=style[6], linewidth=1, alpha=1, color = palette[6])
   ax[0][0].plot(time_Myr_run[4],temp200_run[4]+SD200_run[4], linestyle=style[7], linewidth=1, alpha=1, color = palette[7])
   ax[0][0].plot(time_Myr_run[7],temp200_run[7]+SD200_run[7], linestyle=style[8], linewidth=1, alpha=1, color = palette[8])    
   ax[0][0].fill_between(time_interp4,minTin_disl_f005,maxTin_disl_f005,linewidth=0,color = palette[3],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTin_disl_f01 ,maxTin_disl_f01 ,linewidth=0,color = palette[4],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTin_disl_f02 ,maxTin_disl_f02 ,linewidth=0,color = palette[5],alpha=0.2)
elif batch == 11:
   ax[0][0].plot(time_interp4,meanTin_diff_f005, linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 400, mfc='#ffffff')
   ax[0][0].plot(time_interp4,meanTin_diff_f01,  linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 400, mfc='#ffffff')
   ax[0][0].plot(time_interp4,meanTin_diff_f02,  linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 400, mfc='#ffffff')
   ax[0][0].fill_between(time_interp4,minTin_diff_f005,maxTin_diff_f005,linewidth=0,color = palette[0],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTin_diff_f01 ,maxTin_diff_f01 ,linewidth=0,color = palette[1],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTin_diff_f02 ,maxTin_diff_f02 ,linewidth=0,color = palette[2],alpha=0.2)

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


if batch == 3:
   ax[1][0].plot(time_interp2,meanM_diff,  linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 200, mfc='#ffffff')
   ax[1][0].plot(time_interp2,meanM_disl,  linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 200, mfc='#ffffff')
   ax[1][0].plot(time_interp2,meanM_proxy, linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 200, mfc='#ffffff')
   ax[1][0].plot(time_Myr_run[1],moving_average_mobility_run[1], linestyle=style[3], linewidth=1, alpha=1, color = palette[3])
   ax[1][0].plot(time_Myr_run[5],moving_average_mobility_run[5], linestyle=style[4], linewidth=1, alpha=1, color = palette[4])
   ax[1][0].plot(time_Myr_run[6],moving_average_mobility_run[6], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])
   ax[1][0].fill_between(time_interp2,minM_diff,maxM_diff,linewidth=0,color = palette[0],alpha=0.2)
   ax[1][0].fill_between(time_interp2,minM_disl,maxM_disl,linewidth=0,color = palette[1],alpha=0.2)
   ax[1][0].fill_between(time_interp2,minM_proxy,maxM_proxy,linewidth=0,color = palette[2],alpha=0.2)   
elif batch == 4 or batch == 5:
   ax[1][0].plot(time_interp2,meanM_diff,  linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 200, mfc='#ffffff')
   ax[1][0].plot(time_interp2,meanM_disl,  linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 200, mfc='#ffffff')
   ax[1][0].plot(time_interp2,meanM_proxy, linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 200, mfc='#ffffff')
   ax[1][0].fill_between(time_interp2,minM_diff,maxM_diff,linewidth=0,color = palette[0],alpha=0.2)
   ax[1][0].fill_between(time_interp2,minM_disl,maxM_disl,linewidth=0,color = palette[1],alpha=0.2)
   ax[1][0].fill_between(time_interp2,minM_proxy,maxM_proxy,linewidth=0,color = palette[2],alpha=0.2)   
elif batch == 10:
   ax[1][0].plot(time_interp2,meanM_disl_f005, linestyle=style[3], linewidth=1, alpha=1, color = palette[3],marker=symbol[3], markersize=4, mec=palette[3], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[1][0].plot(time_interp2,meanM_disl_f01,  linestyle=style[4], linewidth=1, alpha=1, color = palette[4],marker=symbol[4], markersize=4, mec=palette[4], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[1][0].plot(time_interp2,meanM_disl_f02,  linestyle=style[5], linewidth=1, alpha=1, color = palette[5],marker=symbol[5], markersize=4, mec=palette[5], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[1][0].plot(time_Myr_run[1],moving_average_mobility_run[1], linestyle=style[6], linewidth=1, alpha=1, color = palette[6])
   ax[1][0].plot(time_Myr_run[4],moving_average_mobility_run[4], linestyle=style[7], linewidth=1, alpha=1, color = palette[7])
   ax[1][0].plot(time_Myr_run[7],moving_average_mobility_run[7], linestyle=style[8], linewidth=1, alpha=1, color = palette[8])   
   ax[1][0].fill_between(time_interp2,minM_disl_f005,maxM_disl_f005,linewidth=0,color = palette[3],alpha=0.2)
   ax[1][0].fill_between(time_interp2,minM_disl_f01 ,maxM_disl_f01 ,linewidth=0,color = palette[4],alpha=0.2)
   ax[1][0].fill_between(time_interp2,minM_disl_f02 ,maxM_disl_f02 ,linewidth=0,color = palette[5],alpha=0.2)
elif batch == 11:
   ax[1][0].plot(time_interp2,meanM_diff_f005, linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 400, mfc='#ffffff')
   ax[1][0].plot(time_interp2,meanM_diff_f01,  linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 400, mfc='#ffffff')
   ax[1][0].plot(time_interp2,meanM_diff_f02,  linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 400, mfc='#ffffff')
   ax[1][0].fill_between(time_interp2,minM_diff_f005,maxM_diff_f005,linewidth=0,color = palette[0],alpha=0.2)
   ax[1][0].fill_between(time_interp2,minM_diff_f01 ,maxM_diff_f01 ,linewidth=0,color = palette[1],alpha=0.2)
   ax[1][0].fill_between(time_interp2,minM_diff_f02 ,maxM_diff_f02 ,linewidth=0,color = palette[2],alpha=0.2)

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


if batch == 3:
   ax[2][0].plot(time_interp3,meanTTG_diff, linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 2000, mfc='#ffffff')
   ax[2][0].plot(time_interp3,meanTTG_disl, linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 2000, mfc='#ffffff')
   ax[2][0].plot(time_interp3,meanTTG_proxy, linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 2000, mfc='#ffffff')
   ax[2][0].plot(time_Myr_step_run[1],TotalTTGmass_run[1], linestyle=style[3], linewidth=1, alpha=1, color = palette[3])
   ax[2][0].plot(time_Myr_step_run[5],TotalTTGmass_run[5], linestyle=style[4], linewidth=1, alpha=1, color = palette[4])
   ax[2][0].plot(time_Myr_step_run[6],TotalTTGmass_run[6], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])
   ax[2][0].fill_between(time_interp3,minTTG_diff,maxTTG_diff,linewidth=0,color = palette[0],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTG_disl,maxTTG_disl,linewidth=0,color = palette[1],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTG_proxy,maxTTG_proxy,linewidth=0,color = palette[2],alpha=0.2)      
elif batch == 4 or batch == 5:
   ax[2][0].plot(time_interp3,meanTTG_diff, linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 2000, mfc='#ffffff')
   ax[2][0].plot(time_interp3,meanTTG_disl, linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 2000, mfc='#ffffff')
   ax[2][0].plot(time_interp3,meanTTG_proxy, linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 2000, mfc='#ffffff')   
   ax[2][0].fill_between(time_interp3,minTTG_diff,maxTTG_diff,linewidth=0,color = palette[0],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTG_disl,maxTTG_disl,linewidth=0,color = palette[1],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTG_proxy,maxTTG_proxy,linewidth=0,color = palette[2],alpha=0.2)      
elif batch == 10:
   ax[2][0].plot(time_interp3,meanTTG_disl_f005, linestyle=style[3], linewidth=1, alpha=1, color = palette[3],marker=symbol[3], markersize=4, mec=palette[3], mew='0.5', markevery = 5000, mfc='#ffffff')
   ax[2][0].plot(time_interp3,meanTTG_disl_f01,  linestyle=style[4], linewidth=1, alpha=1, color = palette[4],marker=symbol[4], markersize=4, mec=palette[4], mew='0.5', markevery = 5000, mfc='#ffffff')
   ax[2][0].plot(time_interp3,meanTTG_disl_f02,  linestyle=style[5], linewidth=1, alpha=1, color = palette[5],marker=symbol[5], markersize=4, mec=palette[5], mew='0.5', markevery = 5000, mfc='#ffffff')   
   ax[2][0].plot(time_Myr_step_run[1],TotalTTGmass_run[1], linestyle=style[6], linewidth=1, alpha=1, color = palette[6])
   ax[2][0].plot(time_Myr_step_run[4],TotalTTGmass_run[4], linestyle=style[7], linewidth=1, alpha=1, color = palette[7])
   ax[2][0].plot(time_Myr_step_run[7],TotalTTGmass_run[7], linestyle=style[8], linewidth=1, alpha=1, color = palette[8])   
   ax[2][0].fill_between(time_interp3,minTTG_disl_f005,maxTTG_disl_f005,linewidth=0,color = palette[3],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTG_disl_f01, maxTTG_disl_f01, linewidth=0,color = palette[4],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTG_disl_f02, maxTTG_disl_f02, linewidth=0,color = palette[5],alpha=0.2)
elif batch == 11:
   ax[2][0].plot(time_interp3,meanTTG_diff_f005, linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 4000, mfc='#ffffff')
   ax[2][0].plot(time_interp3,meanTTG_diff_f01,  linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 4000, mfc='#ffffff')
   ax[2][0].plot(time_interp3,meanTTG_diff_f02,  linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 4000, mfc='#ffffff')
   ax[2][0].fill_between(time_interp3,minTTG_diff_f005,maxTTG_diff_f005,linewidth=0,color = palette[0],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTG_diff_f01, maxTTG_diff_f01, linewidth=0,color = palette[1],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTG_diff_f02, maxTTG_diff_f02, linewidth=0,color = palette[2],alpha=0.2)    

ax[2][0].set_xlim(0., evolved_time)    
ax[2][0].set_xscale('linear')
ax[2][0].set_xticks(time_ticks)
ax[2][0].set_ylim(0., 8.e22)
ax[2][0].set_yscale('linear')
ax[2][0].set_ylabel('TTG mass (kg)',labelpad=23)
ax[2][0].set_box_aspect(1/box_ratio)    

labels = labels_a
ax[legend_panel][0].legend(ncol = legend_columns, labels = (labels), loc=location, labelspacing = 0.5, facecolor = 'none', edgecolor = 'none', fancybox = False)


if batch == 3:
   ax[3][0].plot(time_interp,meanT_diff, linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 200, mfc='#ffffff')
   ax[3][0].plot(time_interp,meanT_disl, linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 200, mfc='#ffffff')
   ax[3][0].plot(time_interp,meanT_proxy, linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 200, mfc='#ffffff')
   ax[3][0].plot(time_Myr_run[1],tmean_run[1], linestyle=style[3], linewidth=1, alpha=1, color = palette[3])
   ax[3][0].plot(time_Myr_run[5],tmean_run[5], linestyle=style[4], linewidth=1, alpha=1, color = palette[4])
   ax[3][0].plot(time_Myr_run[6],tmean_run[6], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])
   ax[3][0].fill_between(time_interp,minT_diff,maxT_diff,linewidth=0,color = palette[0],alpha=0.2)
   ax[3][0].fill_between(time_interp,minT_disl,maxT_disl,linewidth=0,color = palette[1],alpha=0.2)
   ax[3][0].fill_between(time_interp,minT_proxy,maxT_proxy,linewidth=0,color = palette[2],alpha=0.2)   
elif batch == 4 or batch == 5:
   ax[3][0].plot(time_interp,meanT_diff, linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 200, mfc='#ffffff')
   ax[3][0].plot(time_interp,meanT_disl, linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 200, mfc='#ffffff')
   ax[3][0].plot(time_interp,meanT_proxy, linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 200, mfc='#ffffff')
   ax[3][0].fill_between(time_interp,minT_diff,maxT_diff,linewidth=0,color = palette[0],alpha=0.2)
   ax[3][0].fill_between(time_interp,minT_disl,maxT_disl,linewidth=0,color = palette[1],alpha=0.2)
   ax[3][0].fill_between(time_interp,minT_proxy,maxT_proxy,linewidth=0,color = palette[2],alpha=0.2)   
elif batch == 10:
   ax[3][0].plot(time_interp,meanT_disl_f005, linestyle=style[3], linewidth=1, alpha=1, color = palette[3],marker=symbol[3], markersize=4, mec=palette[3], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[3][0].plot(time_interp,meanT_disl_f01,  linestyle=style[4], linewidth=1, alpha=1, color = palette[4],marker=symbol[4], markersize=4, mec=palette[4], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[3][0].plot(time_interp,meanT_disl_f02,  linestyle=style[5], linewidth=1, alpha=1, color = palette[5],marker=symbol[5], markersize=4, mec=palette[5], mew='0.5', markevery = 500, mfc='#ffffff')   
   ax[3][0].plot(time_Myr_run[1],tmean_run[1], linestyle=style[6], linewidth=1, alpha=1, color = palette[6])
   ax[3][0].plot(time_Myr_run[4],tmean_run[4], linestyle=style[7], linewidth=1, alpha=1, color = palette[7])
   ax[3][0].plot(time_Myr_run[7],tmean_run[7], linestyle=style[8], linewidth=1, alpha=1, color = palette[8])   
   ax[3][0].fill_between(time_interp,minT_disl_f005,maxT_disl_f005,linewidth=0,color = palette[3],alpha=0.2)
   ax[3][0].fill_between(time_interp,minT_disl_f01, maxT_disl_f01, linewidth=0,color = palette[4],alpha=0.2)
   ax[3][0].fill_between(time_interp,minT_disl_f02, maxT_disl_f02, linewidth=0,color = palette[5],alpha=0.2)
elif batch == 11:
   ax[3][0].plot(time_interp,meanT_diff_f005, linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 400, mfc='#ffffff')
   ax[3][0].plot(time_interp,meanT_diff_f01,  linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 400, mfc='#ffffff')
   ax[3][0].plot(time_interp,meanT_diff_f02,  linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 400, mfc='#ffffff')
   ax[3][0].fill_between(time_interp,minT_diff_f005,maxT_diff_f005,linewidth=0,color = palette[0],alpha=0.2)
   ax[3][0].fill_between(time_interp,minT_diff_f01, maxT_diff_f01, linewidth=0,color = palette[1],alpha=0.2)
   ax[3][0].fill_between(time_interp,minT_diff_f02, maxT_diff_f02, linewidth=0,color = palette[2],alpha=0.2)    

ax[3][0].set_xlim(0., evolved_time)
ax[3][0].set_xscale('linear')
ax[3][0].set_xticks(time_ticks)   

ax[3][0].set_ylim(minTmean,maxTmean)
ax[3][0].set_yscale('linear')
ax[3][0].set_yticks(meanT_ticks)
ax[3][0].set_ylabel('Mantle T (K)')
ax[3][0].set_xlabel('Time (Myr)')
ax[3][0].set_box_aspect(1/box_ratio) 

ax[0][0].set_title(figure_title)

fig.text(0.01, 0.97, 'a', fontsize=15)
fig.text(0.01, 0.73, 'b', fontsize=15)
fig.text(0.01, 0.50, 'c', fontsize=15)
fig.text(0.01, 0.25, 'd', fontsize=15)

fig.tight_layout()
file_name='fig'+str('{:01d}'.format(batch))+'.pdf'
plt.savefig(here+'/../figures/'+file_name,transparent=True)


