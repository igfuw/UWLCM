from scipy.io import netcdf
import numpy as np
from sys import argv
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from bisect import bisect_left
from matplotlib import rc
import os

import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../latex_labels/")
from latex_labels import *

labeldict = {
 0 : "(a)",
 1 : "(b)",
 2 : "(c)",
 3 : "(d)",
 4 : "(e)",
 5 : "(f)",
 6 : "(g)",
 7 : "(h)",
 8 : "(i)",
 9 : "(j)"
}

def read_my_array(file_obj):
  arr_name = file_obj.readline()
  file_obj.readline() # discarded line with size of the array
  line = file_obj.readline()
  line = line.split(" ")
  del line[0]
  del line[len(line)-1]
  arr = map(float,line)
  return np.array(arr), arr_name

def read_my_var(file_obj, var_name):
  while True:
    arr, name = read_my_array(file_obj)
    if(str(name).strip() == str(var_name).strip()):
      break
  return arr

def plot_my_array(axarr, plot_iter, time, val, nploty, xlabel=None, ylabel=None, varlabel=None , linestyle='--', dashes=(5,2), xlim=None):
  x = int(plot_iter / nploty)
  y = plot_iter % nploty
  if varlabel != None:
    axarr[x, y].plot(time, val, label=varlabel, linestyle=linestyle, linewidth=1, dashes=dashes)
  else:
    axarr[x, y].plot(time, val, linestyle=linestyle, linewidth=1, dashes=dashes)
  if xlabel:
    axarr[x, y].set_xlabel(xlabel)
  if ylabel:
    axarr[x, y].set_ylabel(ylabel)
  if xlim:
    axarr[x, y].set_xlim(xlim)

def plot_profiles(var_list, nplotx, nploty, axarr, show_bin=False, suffix='', reference=True, ylabel=''):
  #read my results
  profiles_files_names = []
  profiles_labels = []
  file_no = np.arange(1, len(sys.argv)-1 , 2)
  print file_no
  for no in file_no:
    profiles_files_names.append(argv[no]+suffix)
    profiles_labels.append(argv[no+1])
  
  plot_iter = 0
  for var in var_list:
    label_counter = 0
    for file_name in profiles_files_names:
    #  try:
      profiles_file = open(file_name, "r")
      my_pos = read_my_var(profiles_file, "position")
      my_res = read_my_var(profiles_file, var)
   
   #     print 'mean nc in cloud cells: ' , np.mean(my_nc[my_nc>20])
    
      profiles_file.close()
    
      linestyles = ['--', '-.', ':']
      dashList = [(3,1),(1,1),(4,1,1,1),(4,2)] 
      if(var == "base_prflux_vs_clhght"):
        plot_my_array(axarr, plot_iter, my_res, my_pos, nploty, xlabel=var_labels[var], ylabel="cloudy column height [m]", varlabel=profiles_labels[label_counter], dashes = dashList[label_counter % len(dashList)], xlim = (0,600))
      else:
        plot_my_array(axarr, plot_iter, my_res, my_pos, nploty, xlabel=var_labels[var], ylabel=ylabel, varlabel=profiles_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
   #   except:
   #     print 'error opening file: ', file_name
   #     my_pos = 0
   #     my_res = 0
      label_counter = label_counter+1
    plot_iter += 1

def plot_series(var_list, nplotx, nploty, axarr, show_bin=False, suffix='', xlabel=''):
  #read my results
  series_files_names = []
  series_labels = []
  file_no = np.arange(1, len(sys.argv)-1 , 2)
  for no in file_no:
    series_files_names.append(argv[no] + suffix)
    series_labels.append(argv[no+1])

  plot_iter = 0
  for var in var_list:
    label_counter=0
    for file_name in series_files_names:
      print file_name, var
      series_file = open(file_name, "r")
      my_times = read_my_var(series_file, "position")
      my_res = read_my_var(series_file, var)
      
      # rescale time to hours
      my_times = my_times / 3600.
      
      series_file.close()
  
      linestyles = ['--', '-.', ':']
      dashList = [(3,1),(1,1),(4,1,1,1),(4,2)] 
      plot_my_array(axarr, plot_iter, my_times, my_res, nploty, xlabel=xlabel, ylabel=var_labels[var], varlabel=series_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
      label_counter+=1
    plot_iter = plot_iter + 1
