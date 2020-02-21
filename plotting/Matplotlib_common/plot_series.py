import numpy as np
from sys import argv

from latex_labels import var_labels
from read_UWLCM_arrays import *

def plot_series(var_list, plot_iter, nplotx, nploty, axarr, xscaledict, yscaledict, xlimdict, ylimdict, show_bin=False, suffix='', xlabel='', file_names=[], file_labels=[]):
  # if file names are not defined, read them and labels from command line
  if len(file_names)==0:
    file_no = np.arange(1, len(argv)-1 , 2)
    for no in file_no:
      file_names.append(argv[no] + suffix)
      file_labels.append(argv[no+1])

  for var in var_list:
    label_counter=0
    for file_name in file_names:
      print file_name, var
      series_file = open(file_name, "r")
      my_times = read_my_var(series_file, "position")
      my_res = read_my_var(series_file, var)
      
      # rescale time to hours
      my_times = my_times / 3600.
      
      series_file.close()
  
      linestyles = ['--', '-.', ':']
      dashList = [(3,1),(1,1),(4,1,1,1),(4,2)] 
      plot_my_array(axarr, plot_iter, my_times, my_res, nploty, xlabel=xlabel, ylabel=var_labels[var], varlabel=file_labels[label_counter], dashes = dashList[label_counter % len(dashList)], xscale=xscaledict[var], yscale=yscaledict[var], xlim=xlimdict[var], ylim=ylimdict[var])
      label_counter+=1
    plot_iter = plot_iter + 1
  return plot_iter
