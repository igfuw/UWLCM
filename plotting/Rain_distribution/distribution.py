import h5py
import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from collections import OrderedDict

#rain_data = ["u", "v", "w"]
#rain_data = ["u", "v", "w", "cloud_rw_mom3", "rv", "th", "RH", "aerosol_rw_mom3"]
#rain_data = ["cloud_rw_mom3"]
rain_data = ["rain_rw_mom3", "precip_rate"]
layer_thickness = 10
cloud_thresh = 1e-4 # qc at cloud base

CellVol = 50.*50*5 # hardcoded cell volume [m^3]
L_evap = 2264.76e3 # latent heat of evapporation [J/kg]

time_start = int(argv[1])
time_end = int(argv[2])
outfreq = int(argv[3])
#from_lvl = int(argv[4])
#to_lvl = int(argv[5])

directories = argv[4:len(argv):2]
labels = argv[5:len(argv):2]
print directories, labels

levels = ["ground", "cloud_base"]

for lvl in levels:
  total_arr = OrderedDict()
  for data in rain_data:
    total_arr[data] = OrderedDict()
  
  plot_labels = OrderedDict()
  tot_cloud_base_lvl = OrderedDict()
  for lab in labels:
    tot_cloud_base_lvl[lab] = np.zeros(0)
  
  # read in nx, ny, nz
  for directory, lab in zip(directories, labels):
    w3d = h5py.File(directory + "/timestep" + str(time_start).zfill(10) + ".h5", "r")["u"][:,:,:]
    nx, ny, nz = w3d.shape
    plot_labels[lab] = lab
  #  Exy_avg = OrderedDict()
  #  for data in rain_data:
  #    Exy_avg[data] = np.zeros(((nx+1)/2))
  
    
    for data in rain_data:
      total_arr[data][lab] = np.zeros(0)
  
      for t in range(time_start, time_end+1, outfreq):
        filename = directory + "/timestep" + str(t).zfill(10) + ".h5"
        print filename
  
        # find cloud base
  
        # based on cloud rw
        w3d = h5py.File(filename, "r")["cloud_rw_mom3"][:,:,:] * 4. / 3. * 3.1416 * 1e3 
        cloud_base_lvl = np.argmax(np.average(w3d, axis=(0,1)) > cloud_thresh)
        # based on RH
  #      w3d = h5py.File(filename, "r")["RH"][:,:,:] # * 4. / 3. * 3.1416 * 1e3 
  #      cloud_base_lvl = np.argmax(np.average(w3d, axis=(0,1)) > .99)
  
        tot_cloud_base_lvl[lab] = np.append(tot_cloud_base_lvl[lab], cloud_base_lvl) # done for each data, but we dont care - wont affect average
  
        print 'cloud base lvl = ', cloud_base_lvl

        if lvl == "cloud_base":
          total_arr[data][lab] = np.append(total_arr[data][lab], h5py.File(filename, "r")[data][:,:,cloud_base_lvl-layer_thickness : cloud_base_lvl])
        if lvl == "ground":
          total_arr[data][lab] = np.append(total_arr[data][lab], h5py.File(filename, "r")[data][:,:, 0 : layer_thickness ])
  
  #    hists[lab] = np.hist(total_arr, bins=100)
  #    _ = plt.hist(total_arr, bins='auto')

      # convert to typical units
      if data == "rain_rw_mom3":
        total_arr[data][lab] *= 4./3. * 3.1416 * 1e3 * 1e3 # [g/kg]
      if data == "precip_rate":
        total_arr[data][lab] *= 4./3. * 3.1416 * 1e3 / CellVol * L_evap

  
  for lab in labels:
  #  print  np.average(total_arr[lab])
    plot_labels[lab] = plot_labels[lab] + '\n <q_r> = {:.3e}'.format(np.average(total_arr["rain_rw_mom3"][lab])) \
                                        + '\n <precip flux> = {:.3e}'.format(np.average(total_arr["precip_rate"][lab])) \
                                        + '\n <cloud base lvl> = {:.2f}'.format(np.average(tot_cloud_base_lvl[lab] * 5))
    
  plt.figure(0)
  plt.rcParams.update({'font.size': 40})
  plt.figure(figsize=(40,40))
  #_ = plt.hist(total_arr["rain_rw_mom3"].values(), bins='auto', label=plot_labels.values(), density=True)
  data = total_arr["rain_rw_mom3"].values()

#  _ = plt.hist(data, bins=100, label=plot_labels.values(), density=False, histtype='step', linewidth=2)
  _ = plt.hist(data, bins=np.logspace(np.log10(1e-6), np.log10(np.amax(data)), 100), label=plot_labels.values(), density=False, histtype='step', linewidth=6)
  plt.xscale('log')
  plt.legend(loc = 'lower center')
  plt.yscale('log')
  plt.xlabel('q_r [g/kg]')
  plt.ylabel('# of cells')
  plt.savefig('rain_histo_' + lvl + '_' + str(time_start) + '_' + str(time_end) +'.png')
  
  plt.figure(1)
  plt.rcParams.update({'font.size': 40})
  plt.figure(figsize=(40,40))
  data = total_arr["precip_rate"].values()
  #_ = plt.hist(data, bins=100, label=plot_labels.values(), density=False, histtype='step', linewidth=2)
  _ = plt.hist(data, bins=np.logspace(np.log10(1e-3), np.log10(np.amax(data)), 100), label=plot_labels.values(), density=False, histtype='step', linewidth=6)
  plt.xscale('log')
  plt.legend(loc = 'lower center')
  plt.yscale('log')
  plt.xlabel('precipitation flux [W / m^2]')
  plt.ylabel('# of cells')
  plt.savefig('precip_rate_histo_' + lvl + '_' + str(time_start) + '_' + str(time_end) +'.png')
