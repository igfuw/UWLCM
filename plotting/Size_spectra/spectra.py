import h5py
import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from collections import OrderedDict

plt.rcParams.update({'font.size': 40})
plt.figure(figsize=(40,40))

size_data = {"rw" : 29, "rd" : 21}

#hardcoded bin edges, need to match UWLCM
left_edges = {"rw": np.zeros(30), "rd": np.zeros(22)}
bin_centers = {"rw": np.zeros(29), "rd": np.zeros(21)}
bin_width = {"rw": np.zeros(29), "rd": np.zeros(21)}

for i in np.arange(30):
  left_edges["rw"][i] = 10**(-3 + i * .2) * 1e-6 # [m]
for i in np.arange(29):
  bin_centers["rw"][i] = 0.5 * (left_edges["rw"][i] + left_edges["rw"][i+1])
  bin_width["rw"][i] = (left_edges["rw"][i+1] - left_edges["rw"][i])

for i in np.arange(22):
  left_edges["rd"][i] = 10**(-3 + i * .2) * 1e-6 # [m]
for i in np.arange(21):
  bin_centers["rd"][i] = 0.5 * (left_edges["rd"][i] + left_edges["rd"][i+1])
  bin_width["rd"][i] = (left_edges["rd"][i+1] - left_edges["rd"][i])

print left_edges
print bin_centers
print bin_width

data_names = {}
for rwrd in size_data:
  bin_no = np.arange(0,size_data[rwrd])
  data_names[rwrd] = []
  for no in bin_no:
    data_names[rwrd] = np.append(data_names[rwrd], rwrd + "_rng" + str(no).zfill(3) + "_mom0")

print data_names

layer_thickness = 10
cloud_thresh = 1e-8

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
  for data in np.append(data_names["rw"], data_names["rd"]):
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
  #  for data in data_names:
  #    Exy_avg[data] = np.zeros(((nx+1)/2))
  
    for data in np.append(data_names["rw"], data_names["rd"]):
      total_arr[data][lab] = np.zeros(0)
  
      for t in range(time_start, time_end+1, outfreq):
        filename = directory + "/timestep" + str(t).zfill(10) + ".h5"
        print filename
  
        # find cloud base
  
        # based on cloud rw
        w3d = h5py.File(filename, "r")["cloud_rw_mom3"][:,:,:] # * 4. / 3. * 3.1416 * 1e3 
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
        print total_arr[data][lab].shape
  
  #    hists[lab] = np.hist(total_arr, bins=100)
  #    _ = plt.hist(total_arr, bins='auto')
  
  for lab in labels:
  #  print  np.average(total_arr[lab])
    plot_labels[lab] = plot_labels[lab] + ' <cloud base lvl> = {:.2f}'.format(np.average(tot_cloud_base_lvl[lab]))
    
  plt.rcParams.update({'font.size': 30})
  plt.figure(figsize=(40,40))
  #_ = plt.hist(total_arr["rain_rw_mom3"].values(), bins='auto', label=plot_labels.values(), density=True)

  avg_conc = {}
  avg_conc_arr = {}
  for rwrd in size_data:
    avg_conc_arr[rwrd] = {}
    for lab in labels:
      avg_conc_arr[rwrd][lab] = np.zeros(size_data[rwrd])
    for name,it in zip(data_names[rwrd], np.arange(0, size_data[rwrd])):
      avg_conc[name] = {}
      for lab in labels:
        avg_conc[name][lab] = np.average(total_arr[name][lab])
        avg_conc_arr[rwrd][lab][it] =  avg_conc[name][lab]

# avg_conc should be divided by rhod?

    print avg_conc
    print avg_conc_arr[rwrd]

    for lab in labels:
      plt.plot(bin_centers[rwrd] * 1e6, avg_conc_arr[rwrd][lab] / bin_width[rwrd], label=rwrd + '_' + lab, linewidth=6)

#    data = total_arr["rain_rw_mom3"].values()

  plt.xlabel('radius [um]')
  plt.ylabel('PDF of concentration')
  plt.xscale('log')
  plt.legend()
  plt.yscale('log')
  plt.savefig('size_spectra_' + lvl + '_' + str(time_start) + '_' + str(time_end) +'.png')
  
