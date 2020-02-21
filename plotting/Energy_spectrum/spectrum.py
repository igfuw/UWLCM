# (C) Maciej Waruszewski

import h5py
import numpy as np
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


#velocities = ["u", "v", "w"]
#velocities = ["u", "v", "w", "cloud_rw_mom3", "rv", "th", "RH", "aerosol_rw_mom3"]
velocities = ["u", "cloud_rw_mom3", "w"]
#velocities = ["cloud_rw_mom3"]
#velocities = ["w"]

time_start = int(argv[1])
time_end = int(argv[2])
outfreq = int(argv[3])
from_lvl = int(argv[4])
to_lvl = int(argv[5])

directories = argv[6:len(argv):2]
labels = argv[7:len(argv):2]
print directories, labels

# read in nx, ny, nz
for directory, lab in zip(directories, labels):
  w3d = h5py.File(directory + "/timestep" + str(time_start).zfill(10) + ".h5", "r")["u"][:,:,:]
  nx, ny, nz = w3d.shape
  Exy_avg = {}
  for vel in velocities:
    Exy_avg[vel] = np.zeros(((nx+1)/2))
  
  for t in range(time_start, time_end+1, outfreq):
    filename = directory + "/timestep" + str(t).zfill(10) + ".h5"
    print filename
  
    for vel in velocities:
  
      w3d = h5py.File(filename, "r")[vel][:,:,:] # * 4. / 3. * 3.1416 * 1e3 
      
      for lvl in range(from_lvl, to_lvl+1):
        w2d = w3d[:, :, lvl]
        
        wkx = 1.0 / np.sqrt(nx - 1) * np.fft.rfft(w2d, axis = 0)
        wky = 1.0 / np.sqrt(ny - 1) * np.fft.rfft(w2d, axis = 1)
        
        #Ex = 0.5 * (np.abs(wkx) ** 2)
        Ex = (np.abs(wkx) ** 2)
        Ex = np.mean(Ex, axis = 1)
        #Ey = 0.5 * (np.abs(wky) ** 2)
        Ey = (np.abs(wky) ** 2)
        Ey = np.mean(Ey, axis = 0)
        
        Exy = 0.5 * (Ex + Ey)
        Exy_avg[vel] += Exy
  
      K = np.fft.rfftfreq(nx - 1)
      print nx
      print K
  #    plt.loglog(K, Exy)
#      lmbd = 50. / K # assume dx=50m
    
#    if (t == time_start and lab==labels[0]):
#      plt.loglog(lmbd, 2e-1* K**(-5./3.) )

  f1=plt.figure(1)  
  for vel in velocities:
    Exy_avg[vel] /= (time_end - time_start) / outfreq + 1
    Exy_avg[vel] /= to_lvl+1 - from_lvl
  #crudely scale
  #  Exy_avg[vel] /= Exy_avg[vel][len(Exy_avg[vel])-1]
  #  plt.loglog(lmbd, Exy_avg[vel] , linewidth=2, label=lab+"_"+vel)
    plt.loglog(K, Exy_avg[vel] , linewidth=2, label=lab+"_"+vel)
#plt.xlim(10**4,10**2)
plt.xlabel("l[m]")
plt.ylabel("PSD")
plt.legend()
plt.grid(True, which='both', linestyle='--')
plt.title("Mean PSD of w 322m<z<642m @3h")
f1.savefig('spectrum.png')
