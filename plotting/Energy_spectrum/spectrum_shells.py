# (C) Maciej Waruszewski

import h5py
import numpy as np
from sys import argv
import matplotlib.pyplot as plt

N = 65

h = h5py.File(argv[1], "r")

NN = N / 2 + 1

u = h["u"][0:N,0:N,0:N]
v = h["v"][0:N,0:N,0:N]
w = h["w"][0:N,0:N,0:N]


o = N * np.fft.fftfreq(N)

dk = 2.0
dkh = dk / 2

K = N * np.fft.rfftfreq(N)
print K
E = np.zeros(len(K))


uk = np.fft.fftn(u) / N ** 3
vk = np.fft.fftn(v) / N ** 3
wk = np.fft.fftn(w) / N ** 3

for i in xrange(N):
    for j in xrange(N):
        for k in xrange(N):
            q = np.sqrt(o[i] ** 2 + o[j] ** 2 + o[k] ** 2)
            for ix in xrange(len(K)):
                if abs(K[ix] - q) < 0.5:
                    e = abs(uk[i,j,k]) ** 2 + abs(vk[i,j,k]) ** 2 + abs(wk[i,j,k]) ** 2
                    E[ix] += 0.5 * e

K = K[2:]
E = E[2:]

plt.loglog(K, E, 'k-', lw = 2)
plt.savefig('spectrum_shells.svg')
