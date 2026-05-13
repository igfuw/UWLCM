UWLCM requires libmpdata++ and libcloudph++ (both also developed at the Institute of Geophysics, Faculty of Physics, University of Warsaw), and several external dependencies. General workflow to build the model is to:

1. Obtain source code of libmpdata++, libcloudph++ and UWLCM from https://github.com/igfuw, e.g.:

```bash
$ git clone https://github.com/igfuw/libmpdataxx.git
```

2. Install external dependencies, three options:
	- manual installation (e.g. as in the Singularity image definition files found in UWLCM/singularity/ folder; uwlcm_ubuntu_24_04_cuda_12_9.def for non-MPI and uwlcm_ubuntu_24_04_cuda_12_9_mvapich2.def for MPI build)
	- download a pre-build singularity image from https://zenodo.org/records/15591478 and do all subsequent operations within that image, see [UWLCM/singularity/README_singularity.md](https://github.com/igfuw/UWLCM/blob/master/singularity/README_singularity.md) for details; Note that using this pre-build image will probably not work for MPI runs or will give poor performance
	- installation with Spack, see [UWLCM/spack/README_spack.md](https://github.com/igfuw/UWLCM/blob/master/spack/README_spack.md) for details

3. Build and install libmpdata++, e.g.:
- create “build” directory: libmpdataxx/libmpdata++/build
- in the build direcory issue:

```bash
$ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=gcc -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} .. ; make install
```

4. Build and install libcloudph++, e.g.:
- create “build” directory: libcloudphxx/build
- in the build direcory issue:

```bash
$ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=gcc -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DLIBCLOUDPHXX_FORCE_MULTI_CUDA=True .. ; make -j4 install
```

5. Build and install UWLCM using the Singularity image:
- create “build” directory: UWLCM/build
- in the build direcory issue:

```bash
$ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=gcc
-DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
-Dlibmpdata++_DIR=${INSTALL_DIR}/share/libmpdata++
-Dlibcloudph++_DIR=${INSTALL_DIR}/share/libcloudph++ .. ; make -j4 install
```

Note that when building for MPI runs, tell CMake to use the MPI compiler instead of gcc, e.g.:

```bash
-DCMAKE_CXX_COMPILER=mpic++
```
