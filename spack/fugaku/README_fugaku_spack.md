The spack.yaml file contains definition of the packages required by UWLCM.
This definition is used to build a Spack environment with these packages.
Once it's built, UWLCM can be ran within this environment.

# Building Spack environment on Fugaku

1. set up private spack instance (based on: https://riken-rccs.github.io/fugaku-doc/docs/user-guide/sys-use/fugakuspackguide/build/en/intro.html):
- $ git clone https://github.com/RIKEN-RCCS/spack.git
- $ cd spack
- $ git checkout fugaku-v1.0.1
- in `/.spack/upstreams.yaml` add:
	upstreams:
	 spack-public-instance:
	    install_tree: /vol0004/apps/oss/spack/opt/spack
- $ spack repo add /vol0004/apps/oss/spack/var/spack/fugaku-packages/repos/spack_repo/fugaku/local
- $ spack repo add /vol0004/apps/oss/spack/var/spack/fugaku-packages/repos/spack_repo/fugaku/update
- $ spack repo add /vol0004/apps/oss/spack/var/spack/fugaku-packages/repos/spack_repo/fugaku/rist
- $ spack repo add /vol0004/apps/oss/spack/var/spack/fugaku-packages/repos/spack_repo/fugaku/rccs
- $ cp /vol0004/apps/oss/spack/etc/spack/packages.yaml ~/.spack/

2. build the environment from the spack.yaml file:
    2.1 run a job on compute node mounting /vol0004, e.g. an interactive job:
	 pjsub --interact --sparam wait-time=90 -L "rscunit=rscunit_ft01,rscgrp=int,node=1,elapse=6:00:00" -g hp250099 --no-check-directory -x PJM_LLIO_GFSCACHE=/vol0004
    2.2 enable private spack instance:
     $ . PRIVATE_SPACK_INSTANCE/spack/share/spack/setup-env.sh
    2.3 activate UWLCM spack environment:
	 $ spack env activate PATH_TO_FOLDER_WITH_THE_ENVIRONMENT
    2.4 generate build plan:
	  $ spack concretize --force
    2.5 build missing packages:
	  $ spack install

# How to run something in the environment:
- run a job on compute node mounting /vol0004, e.g. an interactive job (see 2.1)
- enable private spack instance (see 2.2)
- activate UWLCM spack environment (see 2.3)
- link gcc to gcc12 (optional, needed only for building libs and uwlcm, because it is done for  CMake try_compile to work):
	```
	$ GCC_PATH=$(spack location -i gcc@12.2.0)
		$ export LDFLAGS="-L${GCC_PATH}/lib64 -Wl,-rpath,${GCC_PATH}/lib64 -lstdc++"
	$ export CXXFLAGS="-Wl,-rpath,${GCC_PATH}/lib64"
	$ export CFLAGS="-Wl,-rpath,${GCC_PATH}/lib64"
	```


Within the environment, build libmpdata++, libcloudph++ and UWLCM using CMake, e.g.:
- limpdata++:
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_INSTALL_PREFIX=/vol0001/hp250099/u14261/builds_with_spack_mpi -DLIBMPDATAXX_MPI_THREAD_MULTIPLE=0

- libcloud:
cmake .. -DCMAKE_INSTALL_PREFIX=/vol0001/hp250099/u14261/builds_with_spack_mpi -DCMAKE_CXX_COMPILER=mpic++  -DLIBCLOUDPHXX_DISABLE_BINDINGS=1 -DCMAKE_BUILD_TYPE=Release

- UWLCM:
cmake .. -DCMAKE_INSTALL_PREFIX=/vol0001/hp250099/u14261/builds_with_spack_mpi -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_BUILD_TYPE=Release -DLIBMPDATAXX_MPI_THREAD_MULTIPLE=0 -DUWLCM_TIMING=1 -DUWLCM_DISABLE="PIGGYBACKER;2D_BLK_1M;2D_BLK_2M;3D_BLK_2M;2D_LGRNGN;2D_NONE;3D_NONE"
